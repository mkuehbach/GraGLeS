/*
	IMM_MicrostructureGenerator
	A program to instantiate two-staged Poisson-Voronoi tessellation microstructures of
	parent grains and their sub-grains with adjustable properties such as orientation, and dislocation density
	Copyright (C) 2016
	Christian Miessen (data structures), Markus Kühbach (physical metallurgy functionalities, PRNGs),
	Nikola Velinov (data structures), Luis Antonio Barrales-Mora (PRNGs, Math), Jonathan Nierengarten

	The work was funded by the DFG Reinhart-Koselleck project GO 335/44-1

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "microStructureHdl.h"
#include "Grains.h"
#include "../../../thirdparty/mandatory/voroxx/voro/src/voro++.hh"
#include "utilities.h"
#include <stdexcept>
#include "IterativeGrainScheduler.h"
#include <fstream>
#include <iostream>
#include <string>
#include <omp.h>
#include <numa.h>
#include "mymath.h"
#include "RTree.h"
#include "SubGrain.h"
#include "random.h"
#include <iomanip>
#include "lodepng.h"
//#include "hdf5.h"


using namespace voro;
using namespace std;
using namespace Eigen;

//Definition of the needed functions **********************************************************************************************************************************************************

timeLogger::timeLogger()
{
 entries = 0;
}

timeLogger::~timeLogger()
{
	titles.clear();
	times.clear();
}

void timeLogger::logev(const string title, double time)
{
	titles.push_back(title);
	times.push_back(time);
	entries++;
}

unsigned int timeLogger::get_entries( void )
{
	return entries;
}


microStructureHdl::microStructureHdl()
{
	healthy = true;
	m_container = NULL;
	m_grainScheduler = NULL;
	m_OrientationSpace = NULL;
	m_ParticlesPositions = NULL;
	m_seqRND = NULL;

	//number of gridpoints via mean radius argument
	if (Settings::Dimensionality == E_3D) {
		double voxel_subgrain = (4./3. * PI) * CUBE(Settings::NumberOfPointsPerSubGrain);
		double target = 0.5 + pow((Settings::NumberOfGrains * Settings::NumberOfSubgrains * voxel_subgrain), (1./3.));
		Settings::NumberOfGridpoints = ((unsigned int) target < MINIMUM_VE_SIZE) ? MINIMUM_VE_SIZE : (unsigned int) target;
		cout << "3D, NumberOfGridpoints " << Settings::NumberOfGridpoints << "\n";

		m_container = new DimensionalBuffer<unsigned int>(0, 0, 0,
			Settings::NumberOfGridpoints,
			Settings::NumberOfGridpoints,
			Settings::NumberOfGridpoints);
	}
	else { //E_2D
		double pixel_subgrain = PI / 2. * SQR(Settings::NumberOfPointsPerSubGrain);
		double target = 0.5 + pow((Settings::NumberOfGrains * Settings::NumberOfSubgrains * pixel_subgrain), (1./2.));
		Settings::NumberOfGridpoints = ((unsigned int) target < MINIMUM_VE_SIZE) ? MINIMUM_VE_SIZE : (unsigned int) target;
		cout << "2D, NumberOfGridpoints " << Settings::NumberOfGridpoints << "\n";

		m_container = new DimensionalBuffer<unsigned int>(0, 0, 0,
				Settings::NumberOfGridpoints,
				Settings::NumberOfGridpoints,
				1);
	}
	if ( Settings::NumberOfGridpoints == 0 ) {
		cerr << "NumberOfGridpoints must not be 0 !" << "\n";
		return;
	}

	m_h = 1. / Settings::NumberOfGridpoints;
	cout << "m_h " << m_h << "\n";
	first_id = 0 + 1;

	m_seqRND = new randomClass( -3000 );
	//vector<randomClass*>* m_threadlocalRND = new vector<randomClass*> [omp_get_max_threads() - 1];
}


microStructureHdl::~microStructureHdl() {
	delete m_container;
	delete m_grainScheduler;
	delete m_OrientationSpace;
	delete m_ParticlesPositions;
	delete m_seqRND;
}


void microStructureHdl::GeneratePolycrystallineStructureOfGrains()
{
	cout << __func__ << "\n";
	double tic = omp_get_wtime();

	switch (Settings::MicroGenMode) {
		case E_VORONOI: {
			GeneratorVoro();
			break;
		}
		default: {
			break;
		}
	}
	CopyContainer();

	myprofiler.logev( "GenerateGrainHull", (omp_get_wtime() - tic) );
}


void microStructureHdl::GenerateSubgrainStructureInEachGrain()
{
	cout << "GenerateSubgrainStructureInEachGrain..." << "\n";
	double tic = omp_get_wtime();

	ConstructSubgrains();

	myprofiler.logev( "GenerateSubgrainStructureInEachGrain", (omp_get_wtime() - tic) );
}


void microStructureHdl::ReadDiscreteOrientationSet() {
	cout << __func__ << "\n";
	double tic = omp_get_wtime();

	FILE * OriFromFile = NULL;
	OriFromFile = fopen(Settings::ReadFromFilename.c_str(), "r");
	if ( OriFromFile == NULL ) {
		cerr << "ParticleFile not existent or unreadable!" << "\n"; return;
	}

	int id, N = 0;
	char c;
	// count number of orientations
	do {
		c = fgetc(OriFromFile);
		if (c == '\n') {
			N++;
		}
	} while (c != EOF);
	rewind(OriFromFile);

	// read over header
	double vol, x, y, z, boxSize;
	double q[4];
	m_OrientationSpace = new vector<myQuaternion>;
	m_ParticlesPositions = new vector<Vector3d>;
	for (int i = 0; i < N-1; i++) {
		if ( i > 1 ) {
			fscanf(OriFromFile,
					"%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",
					&id, &x, &y, &z, &vol, &q[0], &q[1], &q[2], &q[3]);
			m_OrientationSpace->push_back(myQuaternion(q[0], q[1], q[2], q[3]));
			m_ParticlesPositions->push_back(
					Vector3d(x / boxSize, y / boxSize, z / boxSize));
		}
		else if ( i == 1 ) {
			do {
				c = fgetc(OriFromFile);
			} while (c != '\n');
			continue;
		}
		else if ( i == 0 ) {
			fscanf(OriFromFile, "%lf \n", &boxSize);
			do {
				c = fgetc(OriFromFile);
			} while (c != '\n');
			continue;
		}
		else {}
	}
	fclose(OriFromFile);
	
	cout << "ParticleFile read successfully with m_OrientationSpace.size "  << m_OrientationSpace->size() << "\n";
	//for ( unsigned int j = 0; j < m_OrientationSpace->size(); j++ ) cout << (*m_OrientationSpace)[j].get_q0() << ";" <<  (*m_OrientationSpace)[j].get_q1() << endl;
}


void microStructureHdl::GeneratorVoro() {
	cout << __func__ << "\n";

	bool is_periodic = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet

	if ( Settings::NumberOfGrains == 1 && Settings::VoronoiPeriodic == true ) {
		is_periodic = true;
	}

	voronoicell_neighbor c; //Voronoi Zelle die ihre Nachbaren kennt
//Erstellung eines Containers mit gewissen Dimensionen und Randbedingungen
	int blocks = (int) (pow((Settings::NumberOfGrains / 8), (1 / 3.)) + 1);
	if (blocks < 1) {
		blocks = 1;
	}
	voro::container con(0., 1., 0., 1., 0., 1.,
		blocks, blocks, blocks,
		is_periodic, is_periodic, is_periodic,
		8);

	c_loop_all vl(con);

	/**********************************************************/
	if ( Settings::GrainAggregation == E_POLYCRYSTAL ) {
		//JitterGrid process
		if ( Settings::GrainShape == E_FLAT ) {
			int Nx, Ny, Nz;
			double scaler = pow( (((double) Settings::NumberOfGrains) * Settings::DefGrainRelDimensionY ) , (1./2.) );
			if ( Settings::Dimensionality == E_3D ) {
				scaler = pow( (((double) Settings::NumberOfGrains)
					* Settings::DefGrainRelDimensionY * Settings::DefGrainRelDimensionZ ) , (1./3.) );
			}

			Nx = scaler;
			Ny = scaler / Settings::DefGrainRelDimensionY;
			Nz = 1;
			if ( Settings::Dimensionality == E_3D ) {
				Nz = scaler / Settings::DefGrainRelDimensionZ;
			}

			double dx = 1. / scaler;
			double dy = 1. / ( scaler * (Settings::DefGrainRelDimensionX / Settings::DefGrainRelDimensionY) );
			double dz = 1. / ( scaler * (Settings::DefGrainRelDimensionX / Settings::DefGrainRelDimensionZ) );
			cout << "JitterGrid FlatGrainMorphology Nx/Ny/Nz = " << Nx << ";" << Ny << ";" << Nz << "---" << dx << ";" << dy << ";" << dz << endl;

			int N = 0; //randomly distored point grid on unit cube
			if ( Settings::Dimensionality == E_2D ) {
				double z = 0.;
				for( double j = (dy / 2); j < 1.; j += dy ) {
					for( double i = (dx / 2); i < 1.; i += dx ) { // 0.5 * dx * m_seqRND->MersenneTwister(), shear band structures
						double x = dx * Settings::RandomnessX * (2*m_seqRND->MersenneTwister() - 1.) + i;
						double y = dy * Settings::RandomnessY * (2*m_seqRND->MersenneTwister() - 1.) + j;
						con.put( N, x, y, z);
						N++;
						cout << N << ";" << x << ";" << y << ";" << z << endl;
					}
				}
			}
			else if ( Settings::Dimensionality == E_3D ) {
				for ( double k = (dz / 2); k < 1.; k += dz ) {
					for( double j = (dy / 2); j < 1.; j += dy ) {
						for( double i = (dx / 2); i < 1.; i += dx ) { // 0.5 * dx * m_seqRND->MersenneTwister(), shear band structures
							double x = dx * Settings::RandomnessX * (2*m_seqRND->MersenneTwister() - 1.) + i;
							double y = dy * Settings::RandomnessY * (2*m_seqRND->MersenneTwister() - 1.) + j;
							double z = dz * Settings::RandomnessZ * (2*m_seqRND->MersenneTwister() - 1.) + k;
							con.put( N, x, y, z);
							N++;
							cout << N << ";" << x << ";" << y << ";" << z << endl;
						}
					}
				}
			}

			Settings::NumberOfGrains = N;
			cout << "JitterGrid-FlatGrainMorphology - I am resetting the number of grains to N/NumberOfGrains = " << N << ";" << Settings::NumberOfGrains << endl;
		}
		else { //E_GLOBULITIC
			// Randomly add particles into the container
			for (int i = 0; i < Settings::NumberOfGrains; i++) {
				double x = m_seqRND->MersenneTwister();
				double y = m_seqRND->MersenneTwister();
				double z = 0.;
				if (Settings::Dimensionality == E_3D) {
					z = m_seqRND->MersenneTwister();
				}
				con.put(i, x, y, z); //assign each Voronoi cell i its coordinates x,y,z
			}
		}
	}
	else if ( Settings::GrainAggregation == E_BICRYSTAL ) {
		//MK::model utilized for M. K\"uhbach, F. Roters, et. al. nonlocal solving of W.B. Hutchinson's old GB nucleation bicrystal experiments on iron
		//the spectral solver applies always periodic boundary conditions so strictly speaking the lower parts see
		//x || RD, y || TD, z || ND
		int i = 0;
		//place upper container for the bicrystal
		double x = 0.5;
		double y = 0.5;
		double z = 0.75;
		con.put(i, x, y, z);

		//place lower container for the bicrystal
		x = 0.5;
		y = 0.5;
		z = 0.25;
		i++;
		con.put(i, x, y, z);
	}
	else {
		cerr << "E_SINGLECRYSTAL currently not implemented!" << "\n";
	}

	//MK::from now on the number of grains is required to be fixed!

	//check Voronoi initialization
	m_grains.resize(Settings::NumberOfGrains + 1);

	/**********************************************************/

	vector<vector<Eigen::Vector3d> > initialHulls;
	vector<double> grainVolume;
	vector<double> cellCoordinates;
	if (vl.start()) {
		initialHulls.resize(Settings::NumberOfGrains + 1);
		grainVolume.resize(Settings::NumberOfGrains + 1);
		do {
			double cur_x, cur_y, cur_z;
			con.compute_cell(c, vl);
			//new: get the grain_id
			unsigned int box_id = vl.pid() + 1;
			vl.pos(cur_x, cur_y, cur_z);
			grainVolume[box_id] = c.volume();
			c.vertices(cur_x, cur_y, cur_z, cellCoordinates); //Erstellung des Polyheders
			for (unsigned int i = 0; i < cellCoordinates.size() / 3; i++) {
				initialHulls.at(box_id).push_back(
						Eigen::Vector3d(cellCoordinates.at(3 * i + 1),
								cellCoordinates.at(3 * i),
								cellCoordinates.at(3 * i + 2)));
			}
		} while (vl.inc());

		double x, y, z, rx, ry, rz;
		int cell_id;
		int z_limit = Settings::NumberOfGridpoints;
		if (Settings::Dimensionality == E_2D) {
			z_limit = 1;
		}
		for (int k = 0; k < z_limit; k++) {
			for (int j = 0; j < Settings::NumberOfGridpoints; j++) {
				for (int i = 0; i < Settings::NumberOfGridpoints; i++) {
					x = i * m_h;
					y = j * m_h;
					z = k * m_h;
					if (con.find_voronoi_cell(x, y, z, rx, ry, rz, cell_id)) {
						unsigned int box_id = cell_id + 1;
						m_container->setValueAt(j, i, k, box_id);
					}
					else {
						m_container->setValueAt(j, i, k, 0);
					}
				}
			}
		}
	} 
	else {
		throw runtime_error("Voronoy container error at start() method!");
	}

	//generate grainboxes by evaluating the initial hulls objects
	InitializeGrains(initialHulls, grainVolume);

	cout << "Parent grain geometry constructed " << endl;
}


void microStructureHdl::ConstructSubgrains()
{
	cout << __func__ << "\n";
	#pragma omp parallel
	{
		randomClass* threadlocalRNG = new randomClass();
		//threadlocalRNG->initPM(-1 * (omp_get_thread_num()) - 1);
		uint32_t overflowint = (uint32_t) pow(2.0, 31);
		threadlocalRNG->initMT( overflowint - omp_get_thread_num() - 1 );
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfGrains) {
				m_grains[id]->SubGrainConstructor(*threadlocalRNG);
			}
		}
		delete threadlocalRNG;
	}
}


void microStructureHdl::CopyContainer()
{
	for (auto it : m_grains) {
		if (it != NULL) {
			it->copyContainerToGrain(m_container);
		}
	}
}


void microStructureHdl::FindNeighbors()
{
	RTree<unsigned int, int, 3, float> tree;
	int min[3], max[3];
	for (unsigned int i = 1; i <= Settings::NumberOfGrains; i++) {
		if (m_grains[i] == NULL)
			continue;
		min[0] = m_grains[i]->getMinX();
		min[1] = m_grains[i]->getMinY();
		min[2] = m_grains[i]->getMinZ();
		max[0] = m_grains[i]->getMaxX();
		max[1] = m_grains[i]->getMaxY();
		max[2] = m_grains[i]->getMaxZ();
		tree.Insert(min, max, i);
	}
	for (unsigned int id = 1; id <= Settings::NumberOfGrains; id++) {
		m_grains[id]->computeDirectNeighbours(tree);
	}
}


void microStructureHdl::DistributeGrainOriAndSee()
{
	cout << __func__ << "\n";
	double tic = omp_get_wtime();

	//TODO::readPreferenceOrientationFromFile();
	for (vector<Grains*>::iterator it = ++m_grains.begin(); it != m_grains.end(); it++) {
		myQuaternion ori;
		if (Settings::MicroGenMode == E_VORONOI) {
			if (Settings::TextureSampling == E_PICK_RANDOMLY) {
				unsigned int mOrientations = m_OrientationSpace->size();
				unsigned int randomOri = m_OrientationSpace->size();
				while ( randomOri >= mOrientations ) {
					randomOri = m_seqRND->MersenneTwister() * m_OrientationSpace->size();
				}
				ori = (*m_OrientationSpace)[randomOri]; 								//pick randomly from list of predefined orientations
				cout << "Grain " << (*it)->get_ID() << " becomes mapped on randomly picked UserDefinedOrientation " << randomOri << endl;
				//ori.randomOriShoemakeQuat(*m_seqRND);
			}
			else { //##MK::at the moment equivalent to == E_DEFAULT_SAMPLING
				unsigned int mOrientations = m_OrientationSpace->size();
				unsigned int linearOri = (*it)->get_ID() % mOrientations;
				ori = (*m_OrientationSpace)[linearOri];
				cout << "Grain " << (*it)->get_ID() << " becomes mapped on linearly identified UserDefinedOrientation " << linearOri << endl;
			}
		}
		else {
			ori.randomOriShoemakeQuat(*m_seqRND); //random orientations
		}

		//assign grain the orientation and scattering properties
		(*it)->set_Orientation(ori);
		(*it)->set_PreforiProperties( FindNextPreferenceOrientation(ori) );
		double gr_see_mu = (*it)->get_SEEFromPrefori();
		double gr_see_sig = (*it)->get_SEEGrainScatterFromPrefOri();
		//cout << "Grain = " << (*it)->get_ID() << " prefori mu/sig/orisig = " << gr_see_mu << ";" << gr_see_sig << ";" << (*it)->get_OriScatterFromPrefOri() << "----prefori q0 " <<  (*it)->get_PrefOriQuatQ0() << ";" << (*it)->get_PrefOriQuatQ1() << ";" << (*it)->get_PrefOriQuatQ2() << ";" << (*it)->get_PrefOriQuatQ3() << endl;
		(*it)->set_SEE( m_seqRND->r4_nor(gr_see_mu, gr_see_sig) ); //set stored elastic energy of grain to a specific value from a normal distribution
	}
	cout << "Sample sub-grain orientation from reference orientation of grain" << "\n";

	myprofiler.logev("DistriGrainPrefOris", (omp_get_wtime() - tic));
}


void microStructureHdl::DistributeSubgrainOrientations()
{
	cout << __func__ << "\n";
	double tic = omp_get_wtime();

	#pragma omp parallel
	{
		randomClass* threadlocalRNG = new randomClass();
		//threadlocalRNG->initPM(-1 * (omp_get_thread_num()) - 1);
		uint32_t overflowint = (uint32_t) pow(2.0, 31);
		threadlocalRNG->initMT( overflowint - omp_get_thread_num() - 1 );

		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfGrains) {
				m_grains[id]->generateSubGrainOri(*threadlocalRNG);
			}
		}
		delete threadlocalRNG;
	}

	myprofiler.logev("DistrSubgrainOris", (omp_get_wtime() - tic) );
}


void microStructureHdl::ReadAdditionalInputFiles() {
	double tic = omp_get_wtime();

	ReadDiscreteOrientationSet();
	ReadPreferenceOrientationFromFile();

	myprofiler.logev("ReadAdditionalInput", (omp_get_wtime() - tic));
}


void microStructureHdl::ReportConfig()
{
	cout << __func__ << "\n";
	cout << "SimulationId " << Settings::SimulationId << "\n";
	cout << "NumberOfGrains " << Settings::NumberOfGrains << "\n";
	cout << "NumberOfSubgrains " << Settings::NumberOfSubgrains << "\n";
	cout << "NumberOfPointsPerSubGrain " << Settings::NumberOfPointsPerSubGrain << "\n";
	cout << "PhysEdgeLengthMatPoint " << Settings::PhysEdgeLenMatPoint << "\n";
	cout << "NumberOfGridpoints " << Settings::NumberOfGridpoints << "\n";
	cout << "MaximumNumberOfThreads " << Settings::MaximumNumberOfThreads << "\n";
	cout << "SubgrainOriScatter " << Settings::SubgrainOriScatter << "\n";
	cout << "StoredElasticEnergyMin " << Settings::StoredElasticEnergyMin << "\n";
	cout << "StoredElasticEnergyMax " << Settings::StoredElasticEnergyMax << "\n";
	cout << "StoredElasticScatterGrain " << Settings::StoredElasticScatterGrain << "\n";
	cout << "StoredElasticScatterSubgrain " << Settings::StoredElasticScatterSubgrain << "\n";
	cout << "DefGrainRelDimensionX " << Settings::DefGrainRelDimensionX << "\n";
	cout << "DefGrainRelDimensionY " << Settings::DefGrainRelDimensionY << "\n";
	cout << "DefGrainRelDimensionZ " << Settings::DefGrainRelDimensionZ << "\n";
	cout << "RandomnessX " << Settings::RandomnessX << "\n";
	cout << "RandomnessY " << Settings::RandomnessY << "\n";
	cout << "RandomnessZ " << Settings::RandomnessZ << "\n";
	cout << "VoronoiPeriodic " << Settings::VoronoiPeriodic << "\n";
	cout << "ExecuteInParallel " << Settings::ExecuteInParallel << "\n";
	cout << "StatusHealthy " << Settings::StatusHealthy << "\n";
	cout << "BreakPerX " << Settings::BreakPerX << "\n";
	cout << "BreakPerY " << Settings::BreakPerY << "\n";
	cout << "BreakPerZ " << Settings::BreakPerZ << "\n";
	cout << "MicroGenMode " << Settings::MicroGenMode << "\n";
	cout << "CrystalStructure " << Settings::CrystalStructure << "\n";
	cout << "TextureSampling " << Settings::TextureSampling << "\n";
	cout << "TextureGen " << Settings::TextureGen << "\n";
	cout << "Dimensionality " << Settings::Dimensionality << "\n";
	cout << "GrainAggregation " << Settings::GrainAggregation << "\n";
	cout << "GrainShape " << Settings::GrainShape << "\n";
	cout << "ConfigFileName " << Settings::ConfigFileName << "\n";
	cout << "ReadFromFilename " <<  Settings::ReadFromFilename << "\n";
	cout << "AdditionalFilename " << Settings::AdditionalFilename << "\n";
	cout << "ResultsFileName " << Settings::ResultsFileName << "\n";
}


void microStructureHdl::ReadPreferenceOrientationFromFile() {
	FILE* file = NULL;
	file = fopen(Settings::AdditionalFilename.c_str(), "r");
	if ( file == NULL ) {
		std::cout << "AdditionalFile not existent or unreadable!" << std::endl;
		Settings::StatusHealthy = false;
		return;
	}
	int N = 0;
	char c;
// count number of orientations
	do {
		c = fgetc(file);
		if (c == '\n')
			N++;
	} while (c != EOF);
	rewind(file);
	m_PreferenceOrientations = new vector<myPreferenceOri>;
	for (int i = 0; i < N; i++) {
		double phi1, PHI, phi2, SEE, oriScatter, SEEGrainScatter, SEESubgrainScatter, RelSubgrainSizeScaler;
		fscanf(file, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf ", &phi1, &PHI,
				&phi2, &oriScatter, &SEE, &SEEGrainScatter, &SEESubgrainScatter, &RelSubgrainSizeScaler ); //SEE reads as mu and SEEScatter as sigma
		char c;
		//read comments
		do {
			c = fgetc(file);
		} while (c != '\n');
		m_PreferenceOrientations->push_back(
				myPreferenceOri(DEG2RAD(phi1), DEG2RAD(PHI), DEG2RAD(phi2), DEG2RAD(oriScatter),
				SEE, SEEGrainScatter, SEESubgrainScatter, RelSubgrainSizeScaler)); //define a preference orientation with specific instead of default properties in terms of mu, sigma and oriscatter
		//cout << "PrefRefs = " << phi1 << ";" << PHI << ";" << phi2 << ";" << oriScatter << ";" << SEE << ";" << SEEGrainScatter << ";" << SEESubgrainScatter << ";" << RelSubgrainSizeScaler << endl;
	}
	fclose(file);
	//cout << "PrefRefDefault = " << Settings::SubgrainOriScatter << ";" << (Settings::StoredElasticEnergyMax+Settings::StoredElasticEnergyMin)*0.5 << ";" << Settings::StoredElasticScatterGrain << ";" << 1. << endl;
}


myPreferenceOri microStructureHdl::FindNextPreferenceOrientation(
		myQuaternion ori) {
	if (Settings::TextureGen == E_USE_PREFERENCEORI) {
		double min = DEG2RAD(15.);
		unsigned int idx = 0;
		for ( unsigned int p = 0; p < m_PreferenceOrientations->size(); p++ ) {
			double misori = ori.misorientationCubicQxQ( &((*m_PreferenceOrientations)[p].ori) );
			if ( misori <= min ) {
				min = misori;
				idx = p;
			}
		}
		if (min <= DEG2RAD(10.)) {
			//return *i;
//cout << "Finding preference orientation with min/idx = " << min << ";" << idx << endl;
			return myPreferenceOri( ori, (*m_PreferenceOrientations)[idx].subgrainsScatterOri, (*m_PreferenceOrientations)[idx].SEE, (*m_PreferenceOrientations)[idx].SEEGrainScatter, (*m_PreferenceOrientations)[idx].SEESubgrainScatter, (*m_PreferenceOrientations)[idx].RelSubgrainSizeScaling );

		}
		else {
//cout << "Check for preference orientation but too far off = " << min << endl;
			return myPreferenceOri(ori); //in case no close enough preference orientation was found, thus the ori itself is returned as a PreferenceOrientation
		}
	}
	else {
//cout << "Simply returning preference orientation " << endl;
		return myPreferenceOri(ori);
	}
}

void microStructureHdl::DistributeSubgrainSee() {
	double gtime = omp_get_wtime();
	cout << "Sample Stored Elastic Energy" << endl;
#pragma omp parallel
	{
		randomClass* threadlocalRNG = new randomClass();
		//threadlocalRNG->initPM(-1 * (omp_get_thread_num()) - 1);
		uint32_t overflowint = (uint32_t) pow(2.0, 31);
		threadlocalRNG->initSHR3( overflowint - omp_get_thread_num() - 1 );
		threadlocalRNG->initR4Uni( overflowint - omp_get_thread_num() - 1 );
		threadlocalRNG->r4_nor_setup();

		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfGrains) {
				m_grains[id]->generateSubGrainSEE(
						*threadlocalRNG);
			}
		}
		delete threadlocalRNG;
	}

	myprofiler.logev("DistrStoredElasticEnergy", (omp_get_wtime() - gtime) );
}


void microStructureHdl::RehashGrainIds() {
	double timer = omp_get_wtime();
	cout << ">Rehashing grain IDs" << endl;

	//classical (Mie\ss{}en and K\"uhbach 2015,2016,2017 grain IDs start with 1
	int offset = 0;

	//MK::consider a potential breaking of the periodicity by adding pair(s) of "air" grains 1,2,3,4,5,6
	if ( Settings::BreakPerX == true )	offset += 2;
	if ( Settings::BreakPerY == true )	offset += 2;
	if ( Settings::BreakPerZ == true )	offset += 2;

	int newoffset = 0;
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {
			newoffset = (*itG)->copySubgrainsToGlobalContainer(m_container, offset); // implicit rehashing of all ID's
			offset = --newoffset;
		}
		else {
			newoffset = (*itG)->copySubgrainsToGlobalContainer(m_container, offset);
			offset = newoffset;
		}
	}

	myprofiler.logev("RehashingIDs", (omp_get_wtime() - timer) );
}


void microStructureHdl::BreakPeriodicity() {
	cout << __func__ << "\n";
	double timer = omp_get_wtime();
	//add "air" grain pair at the boundary to break the symmetry of the domain
	unsigned int offset = 1;

	if ( Settings::BreakPerX == true ) {
		SetFirstId( true, 2 );
		offset += 2;
		unsigned int x1 = m_container->getMaxX() - 1; //getMaxI are exclusive...
		unsigned int x2 = m_container->getMaxX() - 2;
		//reset all m_container values to "air" grain 1 at (x=0,y,z)
		for ( unsigned int z = m_container->getMinZ(); z < m_container->getMaxZ(); ++z ) {
			for ( unsigned int y = m_container->getMinY(); y < m_container->getMaxY(); ++y ) {
				m_container->setValueAt( y, 0, z, offset-2 );
				m_container->setValueAt( y, 1, z, offset-2 );
				m_container->setValueAt( y, x2, z, offset-1 );
				m_container->setValueAt( y, x1, z, offset-1 );
			}
		}
	}
	if ( Settings::BreakPerY == true ) {
		SetFirstId( true, 2 );
		offset += 2;
		unsigned int y1 = m_container->getMaxY() - 1;
		unsigned int y2 = m_container->getMaxY() - 2;
		for ( unsigned int z = m_container->getMinZ(); z < m_container->getMaxZ(); ++z ) {
			for ( unsigned int x = m_container->getMinX(); x < m_container->getMaxX(); ++x ) {
				m_container->setValueAt( 0, x, z, offset-2 );
				m_container->setValueAt( 1, x, z, offset-2 );
				m_container->setValueAt( y2, x, z, offset-1 );
				m_container->setValueAt( y1, x, z, offset-1 );
			}
		}
	}
	if ( Settings::BreakPerZ == true ) {
		SetFirstId( true, 2 );
		offset += 2;
		unsigned int z1 = m_container->getMaxZ() - 1;
		unsigned int z2 = m_container->getMaxZ() - 2;
		for ( unsigned int y = m_container->getMinY(); y < m_container->getMaxY(); ++y ) {
			for ( unsigned int x = m_container->getMinX(); x < m_container->getMaxX(); ++x ) {
				m_container->setValueAt( y, x, 0, offset-2 );
				m_container->setValueAt( y, x, 1, offset-2 );
				m_container->setValueAt( y, x, z2, offset-1 );
				m_container->setValueAt( y, x, z1, offset-1 );
			}
		}
	}

	if ( Settings::BreakPerX == true || Settings::BreakPerY == true || Settings::BreakPerZ == true ) {
		//update subgrain volume
		vector<Grains*>::iterator itG;
		for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
			if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {
				vector<SubGrain*>::iterator it;
				for (it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++) {
					double oldVolume = (*it)->get_Volume();
					double newVolume = 0.;
					unsigned int cand = (*it)->get_ID();
					for ( int z = (*it)->getMinZ(); z < (*it)->getMaxZ(); z++ ) { //##MK::could be improved by just scanning changed boundary cells..
						for ( int y = (*it)->getMinY(); y < (*it)->getMaxY(); y++ ) {
							for ( int x = (*it)->getMinX(); x < (*it)->getMinX(); x++ ) {
								if ( m_container->getValueAt( y, x, z ) == cand )
									newVolume++;
							}
						}
					}
					//##implement setter!
					std::cout << "Sub-grain " << (*it)->get_ID() << " old volume " << oldVolume << " new " << newVolume << endl;
				}
			} else {
				double oldVolume = (*itG)->get_Volume();
				double newVolume = 0.;
				unsigned int cand = (*itG)->get_ID();
				for ( int z = (*itG)->getMinZ(); z < (*itG)->getMaxZ(); z++ ) { //##MK::could be improve by just scanning changed boundary cells..
					for ( int y = (*itG)->getMinY(); y < (*itG)->getMaxY(); y++ ) {
						for ( int x = (*itG)->getMinX(); x < (*itG)->getMinX(); x++ ) {
							if ( m_container->getValueAt( y, x, z ) == cand )
								newVolume++;
						}
					}
				}
				//##implement setter!
				std::cout << "Grain w/o subgrain " << (*itG)->get_ID() << " old volume " << oldVolume << " new " << newVolume << endl;
			}
		}
	}
	//nothing to do otherwise

	myprofiler.logev( "BreakPeriodicity", (omp_get_wtime() - timer) );
}


void microStructureHdl::SaveNeXus() {
	cout << __func__ << "\n";
	double tic = omp_get_wtime();

	HdfFiveSeqHdl h5w = HdfFiveSeqHdl( Settings::ResultsFileName );
	ioAttributes anno = ioAttributes();
	string grpnm = "";
	string dsnm = "";

	vector<double> f64;
	vector<float> f32;
	vector<unsigned int> u32;
	vector<int> i32;
	vector<unsigned char> u8;

	grpnm = "/entry1";
	anno = ioAttributes();
	anno.add( "NX_class", string("NXentry") );
	if ( h5w.nexus_write_group( grpnm, anno ) != MYHDF5_SUCCESS ) { return; }

	grpnm = "/entry1/ms";
	anno = ioAttributes();
	anno.add( "NX_class", string("NXms_snapshot") );
	if ( h5w.nexus_write_group( grpnm, anno ) != MYHDF5_SUCCESS ) { return; }

	/*
	grpnm = "/entry1/ms/grid";
	anno = ioAttributes();
	anno.add( "NX_class", string("NXcg_grid") );
	if ( h5w.nexus_write_group( grpnm, anno ) != MYHDF5_SUCCESS ) { return; }
	*/

	dsnm = grpnm + "/unknown_offset";
	unsigned long idoff = this->GetFirstId(); //for grains or cells?
	anno = ioAttributes();
	if ( h5w.nexus_write( dsnm, idoff, anno ) != MYHDF5_SUCCESS ) { return; }

	dsnm = grpnm + "/average_subgrain_discretization";
	u32 = vector<unsigned int>( Settings::Dimensionality, Settings::NumberOfPointsPerSubGrain);
	anno = ioAttributes();
	if ( h5w.nexus_write(
		dsnm,
		io_info({u32.size()}, {u32.size()}, MYHDF5_COMPRESSION_NONE, 0x00),
		u32, anno ) != MYHDF5_SUCCESS ) { return; }
	u32 = vector<unsigned int>();

	dsnm = grpnm + "/extent"; //real extent ?
	u32 = vector<unsigned int>( Settings::Dimensionality, Settings::NumberOfGridpoints);
	anno = ioAttributes();
	if ( h5w.nexus_write(
		dsnm,
		io_info({u32.size()}, {u32.size()}, MYHDF5_COMPRESSION_NONE, 0x00), 
		u32, anno ) != MYHDF5_SUCCESS ) { return; }
	u32 = vector<unsigned int>();

	dsnm = grpnm + "/number_of_subgrains";
	unsigned int n_subgr = CountNumberOfSubgrains();
	anno = ioAttributes();
	if ( h5w.nexus_write( dsnm, n_subgr, anno ) != MYHDF5_SUCCESS ) { return; }

	grpnm = "/entry1/ms";
	dsnm = grpnm + "/grain_is_subgrain";
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((*itG) != NULL) {
			if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
				for (vector<SubGrain*>::iterator itS = ++((*itG)->m_SubGrains.begin()); itS != (*itG)->m_SubGrains.end(); itS++) {
					if ((*itS) != NULL) {
						u8.push_back(0x01);
					}
				}
			}
			else {
				u8.push_back(0x00);
			}
		}
	}
	anno = ioAttributes();
	if ( h5w.nexus_write(
		dsnm,
		io_info({u8.size()}, {u8.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		u8, anno ) != MYHDF5_SUCCESS ) { return; }
	u8 = vector<unsigned char>();

	dsnm = grpnm + "/grain_identifier";
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((*itG) != NULL) {
			if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
				for (vector<SubGrain*>::iterator itS = ++((*itG)->m_SubGrains.begin()); itS != (*itG)->m_SubGrains.end(); itS++) {
					if ((*itS) != NULL) {
						u32.push_back((*itS)->get_ID());
					}
				}
			}
			else {
				u32.push_back((*itG)->get_ID());
			}
		}
	}
	anno = ioAttributes();
	if ( h5w.nexus_write(
		dsnm,
		io_info({u32.size()}, {u32.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		u32, anno ) != MYHDF5_SUCCESS ) { return; }
	u32 = vector<unsigned int>();

	dsnm = grpnm + "/grain_barycentre_naive";
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((*itG) != NULL) {
			if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
				for (vector<SubGrain*>::iterator itS = ++((*itG)->m_SubGrains.begin()); itS != (*itG)->m_SubGrains.end(); itS++) {
					if ((*itS) != NULL ) {
						f64.push_back((((*itS)->getMaxX() - (*itS)->getMinX()) / 2 + (*itS)->getMinX()));
						f64.push_back((((*itS)->getMaxY() - (*itS)->getMinY()) / 2 + (*itS)->getMinY()));
						if (Settings::Dimensionality == E_3D) {
							f64.push_back((((*itS)->getMaxZ() - (*itS)->getMinZ()) / 2 + (*itS)->getMinZ()));
						}
					}
				}
			}
			else {
				f64.push_back((((*itG)->getMaxX() - (*itG)->getMinX()) / 2 + (*itG)->getMinX()));
				f64.push_back((((*itG)->getMaxY() - (*itG)->getMinY()) / 2 + (*itG)->getMinY()));
				if (Settings::Dimensionality == E_3D) {
					f64.push_back((((*itG)->getMaxZ() - (*itG)->getMinZ()) / 2 + (*itG)->getMinZ()));
				}
			}
		}
	}
	anno = ioAttributes();
	anno.add( "unit", string("m ?????") );
	size_t n_cols = (Settings::Dimensionality == E_3D) ? 3 : 2;
	if ( h5w.nexus_write(
		dsnm,
		io_info({f64.size() / n_cols, n_cols}, {f64.size() / n_cols, n_cols}, MYHDF5_COMPRESSION_GZIP, 0x01),
		f64, anno ) != MYHDF5_SUCCESS ) { return; }
	f64 = vector<double>();

	dsnm = grpnm + "/grain_orientation";
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((*itG) != NULL) {
			if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
				for (vector<SubGrain*>::iterator itS = ++((*itG)->m_SubGrains.begin()); itS != (*itG)->m_SubGrains.end(); itS++) {
					if ((*itS) != NULL) {
						f64.push_back((*itS)->get_Orientation()->get_q0());
						f64.push_back((*itS)->get_Orientation()->get_q1());
						f64.push_back((*itS)->get_Orientation()->get_q2());
						f64.push_back((*itS)->get_Orientation()->get_q3());
					}
				}
			}
			else {
				f64.push_back((*itG)->getOri().get_q0());
				f64.push_back((*itG)->getOri().get_q1());
				f64.push_back((*itG)->getOri().get_q2());
				f64.push_back((*itG)->getOri().get_q3());
			}
		}
	}
	anno = ioAttributes();
	if ( h5w.nexus_write(
		dsnm,
		io_info({f64.size() / 4, 4}, {f64.size() / 4, 4}, MYHDF5_COMPRESSION_GZIP, 0x01),
		f64, anno ) != MYHDF5_SUCCESS ) { return; }
	f64 = vector<double>();

	dsnm = grpnm + "/grain_aabb";
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((*itG) != NULL) {
			if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
				for (vector<SubGrain*>::iterator itS = ++((*itG)->m_SubGrains.begin()); itS != (*itG)->m_SubGrains.end(); itS++) {
					if ((*itS) != NULL) {
						i32.push_back((*itS)->getMinX());
						i32.push_back((*itS)->getMaxX());
						i32.push_back((*itS)->getMinY());
						i32.push_back((*itS)->getMaxY());
						if (Settings::Dimensionality == E_3D) {
							i32.push_back((*itS)->getMinZ());
							i32.push_back((*itS)->getMaxZ());
						}
					}
				}
			}
			else {
				i32.push_back((*itG)->getMinX());
				i32.push_back((*itG)->getMaxX());
				i32.push_back((*itG)->getMinY());
				i32.push_back((*itG)->getMaxY());
				if (Settings::Dimensionality == E_3D) {
					i32.push_back((*itG)->getMinZ());
					i32.push_back((*itG)->getMaxZ());
				}
			}
		}
	}
	anno = ioAttributes();
	n_cols = (Settings::Dimensionality == E_3D) ? 6 : 4;
	if ( h5w.nexus_write(
		dsnm,
		io_info({i32.size() / n_cols, n_cols}, {i32.size() / n_cols, n_cols}, MYHDF5_COMPRESSION_GZIP, 0x01),
		i32, anno ) != MYHDF5_SUCCESS ) { return; }
	i32 = vector<int>();

	dsnm = grpnm + "/grain_size";
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((*itG) != NULL) {
			if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
				for (vector<SubGrain*>::iterator itS = ++((*itG)->m_SubGrains.begin()); itS != (*itG)->m_SubGrains.end(); itS++) {
					if ((*itS) != NULL) {
						double vol = (*itS)->get_Volume();
						double divisor = (Settings::Dimensionality == E_3D ) ? (m_h * m_h * m_h) : (m_h * m_h);
						if ( divisor > DBL_EPSILON ) {
							f64.push_back(vol / divisor);
						}
						else {
							cerr << "Division by zero!" << "\n";
						}
					}
				}
			}
			else {
				double vol = (*itG)->get_Volume();
				double divisor = (Settings::Dimensionality == E_3D) ? (m_h * m_h * m_h) : (m_h * m_h);
				if ( divisor > DBL_EPSILON ) {
					f64.push_back( vol / divisor );
				}
				else {
					cerr << "Division by zero!" << "\n";
				}
			}
		}
	}
	anno = ioAttributes();
	/*
	if (Settings::Dimensionality == E_3D ) {
		anno.add( "unit", string("m^3 ?????") );
	}
	else {
		anno.add( "unit", string("m^2 ?????") );
	}
	*/
	if ( h5w.nexus_write(
		dsnm,
		io_info({f64.size()}, {f64.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		f64, anno ) != MYHDF5_SUCCESS ) { return; }
	f64 = vector<double>();
	
	dsnm = grpnm + "/stored_elastic_energy";
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ( (*itG) != NULL ) {
			if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
				for (vector<SubGrain*>::iterator itS = ++((*itG)->m_SubGrains.begin()); itS != (*itG)->m_SubGrains.end(); itS++) {
					if ((*itS) != NULL ) {
						f64.push_back((*itS)->get_SEE());
					}
				}
			}
			else {
				f64.push_back((*itG)->getSEE());
			}
		}
	}
	anno = ioAttributes();
	anno.add( "unit", string("1/m^2 ?????") );
	if ( h5w.nexus_write( 
		dsnm,
		io_info({u32.size()}, {u32.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		f64, anno ) != MYHDF5_SUCCESS ) { return; }
	f64 = vector<double>();

	dsnm = grpnm + "/grain_parent_identifier";
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ( (*itG) != NULL ) {
			if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
				for (vector<SubGrain*>::iterator itS = ++((*itG)->m_SubGrains.begin()); itS != (*itG)->m_SubGrains.end(); itS++) {
					if ((*itS) != NULL ) {
						u32.push_back((*itG)->get_oldID());
					}
				}
			}
			else {
				u32.push_back((*itG)->get_oldID());
			}
		}
	}
	anno = ioAttributes();
	if ( h5w.nexus_write( 
		dsnm,
		io_info({u32.size()}, {u32.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		u32, anno ) != MYHDF5_SUCCESS ) { return; }
	u32 = vector<unsigned int>();

	dsnm = grpnm + "/grain_disorientation_to_parent";
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ( (*itG) != NULL ) {
			myQuaternion parentori = (*itG)->getOri();
			if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
				for (vector<SubGrain*>::iterator itS = ++((*itG)->m_SubGrains.begin()); itS != (*itG)->m_SubGrains.end(); itS++) {
					if ((*itS) != NULL ) {
						myQuaternion* ori = (*itS)->get_Orientation();
						f64.push_back(parentori.misorientationCubicQxQ(ori));
					}
				}
			}
			else {
				f64.push_back(0.);
			}
		}
	}
	anno = ioAttributes();
	anno.add( "unit", string("° ?????") );
	if ( h5w.nexus_write( 
		dsnm,
		io_info({f64.size()}, {f64.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		f64, anno ) != MYHDF5_SUCCESS ) { return; }
	f64 = vector<double>();

	/*
	dsnm = grpnm + "/grain_target_disorientation_scatter"; //?????
	for (vector<Grains*>::iterator itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ( (*itG) != NULL ) {
			double oriscatter = (*itG)->get_OriScatterFromPrefOri();
			if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
				for (vector<SubGrain*>::iterator itS = ++((*itG)->m_SubGrains.begin()); itS != (*itG)->m_SubGrains.end(); itS++) {
					if ((*itS) != NULL ) {
						myQuaternion* ori = (*itS)->get_Orientation();
						f64.push_back(oriscatter);
					}
				}
			}
			else {
				f64.push_back((*itG)->get_OriScatterFromPrefOri());
			}
		}
	}
	//###
	*/

	dsnm = grpnm + "/grain_identifier";
	//##MK::write m_container->getRawData() in-place
	u32 = m_container->getCopy();
	cout << "u32.size() " << u32.size() << "\n";
	anno = ioAttributes();
	if ( h5w.nexus_write(
		dsnm,
		io_info({m_container->getSize()}, {m_container->getSize()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		u32, anno ) != MYHDF5_SUCCESS ) { return; }
	u32 = vector<unsigned int>();

	double toc = omp_get_wtime();
	myprofiler.logev("WritingNeXus", (toc - tic));
}


unsigned int microStructureHdl::CountNumberOfSubgrains() {
	cout << __func__ << "\n";
	unsigned int NumberOfSubgrains = 0;
	vector<Grains*>::iterator itG;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {
			vector<SubGrain*>::iterator it;
			for (it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++)
				NumberOfSubgrains++;

		} else {
			NumberOfSubgrains++;
		}
	}
	cout << "Total number of subgrains " << NumberOfSubgrains << endl;
	return NumberOfSubgrains;
}


void microStructureHdl::ReportProfile() {
	cout << endl << "Approximate OpenMP omp_get_wtime report (seconds)" << endl;
	for ( unsigned int e = 0; e < myprofiler.get_entries(); e++ ) {
		cout << myprofiler.titles[e] << ";" << setprecision(6) << myprofiler.times[e] << endl; //separator was "\t\t\t\t"
	}
}


struct NUMANode {
	int num_cpus;
	int numa_cpus[MYMAX_NUMA_NODES];
};


unsigned int my_numa_bitmask_weight(const struct bitmask *mask) {
	unsigned int weight = 0;
	for (unsigned int j = 0; j < mask->size; j++) {
		if (numa_bitmask_isbitset(mask, j)) {
			weight++;
		}
	}
	return weight;
}


void microStructureHdl::InitEnvironment()
{
	cout << __func__ << "\n";
	if (Settings::ExecuteInParallel) {
		Settings::MaximumNumberOfThreads = omp_get_max_threads();
	}
	else {
		Settings::MaximumNumberOfThreads = 1;
		omp_set_num_threads(Settings::MaximumNumberOfThreads);
	}
	m_grainScheduler = new IterativeGrainScheduler(
			Settings::MaximumNumberOfThreads, Settings::NumberOfGrains);

	InitNumaBinding();
}


void microStructureHdl::InitNumaBinding()
{
	cout << __func__ << "\n";
	vector<NUMANode> nodes;
	nodes.reserve(16);
	numa_available();
	// returns a mask of CPUs on which the current task is allowed to run.
	bitmask* mask = numa_get_run_node_mask();
	bitmask* cpus = numa_allocate_cpumask();
	for (unsigned int j = 0; j < mask->size; j++) {
		if (numa_bitmask_isbitset(mask, j)) {
			printf("We are allowed to used node %d\n", j);
			NUMANode node;
			memset(&node, 0xFF, sizeof(node));
			//converts a node number to a bitmask of CPUs.
			//The user must pass a bitmask structure with a mask buffer long enough to represent all possible cpu's
			numa_node_to_cpus(j, cpus);
			node.num_cpus = my_numa_bitmask_weight(cpus);
			int cpuCounter = 0;
			for (unsigned int i = 0; i < cpus->size; i++) {

				if (numa_bitmask_isbitset(cpus, i)
						&& numa_bitmask_isbitset(numa_all_cpus_ptr, i)) {
					node.numa_cpus[cpuCounter] = i;
					cpuCounter++;
				}
			}
			nodes.push_back(node);
		}
	}
	numa_free_cpumask(cpus);

	#pragma omp parallel
	{
		int threadID = omp_get_thread_num();
		for (unsigned int i = 0; i < nodes.size(); i++) {
			if (threadID < nodes.at(i).num_cpus) {
				#pragma omp critical
				{
					cout << "Will bind thread " << omp_get_thread_num() << " to node " << nodes.at(i).numa_cpus[threadID];
					cpu_set_t set;
					CPU_ZERO(&set);
					CPU_SET(nodes.at(i).numa_cpus[threadID], &set);
					int res = sched_setaffinity(0, sizeof(set), &set);
					if ( res == 0 ) {
						cout << ", success" << "\n";
					}
					else {
						cout << "\n";
						cerr << ", failed !" << "\n";
					}
				}
				break;
			}
			threadID -= nodes.at(i).num_cpus;
		}
	}
}


void microStructureHdl::InitializeGrains(
	vector<vector<Eigen::Vector3d>> hulls, vector<double> grainVolume)
{
	cout << __func__ << "\n";
	m_grainScheduler->buildThreadWorkloads(hulls, Settings::NumberOfGridpoints);

	bool exceptionHappened = false;
	string error_message;
	m_grains.resize(Settings::NumberOfGrains + 1);

	#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfGrains) {
				try {
					Grains* newGrain = new Grains(id, hulls[id], m_container,
							grainVolume[id]);
					m_grains[id] = newGrain;

				}
				catch (exception& e) {
					#pragma omp critical
					{
						exceptionHappened = true;
						error_message += string("Grain ")
								+ to_string((unsigned long long) id)
								+ " in its constructor! Reason : " + e.what()
								+ string("\n");
					}
				}
			}
			if (exceptionHappened) {
				throw runtime_error(error_message);
			}
		}
	} //end of parallel region
}
