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

#include "hdf5.h"


using namespace voro;
using namespace std;
using namespace Eigen;

//Definition of the needed functions **********************************************************************************************************************************************************

timeLogger::timeLogger(){
 entries = 0;
}

timeLogger::~timeLogger(){
	titles.clear();
	times.clear();
}

void timeLogger::logev(const string title, double time){
	titles.push_back(title);
	times.push_back(time);
	entries++;
}

unsigned int timeLogger::get_entries( void ){
	return entries;
}




microStructureHdl::microStructureHdl() {

	m_container = NULL;
	m_grainScheduler = NULL;
	m_OrientationSpace = NULL;
	m_ParticlesPositions = NULL;
	m_seqRND = NULL;

	int subgrains = Settings::NumberOfSubgrains;
	if (subgrains <= 0)
		subgrains = 1;

	if (Settings::PlotDimension != E_2D) {
		if (Settings::NumberOfGridpoints == 0) { //perform autodetection of computational grid
			if (Settings::MicroGenMode == E_PREDEFINEDSTRUCTURE) {
				Settings::NumberOfGridpoints = pow(
						(subgrains * Settings::NumberOfGrains * 27), (1 / 3.0))
						* pow(4 / 3 * PI, 1. / 3) / 2
						* Settings::NumberOfPointsPerSubGrain;
			} else {
				Settings::NumberOfGridpoints = pow(
						(subgrains * Settings::NumberOfGrains), (1 / 3.0))
						* pow(4 / 3 * PI, 1. / 3) / 2
						* Settings::NumberOfPointsPerSubGrain;
			}
		} else {
			//no autodetection user-defined hard setting of computational grid
			//already performed during reading XML
		}
		m_container = new DimensionalBuffer<unsigned int>(0, 0, 0,
				Settings::NumberOfGridpoints, Settings::NumberOfGridpoints,
				Settings::NumberOfGridpoints);

	} else {
		if (Settings::NumberOfGridpoints == 0) {
			if (Settings::MicroGenMode == E_PREDEFINEDSTRUCTURE) {
				Settings::NumberOfGridpoints = pow(
						(subgrains * Settings::NumberOfGrains * 9), (1 / 2.0))
						* sqrt(PI) / 2 * Settings::NumberOfPointsPerSubGrain;
			} else {
				Settings::NumberOfGridpoints = pow(
						(subgrains * Settings::NumberOfGrains), (1 / 2.0))
						* sqrt(PI) / 2 * Settings::NumberOfPointsPerSubGrain;

			}
		} else {
			//nothing to do
		}

		m_container = new DimensionalBuffer<unsigned int>(0, 0, 0,
				Settings::NumberOfGridpoints, Settings::NumberOfGridpoints, 1);

	}
	m_h = 1. / Settings::NumberOfGridpoints;
	first_id = 0 + 1;

	m_seqRND = new randomClass( -3000 );
	cout << "Number of Gridpoints: " << Settings::NumberOfGridpoints << " h:" << m_h << endl;
	//vector<randomClass*>* m_threadlocalRND = new vector<randomClass*> [omp_get_max_threads() - 1];
}


microStructureHdl::~microStructureHdl() {
	delete m_container;
	delete m_grainScheduler;
	delete m_OrientationSpace;
	delete m_ParticlesPositions;
	delete m_seqRND;
}


void microStructureHdl::GeneratePolycrystallineStructureOfGrains() {
	double gtime = omp_get_wtime();

	switch (Settings::MicroGenMode) {
		case E_VORONOI: {
			VoroGEN();
			break;
		}
		case E_PREDEFINEDSTRUCTURE: {
			VoroGenPseudoPeriodic();
			break;
		}
		default: {
			break;
		}
	}
	copyContainer();

	myprofiler.logev( "GenerateGrainHull", (omp_get_wtime() - gtime) );

//plotGrains();
}

void microStructureHdl::GenerateSubgrainStructureInEachGrain() {
	double gtime = omp_get_wtime();

	Execute_SubgrainConstruction();

	myprofiler.logev( "GenerateSubgrainHull", (omp_get_wtime() - gtime) );
}


void microStructureHdl::readParticleFile() {
	FILE * OriFromFile = NULL;
	OriFromFile = fopen(Settings::ReadFromFilename.c_str(), "r");
	if ( OriFromFile == NULL ) {
		std::cout << "ParticleFile not existent or unreadable!" << std::endl;
		Settings::StatusHealthy = false;
		return;
	}

	int id, N = 0;
	char c;
// count number of orientations
	do {
		c = fgetc(OriFromFile);
		if (c == '\n')
			N++;
	} while (c != EOF);
	rewind(OriFromFile);

// read over header
	double vol, x, y, z, boxSize;
	double q[4];
	m_OrientationSpace = new vector<myQuaternion>;
	m_ParticlesPositions = new vector<Vector3d>;
	for (int i = 0; i < N-1; i++) {
		if (i == 0) {
			fscanf(OriFromFile, "%lf \n", &boxSize);
			do {
				c = fgetc(OriFromFile);
			} while (c != '\n');
			continue;
		}
		if (i == 1) {
			do {
				c = fgetc(OriFromFile);
			} while (c != '\n');
			continue;
		}
		//else {
			fscanf(OriFromFile,
					"%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n",
					&id, &x, &y, &z, &vol, &q[0], &q[1], &q[2], &q[3]);
			m_OrientationSpace->push_back(myQuaternion(q[0], q[1], q[2], q[3]));
			m_ParticlesPositions->push_back(
					Vector3d(x / boxSize, y / boxSize, z / boxSize));
		//}
	}
	fclose(OriFromFile);

cout << "ParticleFile read successfully with m_OrientationSpace.size = "  << m_OrientationSpace->size() << endl;
//for ( unsigned int j = 0; j < m_OrientationSpace->size(); j++ ) cout << (*m_OrientationSpace)[j].get_q0() << ";" <<  (*m_OrientationSpace)[j].get_q1() << endl;
}

void microStructureHdl::VoroGenPseudoPeriodic() {
	cout << "Started Voro Gen Pseudo Periodic" << endl;
	m_grains.resize((Settings::NumberOfGrains * 27) + 1);
	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
	voronoicell_neighbor c;
	int blocks = (int) (pow((Settings::NumberOfGrains * 27 / 8), (1 / 3.)) + 1);
	if (blocks < 1)
		blocks = 1;
	voro::container con(0, 1, 0, 1, 0, 1, blocks, blocks, blocks, randbedingung,
			randbedingung, randbedingung, 8);

	c_loop_all vl(con);

	/**********************************************************/
	int k_limit = 3;
	if (Settings::PlotDimension == E_2D)
		k_limit = 1;
	int region = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < k_limit; k++) {
				for (int id = 0; id < Settings::NumberOfGrains; id++) {
					double x = (*m_ParticlesPositions)[id][0] / 3.
							+ (i * 1. / 3);
					double y = (*m_ParticlesPositions)[id][1] / 3.
							+ (j * 1. / 3);
					double z = (*m_ParticlesPositions)[id][2] / 3.
							+ (k * 1. / 3);
					if (Settings::PlotDimension == E_2D)
						z = 0.;
					int mid = id + (region * Settings::NumberOfGrains);
					con.put(mid, x, y, z);
				}
				region++;
			}
		}
	}

	/**********************************************************/
	Settings::NumberOfGrains *= 27;
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

		con.draw_particles("VoronoyP.gnu");
		con.draw_cells_gnuplot("VoronoyC.gnu");
		double x, y, z, rx, ry, rz;
		int cell_id;
		int z_limit = Settings::NumberOfGridpoints;
		if (Settings::PlotDimension == E_2D)
			z_limit = 1;
		for (int k = 0; k < z_limit; k++) {
			for (int j = 0; j < Settings::NumberOfGridpoints; j++) {
				for (int i = 0; i < Settings::NumberOfGridpoints; i++) {
					y = j * m_h;
					x = i * m_h;
					z = k * m_h;
					if (con.find_voronoi_cell(x, y, z, rx, ry, rz, cell_id)) {
						unsigned int box_id = cell_id + 1;
						m_container->setValueAt(j, i, k, box_id);
					} else {
						m_container->setValueAt(j, i, k, 0);
					}
				}
			}
		}
	} else {
		throw runtime_error("Voronoi container error at start() method!");
	}

//generate grainboxes by evaluating the initial hulls objects
	initializeGrains(initialHulls, grainVolume);
}


void microStructureHdl::VoroGEN() {
	cout << "Started Voro Gen" << endl;


	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet

	if ( Settings::NumberOfGrains == 1 && Settings::VoronoiPeriodic == true )
		randbedingung = true;

		//	if (randbedingung == false)
//		Settings::NumberOfGridpoints = 1;

	voronoicell_neighbor c; //Voronoi Zelle die ihre Nachbaren kennt
//Erstellung eines Containers mit gewissen Dimensionen und Randbedingungen
	int blocks = (int) (pow((Settings::NumberOfGrains / 8), (1 / 3.)) + 1);
	if (blocks < 1)
		blocks = 1;
	voro::container con(0, 1, 0, 1, 0, 1, blocks, blocks, blocks, randbedingung,
			randbedingung, randbedingung, 8);

	c_loop_all vl(con); //Schleife (Iteration)

	/**********************************************************/
	if ( Settings::GrainAggregation == E_POLYCRYSTAL ) {
	//JitterGrid process
		if ( Settings::GrainShape == E_FLAT ) {
			int Nx, Ny, Nz;
			double scaler = pow( (((double) Settings::NumberOfGrains) * Settings::DefGrainRelDimensionY ) , (1.0/2.0) );
			if ( Settings::PlotDimension == E_3D )
				scaler = pow( (((double) Settings::NumberOfGrains) * Settings::DefGrainRelDimensionY * Settings::DefGrainRelDimensionZ ) , (1.0/3.0) );

			Nx = scaler;
			Ny = scaler / Settings::DefGrainRelDimensionY;
			Nz = 1;
			if ( Settings::PlotDimension == E_3D )
				Nz = scaler / Settings::DefGrainRelDimensionZ;

			double dx = 1.0 / scaler;
			double dy = 1.0 / ( scaler * (Settings::DefGrainRelDimensionX / Settings::DefGrainRelDimensionY) );
			double dz = 1.0 / ( scaler * (Settings::DefGrainRelDimensionX / Settings::DefGrainRelDimensionZ) );
	cout << "JitterGrid FlatGrainMorphology Nx/Ny/Nz = " << Nx << ";" << Ny << ";" << Nz << "---" << dx << ";" << dy << ";" << dz << endl;

			int N = 0; //randomly distored point grid on unit cube
			if ( Settings::PlotDimension == E_2D ) {
				double z = 0.0;
				for( double j = (dy / 2); j < 1.0; j += dy ) {
					for( double i = (dx / 2); i < 1.0; i += dx ) { // 0.5 * dx * m_seqRND->MersenneTwister(), shear band structures
						double x = dx * Settings::RandomnessX * (2*m_seqRND->MersenneTwister() - 1.0) + i;
						double y = dy * Settings::RandomnessY * (2*m_seqRND->MersenneTwister() - 1.0) + j;
						con.put( N, x, y, z);
						N++;
	cout << N << ";" << x << ";" << y << ";" << z << endl;
					}
				}
			}
			else if ( Settings::PlotDimension == E_3D ) {
				for ( double k = (dz / 2); k < 1.0; k += dz ) {
					for( double j = (dy / 2); j < 1.0; j += dy ) {
						for( double i = (dx / 2); i < 1.0; i += dx ) { // 0.5 * dx * m_seqRND->MersenneTwister(), shear band structures
							double x = dx * Settings::RandomnessX * (2*m_seqRND->MersenneTwister() - 1.0) + i;
							double y = dy * Settings::RandomnessY * (2*m_seqRND->MersenneTwister() - 1.0) + j;
							double z = dz * Settings::RandomnessZ * (2*m_seqRND->MersenneTwister() - 1.0) + k;
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

	/*
	//Bands, only implemented for 2D
		if ( Settings::GrainShape == E_FLAT ) {
			int Nx, Ny;
			double scaler = pow( (((double) Settings::NumberOfGrains) * Settings::DefGrainRelDimensionY ) , (1.0/2.0) );

			Nx = scaler;
			Ny = scaler / Settings::DefGrainRelDimensionY;
			double dx = 1.0 / scaler;
			double dy = 1.0 / ( scaler * (Settings::DefGrainRelDimensionX / Settings::DefGrainRelDimensionY) );
	cout << "JitterGrid FlatGrainMorphology Nx/Ny = " << Nx << ";" << Ny << ";" << "---" << dx << ";" << dy << endl;

			int N = 0; //randomly distored point grid on unit cube
			double z = 0.0;

			for( double j = (dy / 2); j < 1.0; j += dy ) {
				for( double i = (dx / 2); i < 1.0; i += dx ) { // 0.5 * dx * m_seqRND->MersenneTwister(), shear band structures
					double x = dx * Settings::RandomnessX * (2*m_seqRND->MersenneTwister() - 1.0) + i;
					double y = dy * Settings::RandomnessY * (2*m_seqRND->MersenneTwister() - 1.0) + j;
					//double z = dz * Settings::RandomnessZ * (2*m_seqRND->MersenneTwister() - 1.0) + k;
					con.put( N, x, y, z);
					N++;

	cout << N << ";" << x << ";" << y << ";" << z << endl;
				}
				//dy =  1.0 / ( scaler * (Settings::DefGrainRelDimensionX / Settings::DefGrainRelDimensionY) ) * (1.0 + (2*m_seqRND->MersenneTwister()-1.0)/4.0);
			}

			Settings::NumberOfGrains = N;
	cout << "JitterGrid-FlatGrainMorphology - I am resetting the number of grains to N/NumberOfGrains = " << N << ";" << Settings::NumberOfGrains << endl;

		}
	*/

		else { //E_GLOBULITIC
			// Randomly add particles into the container
			for (int i = 0; i < Settings::NumberOfGrains; i++) {
				double x = m_seqRND->MersenneTwister();
				double y = m_seqRND->MersenneTwister();
				double z = m_seqRND->MersenneTwister();
				if (Settings::PlotDimension == E_2D)
					z = 0.0;
				con.put(i, x, y, z); //assign each Voronoi cell i its coordinates x,y,z
			}
		}
	}
	else if ( Settings::GrainAggregation == E_BICRYSTAL ) {
		//MK::model utilized for M. K\"uhbach, F. Roters, et. al. nonlocal solving of W.B. Hutchinson's old GB nucleation bicrystal experiments on iron
		//the spectral solver applies always periodic boundary conditions so strictly speaking the lower parts see
		//x || RD, y || TD, z || ND
		double x, y, z;
		int i = 0;
		//place upper container for the bicrystal
		x = 0.5;
		y = 0.5;
		z = 0.75;
		con.put( i, x, y, z);

		//place lower container for the bicrystal
		x = 0.5;
		y = 0.5;
		z = 0.25;
		i++;
		con.put( i, x, y, z);
	}
	else {
		//##MK::E_SINGLECRYSTAL currently not implemented
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
		if (Settings::PlotDimension == E_2D)
			z_limit = 1;
		for (int k = 0; k < z_limit; k++) {
			for (int j = 0; j < Settings::NumberOfGridpoints; j++) {
				for (int i = 0; i < Settings::NumberOfGridpoints; i++) {
					y = j * m_h;
					x = i * m_h;
					z = k * m_h;
					if (con.find_voronoi_cell(x, y, z, rx, ry, rz, cell_id)) {
						unsigned int box_id = cell_id + 1;
						m_container->setValueAt(j, i, k, box_id);
					} else {
						m_container->setValueAt(j, i, k, 0);
					}
				}
			}
		}
	} else {
		throw runtime_error("Voronoy container error at start() method!");
	}

//generate grainboxes by evaluating the initial hulls objects
	initializeGrains(initialHulls, grainVolume);

	cout << "Parent grain geometry constructed " << endl;
}

void microStructureHdl::Execute_SubgrainConstruction() {
	cout << "Start subgrain construction " << endl;
#pragma omp parallel
	{
		randomClass* threadlocalRNG = new randomClass();
		//threadlocalRNG->initPM(-1 * (omp_get_thread_num()) - 1);
		uint32_t overflowint = (uint32_t) pow(2.0, 31);
		threadlocalRNG->initMT( overflowint - omp_get_thread_num() - 1 );
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
//		cout << "Thread: " << omp_get_thread_num() << " has received workload" << endl;
		for (auto id : workload) {
//			cout << id << endl;
			if (id <= Settings::NumberOfGrains)
				m_grains[id]->SubGrainConstructor(*threadlocalRNG);
		}
		delete threadlocalRNG;
	}
}

void microStructureHdl::copyContainer() {
	for (auto it : m_grains)
		if (it != NULL)
			it->copyContainerToGrain(m_container);
}

void microStructureHdl::find_neighbors() {
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

void microStructureHdl::DistributeGrainOriAndSEE() {
	double gtime = omp_get_wtime();

	//###MK::readPreferenceOrientationFromFile();
	cout << "Sample grain orientations" << endl;
	vector<Grains*>::iterator it;
	for (it = ++m_grains.begin(); it != m_grains.end(); it++) {
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
		else if (Settings::MicroGenMode == E_PREDEFINEDSTRUCTURE) {					//from predefined orientations
			ori = (*m_OrientationSpace)[((*it)->get_ID() - 1) % 27];
		}
		else {
			ori.randomOriShoemakeQuat(*m_seqRND);									//random orientations
		}

		//assign grain the orientation and scattering properties
		(*it)->set_Orientation(ori);
		(*it)->set_PreforiProperties( findNextPreferenceOrientation(ori) );
		double gr_see_mu = (*it)->get_SEEFromPrefori();
		double gr_see_sig = (*it)->get_SEEGrainScatterFromPrefOri();
//cout << "Grain = " << (*it)->get_ID() << " prefori mu/sig/orisig = " << gr_see_mu << ";" << gr_see_sig << ";" << (*it)->get_OriScatterFromPrefOri() << "----prefori q0 " <<  (*it)->get_PrefOriQuatQ0() << ";" << (*it)->get_PrefOriQuatQ1() << ";" << (*it)->get_PrefOriQuatQ2() << ";" << (*it)->get_PrefOriQuatQ3() << endl;
		(*it)->set_SEE( m_seqRND->r4_nor(gr_see_mu, gr_see_sig) ); //set stored elastic energy of grain to a specific value from a normal distribution
	}
	cout << "Sample sub-grain orientation from reference orientation of grain" << endl;

	myprofiler.logev("DistriGrainPrefOris", (omp_get_wtime() - gtime));
}

void microStructureHdl::DistributeSubgrainOrientations() {
	double gtime = omp_get_wtime();

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

	myprofiler.logev("DistrSubgrainOris", (omp_get_wtime() - gtime) );
}


void microStructureHdl::ReadAdditionalInputFiles() {
	double gtime = omp_get_wtime();

	readParticleFile();
	readPreferenceOrientationFromFile();

	myprofiler.logev("ReadAdditionalInput", (omp_get_wtime() - gtime));
}


void microStructureHdl::readPreferenceOrientationFromFile() {
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
				myPreferenceOri(phi1 * PI / 180.0 , PHI * PI / 180.0 ,
						phi2 * PI / 180.0 , oriScatter * PI / 180.0 , SEE, SEEGrainScatter, SEESubgrainScatter, RelSubgrainSizeScaler)); //define a preference orientation with specific instead of default properties in terms of mu, sigma and oriscatter
cout << "PrefRefs = " << phi1 << ";" << PHI << ";" << phi2 << ";" << oriScatter << ";" << SEE << ";" << SEEGrainScatter << ";" << SEESubgrainScatter << ";" << RelSubgrainSizeScaler << endl;
	}
	fclose(file);
cout << "PrefRefDefault = " << Settings::SubgrainOriScatter << ";" << (Settings::StoredElasticEnergyMax+Settings::StoredElasticEnergyMin)*0.5 << ";" << Settings::StoredElasticScatterGrain << ";" << 1.0 << endl;
}


myPreferenceOri microStructureHdl::findNextPreferenceOrientation(
		myQuaternion ori) {
	if (Settings::TextureGEN == E_USE_PREFERENCEORI) {
		double min = 15.0 * PI / 180.0;
		myPreferenceOri* i;
		unsigned int idx = 0;
		for ( unsigned int p = 0; p < m_PreferenceOrientations->size(); p++ ) {
			double misori = ori.misorientationCubicQxQ( &((*m_PreferenceOrientations)[p].ori) );
			if ( misori <= min ) {
				min = misori;
				idx = p;
			}
		}
		if (min <= 10.0 * PI / 180.0) {
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

void microStructureHdl::DistributeSubgrainSEE() {
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


void microStructureHdl::RehashGrainIDs() {
	double timer = omp_get_wtime();
	cout << ">Rehashing grain IDs" << endl;

	//classical (Mie\ss{}en and K\"uhbach 2015,2016,2017 grain IDs start with 1
	int offset = 0;

	//MK::consider a potential breaking of the periodicity by adding pair(s) of "air" grains 1,2,3,4,5,6
	if ( Settings::BreakPerX == true )	offset += 2;
	if ( Settings::BreakPerY == true )	offset += 2;
	if ( Settings::BreakPerZ == true )	offset += 2;

	int newoffset = 0;
	vector<Grains*>::iterator itG;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {
			newoffset = (*itG)->copySubgrainsToGlobalContainer(m_container, offset); // implicit rehashing of all ID's
			offset = --newoffset;
			vector<SubGrain*>::iterator it;
		} else {
			newoffset = (*itG)->copySubgrainsToGlobalContainer(m_container, offset);
			offset = newoffset;
		}
	}

	myprofiler.logev("RehashingIDs", (omp_get_wtime() - timer) );
}


void microStructureHdl::BreakPeriodicity() {

	double timer = omp_get_wtime();
	//add "air" grain pair at the boundary to break the symmetry of the domain
	unsigned int offset = 1;

	if ( Settings::BreakPerX == true ) {
		set_first_id( true, 2 );
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
		set_first_id( true, 2 );
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
		set_first_id( true, 2 );
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
					double newVolume = 0.0;
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
				double newVolume = 0.0;
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


void microStructureHdl::SaveDataGraGeLeS() {
	double asciiftime = omp_get_wtime();
	cout << "Pipe data into file" << endl;
	stringstream filename;
	filename << "Microstructure.SimID." << Settings::SimID << ".uds";
	ofstream file;
	file.open(filename.str().c_str());

	//file header
	file << "|| Settings || *(sf) || Parameter, Value" << endl;
	file << "FirstID" << "\t" << this->get_first_id() << endl;
	file << "DX" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	file << "DY" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	if (Settings::PlotDimension != E_2D)
		file << "DZ" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	file << "NX" << "\t" << Settings::NumberOfGridpoints << endl;
	file << "NY" << "\t" << Settings::NumberOfGridpoints << endl;
	if (Settings::PlotDimension != E_2D)
		file << "NZ" << "\t" << Settings::NumberOfGridpoints << endl;
	file << "NGrains" << "\t" << this->CountNumberOfSubgrains() << endl;

	if (Settings::PlotDimension != E_2D) {
		file
				<< "|| Points || *(iffffffiiiiiiff) || ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, vol, stored"
				<< endl;
	} else {
		file
				<< "|| Points || *(ifffffiiiiff) || ID, x, y, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, vol, stored"
				<< endl;
	}

	//##int offset = 0;
	//##int newoffset = 0;
	vector<Grains*>::iterator itG;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0)
				&& (*itG)->m_SubGrains.size() > 1) {

			//##newoffset = (*itG)->copySubgrainsToGlobalContainer(m_container, offset); // implicit rehashing of all ID's
			//##offset = --newoffset;
			vector<SubGrain*>::iterator it;
			for (it = ++((*itG)->m_SubGrains.begin());
					it != (*itG)->m_SubGrains.end(); it++) {
				// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
				myQuaternion* ori = (*it)->get_Orientation();
				double* bunge = ori->Quaternion2Euler();
				file << (*it)->get_ID() << "\t"
						<< (((*it)->getMaxX() - (*it)->getMinX()) / 2
								+ (*it)->getMinX()) << "\t"
						<< (((*it)->getMaxY() - (*it)->getMinY()) / 2
								+ (*it)->getMinY()) << "\t";
				if (Settings::PlotDimension != E_2D)
					file
							<< (((*it)->getMaxZ() - (*it)->getMinZ()) / 2
									+ (*it)->getMinZ()) << "\t";
				file << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2] << "\t"
						<< (*it)->getMinX() << "\t" << (*it)->getMaxX() << "\t"
						<< (*it)->getMinY() << "\t" << (*it)->getMaxY() << "\t";
				if (Settings::PlotDimension != E_2D)
					file << (*it)->getMinZ() << "\t" << (*it)->getMaxZ()
							<< "\t";
				file << (*it)->get_Volume() << "\t" << (*it)->get_SEE() << endl;
				delete[] bunge;
			}
		} else {
			// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
			myQuaternion ori = (*itG)->getOri();
			double* bunge = ori.Quaternion2Euler();
			//##newoffset = (*itG)->copySubgrainsToGlobalContainer(m_container, offset);
			//##offset = newoffset;
			file << (*itG)->get_ID() << "\t"
					<< (((*itG)->getMaxX() - (*itG)->getMinX()) / 2
							+ (*itG)->getMinX()) << "\t"
					<< (((*itG)->getMaxY() - (*itG)->getMinY()) / 2
							+ (*itG)->getMinY()) << "\t";
			if (Settings::PlotDimension != E_2D)
				file
						<< (((*itG)->getMaxZ() - (*itG)->getMinZ()) / 2
								+ (*itG)->getMinZ()) << "\t";
			file << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2] << "\t"
					<< (*itG)->getMinX() << "\t" << (*itG)->getMaxX() << "\t"
					<< (*itG)->getMinY() << "\t" << (*itG)->getMaxY() << "\t";
			if (Settings::PlotDimension != E_2D)
				file << (*itG)->getMinZ() << "\t" << (*itG)->getMaxZ() << "\t"
						<< int((*itG)->get_Volume() / m_h / m_h / m_h + 0.5)
						<< "\t";
			file << int((*itG)->get_Volume() / m_h / m_h + 0.5) << "\t"
					<< (*itG)->getSEE() << endl;
			delete[] bunge;
		}
	}
//write grain infos

	file.close();

	myprofiler.logev("WritingASCIIUDS", (omp_get_wtime() - asciiftime) );

	//##cout << "created and piped: " << offset << " subgrains to binary file. Move on ... " << endl;
	cout << "created and piped subgrains to binary file. Move on ... " << endl;

	double binftime = omp_get_wtime();

	int size;
	if (Settings::PlotDimension == E_3D) {
		size = pow(Settings::NumberOfGridpoints, 3.0);
		if ( size != CUBE(Settings::NumberOfGridpoints) ) { cout << "Size inconsistency in 3D!" << endl; return; }
	}
	else if (Settings::PlotDimension == E_2D) {
		size = pow(Settings::NumberOfGridpoints, 2.0);
		if ( size != SQR(Settings::NumberOfGridpoints) )  { cout << "Size inconsistency in 2D!" << endl; return; }
	}
	else { cout << "Nothing implemented for other than 2D or 3D microstructures!" << endl; return; }

	stringstream filename2;
	filename2.str("");
	filename2 << "Microstructure.SimID." << Settings::SimID << ".GrainIDs.";
	if ( Settings::PlotDimension == E_2D )
		filename2 << "2D." << Settings::NumberOfGridpoints << ".raw";
	else
		filename2 << "3D." << Settings::NumberOfGridpoints << ".raw";
	FILE* binaryFile = NULL;

	binaryFile = fopen(filename2.str().c_str(), "w");
	if ( binaryFile == NULL ) { cout << "Failure to open binary file!" << endl; return; }
	fwrite(m_container->getRawData(), sizeof(unsigned int), size, binaryFile);
	fclose(binaryFile);

	myprofiler.logev("WritingBINARYContainer", (omp_get_wtime() - binftime));
}


void microStructureHdl::SaveDataDAMASKMatConfig() {
	cout << "Writing material config file for data import into DAMASK NONLOCAL" << endl;

	stringstream damask_matconfig_fname;
	damask_matconfig_fname << "Microstructure.SimID." << Settings::SimID << ".DAMASKMaterial.config";
	ofstream damask_matconfig;
	damask_matconfig.open ( damask_matconfig_fname.str().c_str() );

	//header
	damask_matconfig << "##IMMMicrostructureGenerator auto-generated DAMASK material.config file!\n";
	damask_matconfig << "##Mind that not all parts of the material.config file required to run a successful DAMASK simulation are generated by this tool!!!\n";
	damask_matconfig << "##We assume at the moment only one crystallite section and only one phase section and volume fraction 1.0\n";
	damask_matconfig << "##End of comments for IMMMicrostructureGenerator autogenerated DAMASK material.config file\n";
	damask_matconfig << "\n\n\n";

	//write microstructure section
	damask_matconfig << "#-------------------#\n";
	damask_matconfig << "<microstructure>\n";
	damask_matconfig << "#-------------------#\n";

	//are there any grains of "air"? //MK::1. crystalline is one, phase is one, texture is one
	unsigned int offset = 0;
	unsigned int texoffset = 0;
	if ( Settings::BreakPerX == true ) {
		offset++;
		damask_matconfig << "[Grain" << offset << "]\n";
		damask_matconfig << "crystallite 1\n"; damask_matconfig << "(constituent)  phase 1   texture 1   fraction 1.0\n";
		offset++;
		damask_matconfig << "[Grain" << offset << "]\n";
		damask_matconfig << "crystallite 1\n"; damask_matconfig << "(constituent)  phase 1   texture 1   fraction 1.0\n";
		texoffset += 1; //to assure that in case of not broken periodicity texture becomes 2
	}
	if ( Settings::BreakPerY == true ) {
		offset++;
		damask_matconfig << "[Grain" << offset << "]\n";
		damask_matconfig << "crystallite 1\n"; damask_matconfig << "(constituent)  phase 1   texture 1   fraction 1.0\n";
		offset++;
		damask_matconfig << "[Grain" << offset << "]\n";
		damask_matconfig << "crystallite 1\n"; damask_matconfig << "(constituent)  phase 1   texture 1   fraction 1.0\n";
		texoffset += 2;
	}
	if ( Settings::BreakPerZ == true ) {
		offset++;
		damask_matconfig << "[Grain" << offset << "]\n";
		damask_matconfig << "crystallite 1\n"; damask_matconfig << "(constituent)  phase 1   texture 1   fraction 1.0\n";
		offset++;
		damask_matconfig << "[Grain" << offset << "]\n";
		damask_matconfig << "crystallite 1\n"; damask_matconfig << "(constituent)  phase 1   texture 1   fraction 1.0\n";
		texoffset += 2;
	}

	//write all subgrains as DAMASK grains
	unsigned int grainid, constitutive;
	vector<Grains*>::iterator itG;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
			vector<SubGrain*>::iterator it;
			for (it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++) {
				grainid = (*it)->get_ID();			//IDs take care of potential periodicity
				constitutive = grainid - texoffset;	//##MK::make a new texture component for each grain

				damask_matconfig << "[Grain" << grainid << "]\n";
				damask_matconfig << "crystallite " << constitutive << "\n";
				damask_matconfig << "(constituent)  phase 1   texture " << constitutive << "   fraction 1.0\n";
			}
		} else {
			grainid = (*itG)->get_ID();				//##MK::GraGeLeS ID start at 1...
			constitutive = grainid - texoffset;		//##MK::make a new texture component for each grain

				damask_matconfig << "[Grain" << grainid << "]\n";
				damask_matconfig << "crystallite " << constitutive << "\n";
				damask_matconfig << "(constituent)  phase 1   texture " << constitutive << "   fraction 1.0\n";
		}
	}
	damask_matconfig << endl << endl;


	damask_matconfig << "#-------------------#\n";
	damask_matconfig << "<texture>\n";
	damask_matconfig << "#-------------------#\n";

	if ( texoffset > 0 ) {
		damask_matconfig << "[Texture1]\n";
		damask_matconfig << "(gauss)  phi1 0.0    Phi 0.0    phi2 0.0   scatter 0.0   fraction 1.0\n";
	}
	double rad2deg = 180.0 / _PI_;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if (Settings::NumberOfSubgrains != 0 && (*itG)->m_SubGrains.size() > 1) {
			vector<SubGrain*>::iterator it;
			for (it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++) {
				grainid = (*it)->get_ID();			//##MK::GraGeLeS ID start at 1...
				constitutive = grainid - texoffset;				//##MK::make a new texture component for each grain
				myQuaternion* ori = (*it)->get_Orientation();
				double* bunge = ori->Quaternion2Euler();

				damask_matconfig << "[Grain" << constitutive << "]\n";
				damask_matconfig << "(gauss)  phi1 " << setprecision(8) << rad2deg * bunge[0] << "    Phi " << setprecision(8) << rad2deg * bunge[1] << "    phi2 " << setprecision(8) << rad2deg * bunge[2] << "   scatter 0.0   fraction 1.0\n";
				delete [] bunge;
			}
		} else {
			grainid = (*itG)->get_ID();			//##MK::GraGeLeS ID start at 1...
			constitutive = grainid - texoffset;				//##MK::make a new texture component for each grain
			myQuaternion ori = (*itG)->getOri();
			double* bunge = ori.Quaternion2Euler();

			damask_matconfig << "[Grain" << constitutive << "]\n";
			damask_matconfig << "(gauss)  phi1 " << setprecision(8) << rad2deg * bunge[0] << "    Phi " << setprecision(8) << rad2deg * bunge[1] << "    phi2 " << setprecision(8) << rad2deg * bunge[2] << "   scatter 0.0   fraction 1.0\n";
			delete [] bunge;
		}
	}
	damask_matconfig << endl << endl;


	damask_matconfig.flush();
	damask_matconfig.close();
}


void microStructureHdl::SaveDataDAMASKGeometry() {
	double geomftime = omp_get_wtime();

	int size;
	if (Settings::PlotDimension == E_3D) {
		size = pow(Settings::NumberOfGridpoints, 3.0);
		if ( size != CUBE(Settings::NumberOfGridpoints) ) { cout << "Size inconsistency in 3D!" << endl; return; }
	}
	//else if (Settings::PlotDimension == E_2D) {
	//	cout << "##MK::DAMASK support currently only in 3D!" << endl; return;
	//}
	//	size = pow(Settings::NumberOfGridpoints, 2.0);
	//	if ( size != SQR(Settings::NumberOfGridpoints) )  { cout << "Size inconsistency in 2D!" << endl; return; }
	//}
	else { cout << "Nothing implemented for other than 2D or 3D microstructures!" << endl; return; }

	//write geomFile
	stringstream damask_geom_fname;
	damask_geom_fname << "Microstructure.SimID." << Settings::SimID << ".DAMASKGeometry.geom";
	ofstream damask_geom;
	damask_geom.open ( damask_geom_fname.str().c_str() );

	//header
	damask_geom << "//IMMMicrostructureGenerator auto-generated DAMASK geometry file!\n";
	damask_geom << "//Mind that coordinate systems in IMMMicrostructureGenerator and DAMASK .......!\n";
	damask_geom << "//IMMMicrostructureGenerator is ..........\n";
	damask_geom << "//DAMASK's is right-handed with +x pointing right, +y pointing inwards, and +z pointing upwards\n";
	damask_geom << "//-->This in DAMASK ..... into +x || RD, but ND and TD change to +z || ND and while +y || TD!\n";
	damask_geom << "//End of comments for IMMMicrostructureGenerator autogenerated DAMASK geometry file\n";
	damask_geom << "\n\n\n";

	//write DAMASK header
	damask_geom << "5\theader\n";
	damask_geom << "grid\ta " << Settings::NumberOfGridpoints << "\tb " << Settings::NumberOfGridpoints << "\tc " << Settings::NumberOfGridpoints << "\n";
	damask_geom << "size\tx 1.000000\ty 1.000000\t z 1.000000\n";
	damask_geom << "origin\tx 0.000000\ty 0.000000\t z 0.000000\n";

	//count number of
	unsigned int NSubgrainsInTotal = 0;
	vector<Grains*>::iterator itG;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {
			vector<SubGrain*>::iterator it;
			for (it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++) {
				NSubgrainsInTotal++;
			}
			continue;
		}
		NSubgrainsInTotal++;
	}
	damask_geom << "microstructures " << NSubgrainsInTotal << "\n";
	damask_geom << "homogenization\t1" << endl;
	//MK::these make the output not portable but one can change here the structure accordingly...

	unsigned int* rawdata = NULL;
	rawdata = m_container->getRawData();
	if ( rawdata == NULL ) { cout << "Unable to draw rawdata!" << endl; damask_geom.flush(); damask_geom.close(); return; }

	unsigned int nlines = Settings::NumberOfGridpoints * Settings::NumberOfGridpoints;
	unsigned int roffset = 0; //##MK::all currently possible computational domains are well below 2^32-1 in size...
	for ( unsigned int line = 0; line < nlines; line++ ) {
		roffset = line * Settings::NumberOfGridpoints;
		for ( unsigned int vx = 0; vx < (Settings::NumberOfGridpoints-1); vx++ ) {
			damask_geom << rawdata[roffset+vx] << " ";
		}
		damask_geom << rawdata[roffset+Settings::NumberOfGridpoints-1] << endl; //MK::to avoid additional blank space before LF
	}

	damask_geom.flush();
	damask_geom.close();
}


void microStructureHdl::SaveDataDAMASK(){

	if ( Settings::IO_DAMASK == false )
		return;

	double matconfigftime = omp_get_wtime();

	SaveDataDAMASKMatConfig();

	myprofiler.logev("WritingDAMASKMatConfig", (omp_get_wtime() - matconfigftime) );

	double geomftime = omp_get_wtime();

	SaveDataDAMASKGeometry();

	myprofiler.logev("WritingDAMASKGeom", (omp_get_wtime() - geomftime) );
}


void microStructureHdl::SaveDetailedDiagnosticsASCII() {
	double asciidiagn = omp_get_wtime();
	cout << "Pipe diagnostics into ASCII file" << endl;
	stringstream filename;
	filename << "MicrostructureDiagnostics.SimID." << Settings::SimID << ".uds";
	ofstream file;
	file.open(filename.str().c_str());

	//write header
	file << "|| Settings || *(sf) || Parameter, Value" << endl;
	file << "FirstID" << "\t" << this->get_first_id() << endl;
	file << "DX" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	file << "DY" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	if (Settings::PlotDimension != E_2D)
		file << "DZ" << "\t" << Settings::NumberOfPointsPerSubGrain << endl;
	file << "NX" << "\t" << Settings::NumberOfGridpoints << endl;
	file << "NY" << "\t" << Settings::NumberOfGridpoints << endl;
	if (Settings::PlotDimension != E_2D)
		file << "NZ" << "\t" << Settings::NumberOfGridpoints << endl;
	file << "NGrains" << "\t" <<this->CountNumberOfSubgrains() << endl;

	if (Settings::PlotDimension != E_2D) {
		file
				<< "|| Subgrains || *(iffffffiiiiiiffifffffff) || ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, vol, stored, parent ID, disori to parent meas, disori to parent target, parent bunge1, parent bunge2, parent bunge3, parent stored, parent stored grain scatter, parent stored subgrain scatter, parent size scaler"
				<< endl;
	} else {
		file
				<< "|| Subgrains || *(ifffffiiiiffifffffff) || ID, x, y, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, vol, stored, parent ID, disori to parent meas, disori to parent target, parent bunge1, parent bunge2, parent bunge3, parent stored, parent stored grain scatter, parent subgrain scatter, parent size scaler"
				<< endl;
	}

	int offset = 0;
	vector<Grains*>::iterator itG;

	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {

			vector<SubGrain*>::iterator it;

			myQuaternion parentori = (*itG)->getOri();
			double* parentbunge = parentori.Quaternion2Euler();
			double oriscatter = (*itG)->get_OriScatterFromPrefOri();
			double see = (*itG)->getSEE();
			double seegrscatter = (*itG)->get_SEEGrainScatterFromPrefOri();
			double seesgrscatter = (*itG)->get_SEESubgrainScatterFromPrefOri();
			double sizescaler = (*itG)->get_RelSizeScalingFromPrefori();

			for (it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++) {
				// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
				myQuaternion* ori = (*it)->get_Orientation();
				double* bunge = ori->Quaternion2Euler();
				file << (*it)->get_ID() << "\t" << (((*it)->getMaxX() - (*it)->getMinX()) / 2 + (*it)->getMinX()) << "\t" << (((*it)->getMaxY() - (*it)->getMinY()) / 2 + (*it)->getMinY()) << "\t";
				if (Settings::PlotDimension != E_2D) { file << (((*it)->getMaxZ() - (*it)->getMinZ()) / 2 + (*it)->getMinZ()) << "\t"; }
				file << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2] << "\t" << (*it)->getMinX() << "\t" << (*it)->getMaxX() << "\t" << (*it)->getMinY() << "\t" << (*it)->getMaxY() << "\t";
				if (Settings::PlotDimension != E_2D) { file << (*it)->getMinZ() << "\t" << (*it)->getMaxZ() << "\t"; }
				file << (*it)->get_Volume() << "\t" << (*it)->get_SEE() << "\t" << (*itG)->get_oldID() << "\t" << (double) parentori.misorientationCubicQxQ( ori ) << "\t" << oriscatter << "\t" << parentbunge[0] << "\t" << parentbunge[1] << "\t" << parentbunge[2] << "\t" << see << "\t" << seegrscatter << "\t" << seesgrscatter << "\t" << sizescaler << endl;
				delete[] bunge;
			}

			delete[] parentbunge;

		} else {
			// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
			myQuaternion ori = (*itG)->getOri();
			double* bunge = ori.Quaternion2Euler();

			file << (*itG)->get_ID() << "\t" << (((*itG)->getMaxX() - (*itG)->getMinX()) / 2 + (*itG)->getMinX()) << "\t" << (((*itG)->getMaxY() - (*itG)->getMinY()) / 2 + (*itG)->getMinY()) << "\t";
			if (Settings::PlotDimension != E_2D) { file << (((*itG)->getMaxZ() - (*itG)->getMinZ()) / 2 + (*itG)->getMinZ()) << "\t"; }
			file << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2] << "\t"
					<< (*itG)->getMinX() << "\t" << (*itG)->getMaxX() << "\t"
					<< (*itG)->getMinY() << "\t" << (*itG)->getMaxY() << "\t";
			if (Settings::PlotDimension != E_2D) { file << (*itG)->getMinZ() << "\t" << (*itG)->getMaxZ() << "\t"; }
			file << int((*itG)->get_Volume() / m_h / m_h + 0.5) << "\t" << (*itG)->getSEE() << "\t" << (*itG)->get_oldID() << "\t" << 0.0 << "\t" << (*itG)->get_OriScatterFromPrefOri() << "\t" << bunge[0] << "\t" << bunge[1] << "\t" << bunge[2] << "\t" << (*itG)->getSEE() << "\t" << (*itG)->get_SEEGrainScatterFromPrefOri() << "\t" << (*itG)->get_SEESubgrainScatterFromPrefOri() << "\t" << (*itG)->get_RelSizeScalingFromPrefori() << endl;
			delete[] bunge;
		}
	}

	file.close();

	myprofiler.logev("WritingASCIIDiagnostics", (omp_get_wtime() - asciidiagn) );
}


void microStructureHdl::SaveDetailedDiagnosticsBINARY() {
	//writes one row of doubles per sub-grain to ease the importing of GB-sized diagnostics into MATLAB
	double binarydiagn = omp_get_wtime();
	cout << "Pipe diagnostics into binary file sizeof(BinDiagData) = " << sizeof(BinDiagData) << " bytes." << endl;

	stringstream smartfname;
	smartfname << "Diagnostics.SimID." << Settings::SimID << ".FirstID." << this->get_first_id() << ".";
	if ( Settings::PlotDimension == E_2D ) {
		smartfname << "DIM.2.DXY." << Settings::NumberOfPointsPerSubGrain << ".NXY." << Settings::NumberOfGridpoints;
	}
	else {
		smartfname << "DIM.3.DXYZ." << Settings::NumberOfPointsPerSubGrain << ".NXYZ." << Settings::NumberOfGridpoints;
	}
	smartfname << ".NSGR." << this->CountNumberOfSubgrains() << ".bin";

	FILE* binaryDiagnostics = fopen(smartfname.str().c_str(), "wb");
	//format is as the order in the BinDiagData struct, hence no distinction is made between 2D and 3D, hence readers apply to both files

	vector<Grains*>::iterator itG;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {

			vector<SubGrain*>::iterator it;
			struct BinDiagData asgr;
			asgr.ClearingComplete(); //##MK::safety

			myQuaternion parentori = (*itG)->getOri();
			double* parentbunge = parentori.Quaternion2Euler();

			asgr.dis2parent_targ = (*itG)->get_OriScatterFromPrefOri(); //ascii::oriscatter
			asgr.parent_stored = (*itG)->getSEE();
			asgr.parent_stored_grsc = (*itG)->get_SEEGrainScatterFromPrefOri();
			asgr.parent_stored_sgrsc = (*itG)->get_SEESubgrainScatterFromPrefOri();
			asgr.parent_size_sc = (*itG)->get_RelSizeScalingFromPrefori();
			asgr.parentid = (*itG)->get_oldID();
			asgr.parent_bunge1 = parentbunge[0];
			asgr.parent_bunge2 = parentbunge[1];
			asgr.parent_bunge3 = parentbunge[2];

			for (it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++) {
				asgr.ClearingForSubgrains(); //##MK::safety

				// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
				myQuaternion* ori = (*it)->get_Orientation();
				double* bunge = ori->Quaternion2Euler();

				asgr.id = (*it)->get_ID();
				asgr.x =  (((*it)->getMaxX() - (*it)->getMinX()) / 2 + (*it)->getMinX());
				asgr.y =  (((*it)->getMaxY() - (*it)->getMinY()) / 2 + (*it)->getMinY());
				asgr.z = 0;
				if (Settings::PlotDimension != E_2D) { //global setting
					asgr.z = (((*it)->getMaxZ() - (*it)->getMinZ()) / 2 + (*it)->getMinZ());
				}
				asgr.xmi = (*it)->getMinX();
				asgr.xmx = (*it)->getMaxX();
				asgr.ymi = (*it)->getMinY();
				asgr.ymx = (*it)->getMaxY();
				asgr.zmi = 0;
				asgr.zmx = 0;
				if (Settings::PlotDimension != E_2D) { //global setting
					asgr.zmi = (*it)->getMinZ();
					asgr.zmx = (*it)->getMaxZ();
				}
				asgr.bunge1 = bunge[0];
				asgr.bunge2 = bunge[1];
				asgr.bunge3 = bunge[2];
				asgr.volume = (*it)->get_Volume();
				asgr.stored = (*it)->get_SEE();
				asgr.dis2parent_meas = (double) parentori.misorientationCubicQxQ( ori );

				fwrite(&asgr, sizeof(BinDiagData), 1, binaryDiagnostics);

				delete[] bunge;
			}

			delete[] parentbunge;

		} else {
			struct BinDiagData agr;
			agr.ClearingComplete(); //##MK::safety

			// ID, x, y, z, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, zmin, zmax, volume, stored
			myQuaternion ori = (*itG)->getOri();
			double* bunge = ori.Quaternion2Euler();

			agr.bunge1 = bunge[0];
			agr.bunge2 = bunge[1];
			agr.bunge3 = bunge[2];
			agr.parent_bunge1 = bunge[0]; //because the parent has no sub-grains
			agr.parent_bunge2 = bunge[1];
			agr.parent_bunge3 = bunge[2];

			agr.dis2parent_meas = 0.0;
			agr.dis2parent_targ =  (*itG)->get_OriScatterFromPrefOri();
			agr.parent_stored = (*itG)->getSEE();
			agr.parent_stored_grsc = (*itG)->get_SEEGrainScatterFromPrefOri();
			agr.parent_stored_sgrsc = (*itG)->get_SEESubgrainScatterFromPrefOri();
			agr.parent_size_sc = (*itG)->get_RelSizeScalingFromPrefori();

			agr.id = (*itG)->get_ID();
			agr.x =  (((*itG)->getMaxX() - (*itG)->getMinX()) / 2 + (*itG)->getMinX());
			agr.y =  (((*itG)->getMaxY() - (*itG)->getMinY()) / 2 + (*itG)->getMinY());
			agr.z = 0;
			if (Settings::PlotDimension != E_2D) { //globally set, i.e.  wont change for any grain
				agr.z = (((*itG)->getMaxZ() - (*itG)->getMinZ()) / 2 + (*itG)->getMinZ());
			}
			agr.xmi = (*itG)->getMinX();
			agr.xmx = (*itG)->getMaxX();
			agr.ymi = (*itG)->getMinY();
			agr.ymx = (*itG)->getMaxY();
			agr.zmi = 0;
			agr.zmx = 0;
			if (Settings::PlotDimension != E_2D) { //globally set, i.e. wont change for any grain
				agr.zmi = (*itG)->getMinZ();
				agr.zmx = (*itG)->getMaxZ();
			}
			agr.volume = ((*itG)->get_Volume() / m_h / m_h + 0.5);
			agr.stored = (*itG)->getSEE();
			agr.parentid = (*itG)->get_oldID();
			agr.dis2parent_meas = 0.0;

			fwrite(&agr, sizeof(BinDiagData), 1, binaryDiagnostics);

			delete[] bunge;
		}
	} //next parent grain

	fclose(binaryDiagnostics);
	myprofiler.logev("WritingBINARYDiagnostics", (omp_get_wtime() - binarydiagn) );
}


void microStructureHdl::SaveParenthood() {
	//stores all subgrain IDs and their parent grain IDs to become in a postprocessing of simulation data capable of
	//identifying all sub-grains which are adjacent to a grain-grain boundary even though no neighbors were kept track of
	double ptime = omp_get_wtime();

	//count grains
	vector<Grains*>::iterator itG;
	vector<SubGrain*>::iterator it;
	unsigned int NrGrains = 0;
	if (Settings::NumberOfSubgrains != 0)
		for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
			if ((*itG)->m_SubGrains.size() > 1)
				NrGrains += ((*itG)->m_SubGrains.size() - 1);
			else
				NrGrains++;
		}
	else
		NrGrains = Settings::NumberOfGrains;

	unsigned int* rawdata = NULL;
	rawdata = new unsigned int[2*NrGrains]; //SubgrainID1, ParentID1, SubgrainID2, ParentID2, ..., SubgrainIDN, ParentIDN

	//populate
	unsigned int i = 0;
	for ( itG = ++m_grains.begin(); itG != m_grains.end(); itG++ ) {
		if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {
			for ( it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++ ) {
				rawdata[i] = (*it)->get_ID();
				rawdata[i+1] = (*itG)->get_oldID();
				i += 2;
			}
		} else {
			rawdata[i] = (*itG)->get_ID();
			rawdata[i+1] = (*itG)->get_oldID();
			i += 2;
		}
	}
	if ( i != 2*NrGrains ) {
		cout << "Determination of parenthood resulted in inconsistencies!" << endl; delete [] rawdata; return;
	}

//##MK::Test, can be deleted for ( i = 0; i < 2*NrGrains; i += 2 ) { cout << "SGID/GID = " << rawdata[i] << "\t\t" << rawdata[i+1] << endl; }

	//store in binary file
	stringstream parentfn;
	parentfn << "Parenthood.SimID." << Settings::SimID << ".bin";

	FILE* binaryFile;
	binaryFile = fopen(parentfn.str().c_str(), "wb"); //##MK::was "w"
	fwrite( rawdata , sizeof(unsigned int), 2*NrGrains, binaryFile);
	fclose(binaryFile);

	delete [] rawdata,

	myprofiler.logev("WritingBINARYParenthood", (omp_get_wtime() - ptime));
	cout << "Storing of parenthood in binary file was successful." << endl;
}


void microStructureHdl::DebugGetDistParentGrainBnd( const string udsfn, const string rawfn )
{
	//load udsfn file

	//load rawfn file
}


void microStructureHdl::DebugHDF5() {
		//https:\//support.hdfgroup.org/HDF5/examples/api18-c.html
	hid_t file, space, dset, group, dcpl;
	herr_t status;

	//create a new file using the default properties
	string fname = "Debug.h5";
std::cout << "Debugging HDF5 I/O " << fname << std::endl;

	file = H5Fcreate( fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//create DAMASK dataset groups in the file
	//for further references see M. Diehl, P. Eisenlohr, C. Zhang, J. Nastola, P. Shanthraj and F. Roters:
	//"A flexible and efficient output file format for grain scale multiphysics simulations"
	//Integrating Materials and Manufacturing Innovation 2016

	//generate group
	group = H5Gcreate2( file, "/geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
	status = H5Gclose(group);

	group = H5Gcreate2( file, "/mapping", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
	status = H5Gclose(group);

	group = H5Gcreate2( file, "/mapping/cells", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
	status = H5Gclose(group);

	group = H5Gcreate2( file, "/mapping/cells/constitutive", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
	status = H5Gclose(group);

	group = H5Gcreate2( file, "/mapping/cells/constitutive/plasticity", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
	status = H5Gclose( group );

	/*for ( int slipsystem = 0; slipsystem < 1; slipsystem++ ) {
		string grpname = "/mapping/cells/constitutive/plasticity/" + std::to_string(slipsystem);

		//create the dimensions, space of the dataset and copy data from the SU volume
		hsize_t dims[2] = { this->myCAGeometry.nboxvol_rdndtd, 1 };

		float* rho = NULL;
		rho = new float[this->myCAGeometry.nboxvol_rdndtd];
		QUICKASSERT( rho != NULL );
		h5_return_dislocationdensity( rho, this->myCAGeometry.nboxvol_rdndtd );

		space = H5Screate_simple(2, dims, NULL);

		//create the dataset with default properties
		dset = H5Dcreate(file, grpname.c_str(), H5T_IEEE_F32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		//write dataset
		status = H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho );

		status = H5Sclose (space);
		status = H5Dclose (dset);

		delete [] rho; rho = NULL;
	}*/

	status = H5Fclose(file);
}




void microStructureHdl::DebugHDF5XDMF()
{
	hid_t file, space, dset, group, dcpl;
	herr_t status;

	//create a new file using the default properties
	string fname = "DebugXDMFAuto.h5";
std::cout << "Debugging HDF5 I/O with XDMF " << fname << std::endl;

	stringstream xdmf_fname;
	xdmf_fname << "DebugXDMFTimeAuto.xdmf";
	ofstream xdmf;
	xdmf.open ( xdmf_fname.str().c_str() );


	//XDMF header
	// __ / __ ascii code 47
	// __ " __ ascii code 34
	xdmf <<(char)60<<(char)63<<"xml version"<<(char)61<<(char)34<<"1.0"<<(char)34<< " "<<(char)63<<">\n";
	xdmf <<(char)60<<(char)33<<"DOCTYPE Xdmf SYSTEM "<<(char)34<<"Xdmf.dtd"<<(char)34<<" []>\n";
	xdmf <<(char)60<<"Xdmf xmlns:xi="<<(char)34<<"http:"<<(char)47<<(char)47<<"www.w3.org/2003/XInclude"<<(char)34<<" Version"<<(char)61<<(char)34<<"2.2"<<(char)34<<">\n";
	xdmf <<"\t<Domain>\n";


	file = H5Fcreate( fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//create DAMASK dataset groups in the file
	//for further references see M. Diehl, P. Eisenlohr, C. Zhang, J. Nastola, P. Shanthraj and F. Roters:
	//"A flexible and efficient output file format for grain scale multiphysics simulations"
	//Integrating Materials and Manufacturing Innovation 2016

	//dummy coordinates
#define DIM0	((20)*(20)*(20))
#define DIM1	3
	double p[DIM0][DIM1];
	unsigned int j = 0;
	for ( unsigned int z = 0; z < 20; z++ ) {
		for ( unsigned int y = 0; y < 20; y++ ) {
			for ( unsigned int x = 0; x < 20; x++ ) {
				p[j][0] = x;	p[j][1] = y;		p[j][2] = z;
				j++;
			}
		}
	}

	//once write the initial coordinates
	group = H5Gcreate2( file, "/inc0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
	status = H5Gclose(group);

	hsize_t dims[2] = {DIM0, DIM1};
	space = H5Screate_simple (2, dims, NULL);
	dset = H5Dcreate(file, "/inc0/coordinates", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p);
	status = H5Dclose (dset);
	status = H5Sclose (space);



	//now pretend a history of dislocation density evolution
	//initial values
	double* test = new double[DIM0];
	for ( unsigned int j = 0; j < DIM0; j++ ) {
		if ( j < 4000 )
			test[j] = 1.0e12;
		else
			test[j] = 2.0e12;
	}

	//exemplify a time evolution
	string grpname;
	string grpname1, grpname2;

	for ( unsigned int t = 1; t <= 100; t++ ) {
		//"strain..."
		for ( unsigned int j = 0; j < DIM0; j++ )
			test[j] *= 1.02;

		//generate group
		grpname = "/inc" + std::to_string(t);
		group = H5Gcreate2( file, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
		status = H5Gclose(group);

		//adjust dislocation density artificially for "testing"...
		grpname = grpname + "/constitutive";
		group = H5Gcreate2( file, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
		status = H5Gclose(group);
		grpname = grpname + "/plasticity";
		group = H5Gcreate2( file, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
		status = H5Gclose(group);
		grpname = grpname + "/nonlocal";
		group = H5Gcreate2( file, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
		status = H5Gclose(group);

		hsize_t tdims[1] = {DIM0};
		space = H5Screate_simple (1, tdims, NULL);
		grpname1 = grpname + "/rho1";
		dset = H5Dcreate(file, grpname1.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, test);
		status = H5Dclose (dset);
		status = H5Sclose (space);

		space = H5Screate_simple (1, tdims, NULL);
		grpname2 = grpname + "/rho2";
		dset = H5Dcreate(file, grpname2.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, test);
		status = H5Dclose (dset);
		status = H5Sclose (space);

		/*space = H5Screate_simple (1, tdims, NULL);
		dset = H5Dcreate(file, "/inc1/constitutive/plasticity/nonlocal/rho2", H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, test);
		status = H5Dclose (dset);
		status = H5Sclose (space);*/

		//UPDATE XDMF
		xdmf <<"\t\t<Grid Name=" <<(char)34<<"DAMASK_Results"<<(char)34<<" GridType="<<(char)34<<"Uniform"<<(char)34<<">\n";
		xdmf <<"\t\t\t<Time Type="<<(char)34<<"Single"<<(char)34<<" Value="<<(char)34<<"   "<< setprecision(6) << (double) t<<(char)34<<" />\n";
		xdmf <<"\t\t\t<Topology TopologyType="<<(char)34<<"Polyvertex"<<(char)34<<" Dimensions="<<(char)34<<"20 20 20"<<(char)34<<"/>\n";
		xdmf <<"\t\t\t<Geometry GeometryType="<<(char)34<<"XYZ"<<(char)34<<">\n";
		xdmf <<"\t\t\t\t<DataItem   \n";
		xdmf <<"\t\t\t\t\tFormat="<<(char)34<<"HDF"<<(char)34<<"\n";
		xdmf <<"\t\t\t\t\tNumberType="<<(char)34<<"Float"<<(char)34<<"\n";
		xdmf <<"\t\t\t\t\tPrecision="<<(char)34<<"8"<<(char)34<<"\n";
		xdmf <<"\t\t\t\t\tDimensions="<<(char)34<<"8000 3"<<(char)34<<">\n";
		xdmf <<"\t\t\t\t\tDebugXDMFTime.h5:/inc0/coordinates\n";
		xdmf <<"\t\t\t\t</DataItem>\n";
		xdmf <<"\t\t\t</Geometry>\n";

		xdmf <<"\t\t<Attribute Name="<<(char)34<<"DislocationDensity1"<<(char)34<<" AttributeType="<<(char)34<<"Scalar"<<(char)34<<" Center="<<(char)34<<"Node"<<(char)34<<">\n";
		xdmf <<"\t\t\t<DataItem   \n";
		xdmf <<"\t\t\t\tFormat="<<(char)34<<"HDF"<<(char)34<<"\n";
		xdmf <<"\t\t\t\tNumberType="<<(char)34<<"Float"<<(char)34<<"\n";
		xdmf <<"\t\t\t\tPrecision="<<(char)34<<"8"<<(char)34<<"\n";
		xdmf <<"\t\t\t\tDimensions="<<(char)34<<"8000"<<(char)34<<">\n";
		xdmf <<"\t\t\t\tDebugXDMFTime.h5:/inc"<<(int)t<<"/constitutive/plasticity/nonlocal/rho1\n";
		xdmf <<"\t\t\t\t</DataItem>\n";
		xdmf <<"\t\t\t</Attribute>\n";

		xdmf <<"\t\t<Attribute Name="<<(char)34<<"DislocationDensity2"<<(char)34<<" AttributeType="<<(char)34<<"Scalar"<<(char)34<<" Center="<<(char)34<<"Node"<<(char)34<<">\n";
		xdmf <<"\t\t\t<DataItem   \n";
		xdmf <<"\t\t\t\tFormat="<<(char)34<<"HDF"<<(char)34<<"\n";
		xdmf <<"\t\t\t\tNumberType="<<(char)34<<"Float"<<(char)34<<"\n";
		xdmf <<"\t\t\t\tPrecision="<<(char)34<<"8"<<(char)34<<"\n";
		xdmf <<"\t\t\t\tDimensions="<<(char)34<<"8000"<<(char)34<<">\n";
		xdmf <<"\t\t\t\tDebugXDMFTime.h5:/inc"<<(int)t<<"/constitutive/plasticity/nonlocal/rho2\n";
		xdmf <<"\t\t\t\t</DataItem>\n";
		xdmf <<"\t\t\t</Attribute>\n";
		xdmf <<"\t\t</Grid>"<<endl;

		cout << "Timestep " << t << endl;
	}

	//Finalize XDMF
	xdmf << "\t</Domain>\n";
	xdmf << "</Xdmf>";
	xdmf.flush();
	xdmf.close();

	delete [] test;
	status = H5Fclose (file);
}

//void microStructureHdl::SaveHDF5() {}



void microStructureHdl::SaveHDF5_CollectMetadata( SubgrainMetadata* buf ) {
	//fill existing buffer of known size
	unsigned int i = 0;
	vector<Grains*>::iterator itG;
	for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
		if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {
			vector<SubGrain*>::iterator it;
			myQuaternion parentori = (*itG)->getOri();
			double* parentbunge = parentori.Quaternion2Euler();
			double oriscatter = (*itG)->get_OriScatterFromPrefOri();
			double see = (*itG)->getSEE();
			double seegrscatter = (*itG)->get_SEEGrainScatterFromPrefOri();
			double seesgrscatter = (*itG)->get_SEESubgrainScatterFromPrefOri();
			double sizescaler = (*itG)->get_RelSizeScalingFromPrefori();

			for (it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++) {
				buf[i].gid = (*it)->get_ID();
				buf[i].prid = (*itG)->get_oldID();
				buf[i].x = (((*it)->getMaxX() - (*it)->getMinX()) / 2 + (*it)->getMinX());
				buf[i].y = (((*it)->getMaxY() - (*it)->getMinY()) / 2 + (*it)->getMinY());
				buf[i].z = (((*it)->getMaxZ() - (*it)->getMinZ()) / 2 + (*it)->getMinZ());
				buf[i].xmi = (*it)->getMinX();
				buf[i].xmx = (*it)->getMaxX();
				buf[i].ymi = (*it)->getMinY();
				buf[i].ymx = (*it)->getMaxY();
				buf[i].zmi = (*it)->getMinZ();
				buf[i].zmx = (*it)->getMaxZ();

				myQuaternion* ori = (*it)->get_Orientation();
				double* bunge = ori->Quaternion2Euler();
				buf[i].bunge1 = bunge[0];
				buf[i].bunge2 = bunge[1];
				buf[i].bunge3 = bunge[2];
				delete [] bunge;

				buf[i].size = (*it)->get_Volume(); //##MK::not corrected for periodicity!
				buf[i].see = (*it)->get_SEE();
				buf[i].psee = see;
				buf[i].d2pr_meas = (double) parentori.misorientationCubicQxQ( ori );
				buf[i].d2pr_targ = oriscatter;
				buf[i].psee_grsc = seegrscatter;
				buf[i].psee_sgrsc = seesgrscatter;
				buf[i].pszsc = sizescaler;
				buf[i].pbunge1 = parentbunge[0];
				buf[i].pbunge2 = parentbunge[1];
				buf[i].pbunge3 = parentbunge[2];

				i++;
			}
			delete[] parentbunge;
		}
		else {
			buf[i].gid = (*itG)->get_ID();
			buf[i].prid = (*itG)->get_oldID();
			buf[i].x = (((*itG)->getMaxX() - (*itG)->getMinX()) / 2 + (*itG)->getMinX());
			buf[i].y = (((*itG)->getMaxY() - (*itG)->getMinY()) / 2 + (*itG)->getMinY());
			buf[i].z = (((*itG)->getMaxZ() - (*itG)->getMinZ()) / 2 + (*itG)->getMinZ());
			buf[i].xmi = (*itG)->getMinX();
			buf[i].xmx = (*itG)->getMaxX();
			buf[i].ymi = (*itG)->getMinY();
			buf[i].ymx = (*itG)->getMaxY();
			buf[i].zmi = (*itG)->getMinZ();
			buf[i].zmx = (*itG)->getMaxZ();

			myQuaternion ori = (*itG)->getOri();
			double* bunge = ori.Quaternion2Euler();
			buf[i].bunge1 = bunge[0];
			buf[i].bunge2 = bunge[1];
			buf[i].bunge3 = bunge[2];
			buf[i].pbunge1 = bunge[0];
			buf[i].pbunge2 = bunge[1];
			buf[i].pbunge3 = bunge[2];
			delete [] bunge;

			buf[i].size = (*itG)->get_Volume(); //##MK::not corrected for periodicity!
			buf[i].see = (*itG)->getSEE();
			buf[i].psee = (*itG)->getSEE();
			buf[i].d2pr_meas = 0.0;
			buf[i].d2pr_targ = (*itG)->get_OriScatterFromPrefOri();
			buf[i].psee_grsc = (*itG)->get_SEEGrainScatterFromPrefOri();
			buf[i].psee_sgrsc = (*itG)->get_SEESubgrainScatterFromPrefOri();
			buf[i].pszsc = (*itG)->get_RelSizeScalingFromPrefori();
			i++;
		}
	}
}


void microStructureHdl::SaveHDF5_WriteGeometry( const char* fname ) {
	//handles
	hid_t			file, filetype, memtype, space, dset, group;
	herr_t			status;
	H5D_layout_t	layout;
	hsize_t			dims[1] = { 3 };

	SimMetadata* buffer = NULL;
	buffer = new SimMetadata[1];

	//##MK::implicit conversion from uint to int will become problematic if populations become too large..., i.e. 2^32 -1
	buffer[0].simid = Settings::SimID;
	buffer[0].NGrains = CountNumberOfSubgrains();
	buffer[0].firstgid = this->get_first_id();
	buffer[0].DX = Settings::NumberOfPointsPerSubGrain;
	buffer[0].DY = Settings::NumberOfPointsPerSubGrain;
	buffer[0].DZ = 1;
	buffer[0].NX = Settings::NumberOfGridpoints;
	buffer[0].NY = Settings::NumberOfGridpoints;
	buffer[0].NZ = 1;
	if ( Settings::PlotDimension == E_3D ) {
		buffer[0].DZ = Settings::NumberOfPointsPerSubGrain;
		buffer[0].NZ = Settings::NumberOfGridpoints;
	}

	file = H5Fopen( fname, H5F_ACC_RDWR, H5P_DEFAULT);

	//create first a group sim_meta
	string grpname = "sim_meta";
	group = H5Gcreate(file, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Gclose(group);

	grpname = "sim_meta/population";
	space = H5Screate_simple (1, dims, NULL );
	dset = H5Dcreate (file, grpname.c_str(), H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//SimID; NGrains, FirstID
	int data[3] = { buffer[0].simid, buffer[0].NGrains, buffer[0].firstgid };
	status = H5Dwrite (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0] );
	status = H5Dclose (dset);
	status = H5Sclose (space);

	grpname = "sim_meta/domainsize";
	space = H5Screate_simple (1, dims, NULL );
	dset = H5Dcreate (file, grpname.c_str(), H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	//NX, NY, NZ
	data[0] = buffer[0].NX;
	data[1] = buffer[0].NY;
	data[2] = buffer[0].NZ;
	status = H5Dwrite (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data[0] );
	status = H5Dclose (dset);
	status = H5Sclose (space);

	status = H5Fclose (file);

	delete [] buffer;
}


void microStructureHdl::SaveHDF5_WriteGeometryCompound( const char* fname ) {
	//handles
	hid_t		file, filetype, memtype, space, dset;
	herr_t		status;
	hsize_t		dims[1] = { 1 };

	SimMetadata* buffer = NULL;
	buffer = new SimMetadata[1];

	//##MK::implicit conversion from uint to int will become problematic if populations become too large..., i.e. 2^32 -1
	buffer[0].simid = Settings::SimID;
	buffer[0].NGrains = CountNumberOfSubgrains();
	buffer[0].firstgid = this->get_first_id();
	buffer[0].DX = Settings::NumberOfPointsPerSubGrain;
	buffer[0].DY = Settings::NumberOfPointsPerSubGrain;
	buffer[0].DZ = 1;
	buffer[0].NX = Settings::NumberOfGridpoints;
	buffer[0].NY = Settings::NumberOfGridpoints;
	buffer[0].NZ = 1;
	if ( Settings::PlotDimension == E_3D ) {
		buffer[0].DZ = Settings::NumberOfPointsPerSubGrain;
		buffer[0].NZ = Settings::NumberOfGridpoints;
	}

	file = H5Fopen( fname, H5F_ACC_RDWR, H5P_DEFAULT);

	//Create variable-length string datatype//strtype = H5Tcopy (H5T_C_S1);//status = H5Tset_size (strtype, H5T_VARIABLE);
	//Create the compound datatype for memory
	memtype		= H5Tcreate (H5T_COMPOUND, sizeof (SimMetadata));

	status		= H5Tinsert (memtype, "SimulationID", HOFFSET (SimMetadata, simid), H5T_NATIVE_INT); //H5T_NATIVE_INT);
	status		= H5Tinsert (memtype, "NumberOfSubgrains", HOFFSET (SimMetadata, NGrains), H5T_NATIVE_INT); //H5T_NATIVE_UINT);
	status		= H5Tinsert (memtype, "FirstGrainID", HOFFSET (SimMetadata, firstgid), H5T_NATIVE_INT); //H5T_NATIVE_UINT);

	status		= H5Tinsert (memtype, "GridPtsPerSubGrX", HOFFSET (SimMetadata, DX), H5T_NATIVE_INT); //H5T_NATIVE_UINT);
	status		= H5Tinsert (memtype, "GridPtsPerSubGrY", HOFFSET (SimMetadata, DY), H5T_NATIVE_INT); //H5T_NATIVE_UINT);
	status		= H5Tinsert (memtype, "GridPtsPerSubGrZ", HOFFSET (SimMetadata, DZ), H5T_NATIVE_INT); //H5T_NATIVE_UINT);

	status		= H5Tinsert (memtype, "GridPtsDomainX", HOFFSET (SimMetadata, NX), H5T_NATIVE_INT); //H5T_NATIVE_UINT);
	status		= H5Tinsert (memtype, "GridPtsDomainY", HOFFSET (SimMetadata, NY), H5T_NATIVE_INT); //H5T_NATIVE_UINT);
	status		= H5Tinsert (memtype, "GridPtsDomainZ", HOFFSET (SimMetadata, NZ), H5T_NATIVE_INT); //H5T_NATIVE_UINT);

	//Create the compound datatype for the file with manually calculated offset for each member, ##MK::sizeof (hvl_t) for variable length string types...
	filetype = H5Tcreate (H5T_COMPOUND, (10*4) );

	size_t offset = 0*4;
	status = H5Tinsert (filetype, "SimulationID", offset, H5T_STD_I32LE);				offset = 1*4;
	status = H5Tinsert (filetype, "NumberOfSubgrains", offset, H5T_STD_I32LE);			offset = 2*4;
	status = H5Tinsert (filetype, "FirstGrainID", offset, H5T_STD_I32LE);				offset = 3*4;
	status = H5Tinsert (filetype, "GridPtsPerSubGrX", offset, H5T_STD_I32LE);			offset = 4*4;
	status = H5Tinsert (filetype, "GridPtsPerSubGrY", offset, H5T_STD_I32LE);			offset = 5*4;
	status = H5Tinsert (filetype, "GridPtsPerSubGrZ", offset, H5T_STD_I32LE);			offset = 6*4;
	status = H5Tinsert (filetype, "GridPtsDomainX", offset, H5T_STD_I32LE);				offset = 7*4;
	status = H5Tinsert (filetype, "GridPtsDomainY", offset, H5T_STD_I32LE);				offset = 8*4;
	status = H5Tinsert (filetype, "GridPtsDomainZ", offset, H5T_STD_I32LE);				offset = 9*4;

	//Create dataspace.  Setting maximum size to NULL sets the maximum, size to be the current size
	space = H5Screate_simple (1, dims, NULL);

	//Create the dataset and write the compound data to it.
	dset = H5Dcreate (file, "sim_meta", filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

	//close resources
	status = H5Dclose (dset);
	status = H5Sclose (space);
	status = H5Tclose (filetype);
	status = H5Fclose (file);

	delete [] buffer; buffer = NULL;
}


void microStructureHdl::SaveHDF5_WriteSubgrainMetadataType( const char* fname )
{
	hid_t		file, space, dset, group;
	herr_t		status;

	size_t len = this->CountNumberOfSubgrains();
	hsize_t		dims[1] = { len };

	SubgrainMetadata* buffer = NULL;
	buffer = new SubgrainMetadata[len];

	SaveHDF5_CollectMetadata( buffer );

	file = H5Fopen( fname, H5F_ACC_RDWR, H5P_DEFAULT);

	//create first a group sim_meta
	string grpname = "grains_meta";
	group = H5Gcreate(file, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Gclose(group);

	grpname = "grains_meta/rho";
	space = H5Screate_simple (1, dims, NULL );
	dset = H5Dcreate (file, grpname.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	double* val = NULL;
	val = new double[len];
	for ( unsigned int g = 0; g < len; g++ )
		val[g] = buffer[g].see;

	status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, val );
	status = H5Dclose (dset);
	status = H5Sclose (space);
	delete [] val;

	status = H5Fclose (file);

	delete [] buffer;
}


void microStructureHdl::SaveHDF5_WriteSubgrainMetadataTypeCompound( const char* fname )
{
	//handles
	hid_t		file, filetype, memtype, space, dset, dcpl;
	herr_t		status;

	size_t len = this->CountNumberOfSubgrains();
	hsize_t		dims[1] = { len };

	SubgrainMetadata* buffer = NULL;
	buffer = new SubgrainMetadata[len];

	SaveHDF5_CollectMetadata( buffer );

	file = H5Fopen( fname, H5F_ACC_RDWR, H5P_DEFAULT);

	//Create variable-length string datatype//strtype = H5Tcopy (H5T_C_S1);//status = H5Tset_size (strtype, H5T_VARIABLE);
	//Create the compound datatype for memory
	memtype		= H5Tcreate (H5T_COMPOUND, sizeof (SubgrainMetadata));
	status		= H5Tinsert (memtype, "Subgrain ID", HOFFSET (SubgrainMetadata, gid), H5T_NATIVE_INT);
	status		= H5Tinsert (memtype, "Parentgrain ID", HOFFSET (SubgrainMetadata, prid), H5T_NATIVE_INT);
	status		= H5Tinsert (memtype, "Barycenter X", HOFFSET (SubgrainMetadata, x), H5T_NATIVE_INT);
	status		= H5Tinsert (memtype, "Barycenter Y", HOFFSET (SubgrainMetadata, y), H5T_NATIVE_INT);

	status		= H5Tinsert (memtype, "Barycenter Z", HOFFSET (SubgrainMetadata, z), H5T_NATIVE_INT);
	status		= H5Tinsert (memtype, "AABB XMI", HOFFSET (SubgrainMetadata, xmi), H5T_NATIVE_INT);
	status		= H5Tinsert (memtype, "AABB XMX", HOFFSET (SubgrainMetadata, xmx), H5T_NATIVE_INT);
	status		= H5Tinsert (memtype, "AABB YMI", HOFFSET (SubgrainMetadata, ymi), H5T_NATIVE_INT);

	status		= H5Tinsert (memtype, "AABB YMX", HOFFSET (SubgrainMetadata, ymx), H5T_NATIVE_INT);
	status		= H5Tinsert (memtype, "AABB ZMI", HOFFSET (SubgrainMetadata, zmi), H5T_NATIVE_INT);
	status		= H5Tinsert (memtype, "AABB ZMX", HOFFSET (SubgrainMetadata, zmx), H5T_NATIVE_INT);

	status		= H5Tinsert (memtype, "Bunge phi1 (rad)", HOFFSET (SubgrainMetadata, bunge1), H5T_NATIVE_DOUBLE);
	status		= H5Tinsert (memtype, "Bunge Psi (rad)", HOFFSET (SubgrainMetadata, bunge2), H5T_NATIVE_DOUBLE);
	status		= H5Tinsert (memtype, "Bunge phi2 (rad)", HOFFSET (SubgrainMetadata, bunge3), H5T_NATIVE_DOUBLE);
	status		= H5Tinsert (memtype, "Size", HOFFSET (SubgrainMetadata, size), H5T_NATIVE_DOUBLE);

	status		= H5Tinsert (memtype, "SEE (unit??)", HOFFSET (SubgrainMetadata, see), H5T_NATIVE_DOUBLE);
	status		= H5Tinsert (memtype, "Parent SEE (unit??)", HOFFSET (SubgrainMetadata, psee), H5T_NATIVE_DOUBLE);
	status		= H5Tinsert (memtype, "Disori2Parent measured (rad)", HOFFSET (SubgrainMetadata, d2pr_meas), H5T_NATIVE_DOUBLE);

	status		= H5Tinsert (memtype, "Disori2Parent target (rad)", HOFFSET (SubgrainMetadata, d2pr_targ), H5T_NATIVE_DOUBLE);
	status		= H5Tinsert (memtype, "Parent SEE grain scatter (unit??)", HOFFSET (SubgrainMetadata, psee_grsc), H5T_NATIVE_DOUBLE);
	status		= H5Tinsert (memtype, "Parent SEE subgrain scatter (unit??)", HOFFSET (SubgrainMetadata, psee_sgrsc), H5T_NATIVE_DOUBLE);
	status		= H5Tinsert (memtype, "Parent size scaler", HOFFSET (SubgrainMetadata, pszsc), H5T_NATIVE_DOUBLE);

	status		= H5Tinsert (memtype, "Parent phi1 (rad)", HOFFSET (SubgrainMetadata, pbunge1), H5T_NATIVE_DOUBLE);
	status		= H5Tinsert (memtype, "Parent Psi (rad)", HOFFSET (SubgrainMetadata, pbunge2), H5T_NATIVE_DOUBLE);
	status		= H5Tinsert (memtype, "Parent phi2 (rad)", HOFFSET (SubgrainMetadata, pbunge3), H5T_NATIVE_DOUBLE);

	//Create the compound datatype for the file with manually calculated offset for each member, ##MK::sizeof (hvl_t) for variable length string types...
	filetype = H5Tcreate (H5T_COMPOUND, (11*4 + 14*8) );

	size_t offset = 0*4 + 0*8;
	status = H5Tinsert (filetype, "Subgrain ID", offset, H5T_STD_I32LE);							offset = 1*4 + 0*8;
	status = H5Tinsert (filetype, "Parentgrain ID", offset, H5T_STD_I32LE);							offset = 2*4 + 0*8;
	status = H5Tinsert (filetype, "Barycenter X", offset, H5T_STD_I32LE);							offset = 3*4 + 0*8;
	status = H5Tinsert (filetype, "Barycenter Y", offset, H5T_STD_I32LE);							offset = 4*4 + 0*8;
	status = H5Tinsert (filetype, "Barycenter Z", offset, H5T_STD_I32LE);							offset = 5*4 + 0*8;
	status = H5Tinsert (filetype, "AABB XMI", offset, H5T_STD_I32LE);								offset = 6*4 + 0*8;
	status = H5Tinsert (filetype, "AABB XMX", offset, H5T_STD_I32LE);								offset = 7*4 + 0*8;
	status = H5Tinsert (filetype, "AABB YMI", offset, H5T_STD_I32LE);								offset = 8*4 + 0*8;
	status = H5Tinsert (filetype, "AABB YMX", offset, H5T_STD_I32LE);								offset = 9*4 + 0*8;
	status = H5Tinsert (filetype, "AABB ZMI", offset, H5T_STD_I32LE);								offset = 10*4 + 0*8;
	status = H5Tinsert (filetype, "AABB ZMX", offset, H5T_STD_I32LE);								offset = 11*4 + 0*8;

	status = H5Tinsert (filetype, "Bunge phi1 (rad)", offset, H5T_IEEE_F64LE);						offset = 11*4 + 1*8;
	status = H5Tinsert (filetype, "Bunge Psi (rad)", offset, H5T_IEEE_F64LE);						offset = 11*4 + 2*8;
	status = H5Tinsert (filetype, "Bunge phi2 (rad)", offset, H5T_IEEE_F64LE);						offset = 11*4 + 3*8;
	status = H5Tinsert (filetype, "Size", offset, H5T_IEEE_F64LE);									offset = 11*4 + 4*8;
	status = H5Tinsert (filetype, "SEE (unit??)", offset, H5T_IEEE_F64LE);							offset = 11*4 + 5*8;
	status = H5Tinsert (filetype, "Parent SEE (unit??)", offset, H5T_IEEE_F64LE);					offset = 11*4 + 6*8;
	status = H5Tinsert (filetype, "Disori2Parent measured (rad)", offset, H5T_IEEE_F64LE);			offset = 11*4 + 7*8;
	status = H5Tinsert (filetype, "Disori2Parent target (rad)", offset, H5T_IEEE_F64LE);			offset = 11*4 + 8*8;
	status = H5Tinsert (filetype, "Parent SEE grain scatter (unit??)", offset, H5T_IEEE_F64LE);		offset = 11*4 + 9*8;
	status = H5Tinsert (filetype, "Parent SEE subgrain scatter (unit??)", offset, H5T_IEEE_F64LE);	offset = 11*4 + 10*8;
	status = H5Tinsert (filetype, "Parent size scaler", offset, H5T_IEEE_F64LE);					offset = 11*4 + 11*8;
	status = H5Tinsert (filetype, "Parent phi1 (rad)", offset, H5T_IEEE_F64LE);						offset = 11*4 + 12*8;
	status = H5Tinsert (filetype, "Parent Psi (rad)", offset, H5T_IEEE_F64LE);						offset = 11*4 + 13*8;
	status = H5Tinsert (filetype, "Parent phi2 (rad)", offset, H5T_IEEE_F64LE);

	//Create dataspace.  Setting maximum size to NULL sets the maximum, size to be the current size
	space = H5Screate_simple (1, dims, NULL);

	//Create the dataset and write the compound data to it.
	dset = H5Dcreate (file, "grains_meta", filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

	//##MK::take care about too large datasets...

	//close resources
	status = H5Dclose (dset);
	status = H5Sclose (space);
	status = H5Tclose (filetype);
	status = H5Fclose (file);

	delete [] buffer; buffer = NULL;
}


void microStructureHdl::SaveHDF5_WriteVoxeldata( const char* fname ) {
	//add voxel data as a compact data set
	hid_t file, space, dset;
	herr_t status;
	H5D_layout_t	layout;

	unsigned int len = Settings::NumberOfGridpoints * Settings::NumberOfGridpoints * 1;
	if ( Settings::PlotDimension == E_3D )	len *= Settings::NumberOfGridpoints;

	//hsize_t		dims[2] = { 1, len }; //##MK::for along columns
	hsize_t dims[1] = { len }; //along rows

	file = H5Fopen( fname, H5F_ACC_RDWR, H5P_DEFAULT);

	//create dataspace
	//space = H5Screate_simple (2, daims, NULL); //##MK::along the columns

	space = H5Screate_simple (1, dims, NULL ); //##MK::along the rows

	//create the dataset with default props
	dset = H5Dcreate (file, "volume_grassgn", H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//##MK::be careful here at the moment implict cast fomr uint to int....!
	status = H5Dwrite (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_container->getRawData() );

	status = H5Dclose (dset);
	status = H5Sclose (space);
	status = H5Fclose (file);
}

void microStructureHdl::SaveHDF5( void )
{
	if ( Settings::IO_HDF5 == false )
		return;

	double h5time0 = omp_get_wtime();
	double h5time1 = h5time0;

	stringstream h5fname;
	h5fname << "Microstructure.SimID." << Settings::SimID << ".h5";

	hid_t file; //generate a file such that it exists
	herr_t status;
	file = H5Fcreate ( h5fname.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Fclose( file );

	SaveHDF5_WriteGeometry( h5fname.str().c_str() );
	//SaveHDF5_WriteGeometryCompound( h5fname.str().c_str() );
	h5time1 = omp_get_wtime();
	myprofiler.logev("HDF5WriteGeometry", ( h5time1 - h5time0 ));
	h5time0 = h5time1;

	SaveHDF5_WriteSubgrainMetadataType( h5fname.str().c_str() );
	//SaveHDF5_WriteSubgrainMetadataTypeCompound( h5fname.str().c_str() );
	h5time1 = omp_get_wtime();
	myprofiler.logev("HDF5WriteMetadata", ( h5time1 - h5time0 ));
	h5time0 = h5time1;

	SaveHDF5_WriteVoxeldata( h5fname.str().c_str() );
	h5time1 = omp_get_wtime();
	myprofiler.logev("HDF5AssignmentWritten", ( h5time1 - h5time0 ));

	cout << "HDF5 file " << h5fname.str().c_str() << " written to file!" << endl;
}



void microStructureHdl::CreateColormap( struct RGB* thecolormap )
{
	double gtime = omp_get_wtime();

#pragma omp parallel
	{
		vector<Grains*>::iterator itG;
		int gid = 0;
		for (itG = ++m_grains.begin(); itG != m_grains.end(); itG++) {
			if ( gid % omp_get_num_threads() == omp_get_thread_num() ) {
				if ((Settings::NumberOfSubgrains != 0) && (*itG)->m_SubGrains.size() > 1) {
					vector<SubGrain*>::iterator it;
					for (it = ++((*itG)->m_SubGrains.begin()); it != (*itG)->m_SubGrains.end(); it++) {
						unsigned int id = (*it)->get_ID() - 1;
						myQuaternion* ori = (*it)->get_Orientation();

						unsigned char argb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN, UCHAR_RANGE_MIN };
						ori->quat2ipfz( argb );
//#pragma omp critical
//{
//	std::cout << "thread/sid/R/G/B/colormap[id].ALPHA = " << omp_get_thread_num() << ";" << "\t\t" << id << "\t\t" << (int) argb[0] << ";" << (int) argb[1] << ";" << (int) argb[2] << "\t\t\t" << ori->get_x() << ";" << ori->get_y() << std::endl;
//}

						thecolormap[id].R = argb[REDCHAN];
						thecolormap[id].G = argb[GREENCHAN];
						thecolormap[id].B = argb[BLUECHAN];
						thecolormap[id].ALPHA = UCHAR_RANGE_MAX; //=255;


					} //for all subgrains
				}
				else {
					unsigned int id = (*itG)->get_ID()-1;
					myQuaternion ori = (*itG)->getOri();
					unsigned char argb[3] = { UCHAR_RANGE_MIN, UCHAR_RANGE_MIN, UCHAR_RANGE_MIN };

					ori.quat2ipfz( argb );

					thecolormap[id].R = argb[REDCHAN];
					thecolormap[id].G = argb[GREENCHAN];
					thecolormap[id].B = argb[BLUECHAN];
					thecolormap[id].ALPHA = UCHAR_RANGE_MAX; //=255;
//#pragma omp critical
//{
//std::cout << "thread/gid/R/G/B/colormap[id].ALPHA = " << omp_get_thread_num() << ";" << "\t\t" << id << "\t\t" << (int) argb[0] << ";" << (int) argb[1] << ";" << (int) argb[2] << "\t\t\t" << ori.get_x() << ";" << ori.get_y() << std::endl;
//}
				}
			} //end the work I do
			gid++;
		}
	} //end of parallel region

	myprofiler.logev("ParallelIPFColorMapping", (omp_get_wtime() - gtime));
}


unsigned int microStructureHdl::CountNumberOfSubgrains() {
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
cout << "Total number of subgrains = " << NumberOfSubgrains << endl;
	return NumberOfSubgrains;
}


void microStructureHdl::PlotIPF2DSection() {
	//MK::requires to become executed after an ID rehashing for the subgrains has been performed!

	if ( Settings::PlotIPF2DSection == true ) {
		double gtime = omp_get_wtime();

		double edgelength = Settings::NumberOfGridpoints;
		unsigned int xmi = Settings::PlotWindowXMin * edgelength;
		unsigned int xmx = Settings::PlotWindowXMax * edgelength;
		unsigned int ymi = Settings::PlotWindowYMin * edgelength;
		unsigned int ymx = Settings::PlotWindowYMax * edgelength;
		if ( xmx >= m_container->getMaxX() ) xmx--; //bounds check
		if ( ymx >= m_container->getMaxY() ) ymx--;
		if ( (xmx - xmi + 1) > IPFZMAPPING_MAXSIZE ) {
			xmi = 0;
			xmx = IPFZMAPPING_MAXSIZE;
			ymi = 0;
			ymx = IPFZMAPPING_MAXSIZE;
		}
		unsigned int z = 0;
		unsigned int xlen = xmx - xmi + 1;

		cout << "Generating a IPF mapping of the domain in 2D on " << xmi << ";" << xmx << ";" << ymi << ";" << ymx << " width " << (xmx - xmi + 1) << " x " << (ymx - ymi + 1) << endl;

		unsigned int SubgrainsInTotal = CountNumberOfSubgrains();
		struct RGB* colormap = new struct RGB[SubgrainsInTotal];

//for ( unsigned int id = 0; id < SubgrainsInTotal; id++ ) {cout << "thread/0/R/G/B/colormap[id].ALPHA = " << omp_get_thread_num() << ";" << id << ";" << (int) colormap[id].R << ";" << (int) colormap[id].G << ";" << (int) colormap[id].B << "--" << (int) colormap[id].ALPHA << endl;}cout << endl << endl;

		CreateColormap( colormap );

		cout << "Colormapping of orientations was successful" << endl;

//cout << endl << endl; for ( unsigned int id = 0; id < SubgrainsInTotal; id++ ) { cout << "thread/sid/R/G/B/colormap[id].ALPHA = " << omp_get_thread_num() << ";" << id << ";" << (int) colormap[id].R << ";" << (int) colormap[id].G << ";" << (int) colormap[id].B << "--" << (int) colormap[id].ALPHA << endl;}

		//implicit barrier after pragma omp parallel assures colormap is filled before filling the IPF mapping with cached color values
		gtime = omp_get_wtime();

		unsigned char* ipf = new unsigned char[4*(xmx-xmi+1)*(ymx-ymi+1)];
		unsigned int pxid = 0;
		unsigned int pxc = 0;
		//cout << m_container->getMinX() << ";" << m_container->getMaxX() << "\t\t" << m_container->getMinY() << ";" << m_container->getMaxY() << "\t\t" << m_container->getMinZ() << ";" << m_container->getMaxZ() << endl;
		for ( unsigned int y = ymi; y <= ymx; y++ ) {
			for ( unsigned int x = xmi; x <= xmx; x++ ) {
				//row the y, column the x, and depth the z coordinate of the element
				pxid = (m_container->getValueAt( y, x, z ))-1;
				pxc = 4*((x-xmi) + ((y-ymi) * xlen));
				ipf[pxc + REDCHAN] =  colormap[pxid].R;
				ipf[pxc + GREENCHAN] = colormap[pxid].G;
				ipf[pxc + BLUECHAN] = colormap[pxid].B;
				ipf[pxc + ALPHACHAN] = 255;
			}
		}

		myprofiler.logev("MappingIDs2Colors", (omp_get_wtime() - gtime));

		cout << "Mapping of sub-grain IDs to colors was successful, writing image..." << endl;
		gtime = omp_get_wtime();

		ostringstream fname;
		fname << "Microstructure.SimID." << Settings::SimID << ".IPFZ.2D.png";

		lodepng::encode( fname.str().c_str(), ipf , (xmx - xmi + 1), (ymx - ymi + 1) );	//utilize lodePNG to plot

		delete [] ipf;
		delete [] colormap;

		myprofiler.logev("PlottingStructure2D", (omp_get_wtime() - gtime));
	}
}


void microStructureHdl::Plot3DVolume() {
	if ( Settings::PlotDimension == E_3D ) {
		double gtime = omp_get_wtime();

		double edgelength = Settings::NumberOfGridpoints;
		unsigned int xmi = 0; //Settings::PlotWindowXMin * edgelength;
		unsigned int xmx = m_container->getMaxX(); //Settings::PlotWindowXMax * edgelength;
		unsigned int ymi = 0; //Settings::PlotWindowYMin * edgelength;
		unsigned int ymx =  m_container->getMaxY(); //Settings::PlotWindowYMax * edgelength;
		unsigned int zmi = 0; //Settings::PlotWindowZMin * edgelength;
		unsigned int zmx =  m_container->getMaxZ(); //Settings::PlotWindowZMax * edgelength;
		//if ( xmx >= m_container->getMaxX() ) xmx--; //MK::upper AABB bound now made inclusive!
		//if ( ymx >= m_container->getMaxY() ) ymx--;
		//if ( zmx >= m_container->getMaxZ() ) zmx--;

		//getMaxI are exclusive make them inclusive
		xmx--;
		ymx--;
		zmx--;
		cout << "Generating a grain ID mapping of the domain in 3D on " << xmi << ";" << xmx << ";" << ymi << ";" << ymx << ";" << zmi << ";" << zmx << endl;
		//##MK::see 2D counterpart of the function for how to implement 3D IPF coloring...

		unsigned int* buffer = NULL;
		unsigned int len = xmx - xmi + 1;
		buffer = new unsigned int[len];

		cout << "Buffer length len " << len << endl;
		stringstream fname;
		fname << "Microstructure.SimID." << Settings::SimID << ".GrainID.3D.raw";
		FILE* msvol = fopen(fname.str().c_str(), "wb");
		//filename.str(std::string()); filename.clear();

		//size_t invalid = 0;
		for ( unsigned int z = zmi; z <= zmx; z++ ) {
			for ( unsigned int y = ymi; y <= ymx; y++ ) {
				for( unsigned int i = 0; i < len; ++i )
					buffer[i] = numeric_limits<unsigned int>::max();

				unsigned int valid = 0;
				for ( unsigned int x = xmi; x <= xmx; x++ ) {
					//row the y, column the x, and depth the z coordinate of the element
					buffer[x-xmi] = m_container->getValueAt( y, x, z );
					//if ( m_container->getValueAt( y, x, z ) >= 1 && m_container->getValueAt( y, x, z ) <= 4096 )
					//	valid++;
				}
				//if ( valid != len )
				//	invalid++;

				//##MK::could be improved output line by line, now a small buffer
				fwrite( &buffer[0], sizeof(unsigned int), len, msvol );
			}
		}
		//cout << "Invalid " << invalid << endl;

		delete [] buffer; buffer = NULL;
		cout << "Writing 3D volume representation with " << this->get_first_id() << " as the first ID was successful." << endl;

		fclose(msvol);

		myprofiler.logev( "PlottingStructure3D", (omp_get_wtime() - gtime));
	}
	else {
		cout << "Microstructure case is not three-dimensional!" << endl; return;
	}
}


void microStructureHdl::ReportProfile() {
	cout << endl << "Approximate OpenMP omp_get_wtime report (seconds)" << endl;
	for ( unsigned int e = 0; e < myprofiler.get_entries(); e++ ) {
		cout << myprofiler.titles[e] << ";" << setprecision(6) << myprofiler.times[e] << endl; //separator was "\t\t\t\t"
	}
}



struct NUMANode {
	int num_cpus;
	int numa_cpus[64];
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

void microStructureHdl::initEnvironment() {

//Set up correct Maximum Number of threads
	if (Settings::ExecuteInParallel) {
		Settings::MaximumNumberOfThreads = omp_get_max_threads();

	} else {
		Settings::MaximumNumberOfThreads = 1;
		omp_set_num_threads(Settings::MaximumNumberOfThreads);
	}

//These lines might need to be moved if spatial distribution of grains is utilized
//At best the grain scheduler should be configurable through the parameters file

//m_grainScheduler = new IterativeGrainScheduler(Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);

//choose grain scheduler:
//	if (Settings::GrainScheduler == E_ITERATIVE) {
//		m_grainScheduler = new IterativeGrainScheduler(
//				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
//	} else if (Settings::GrainScheduler == E_SQUARES) {
//		m_grainScheduler = new SquaresGrainScheduler(
//				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
//	} else if (Settings::GrainScheduler == E_DEFAULT_SCHEDULER) {
//		m_grainScheduler = new IterativeGrainScheduler(
//				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
//	}
	m_grainScheduler = new IterativeGrainScheduler(
			Settings::MaximumNumberOfThreads, Settings::NumberOfGrains);
	initNUMABindings();
}

void microStructureHdl::initNUMABindings() {
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
					printf("Will bind thread %d to cpu %d\n",
							omp_get_thread_num(),
							nodes.at(i).numa_cpus[threadID]);
					cpu_set_t set;
					CPU_ZERO(&set);
					CPU_SET(nodes.at(i).numa_cpus[threadID], &set);
					int res = sched_setaffinity(0, sizeof(set), &set);
					printf(res == 0 ? "Successful\n" : "Failed\n");
				}
				break;
			}
			threadID -= nodes.at(i).num_cpus;
		}
	}
}

void microStructureHdl::initializeGrains(vector<vector<Eigen::Vector3d>> hulls,
		vector<double> grainVolume) {
	cout << "Started Initializing Grains " << endl;
	m_grainScheduler->buildThreadWorkloads(hulls, Settings::NumberOfGridpoints);

	bool exceptionHappened = false;
	string error_message;
	m_grains.resize(Settings::NumberOfGrains + 1);

#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
		for (auto id : workload) {
			if (id <= Settings::NumberOfGrains)
				try {
					Grains* newGrain = new Grains(id, hulls[id], m_container,
							grainVolume[id]);
					m_grains[id] = newGrain;

				} catch (exception& e) {
#pragma omp critical
					{
						exceptionHappened = true;
						error_message += string("Grain ")
								+ to_string((unsigned long long) id)
								+ " in its constructor! Reason : " + e.what()
								+ string("\n");
					}
				}

			if (exceptionHappened) {
				throw runtime_error(error_message);
			}
		}
	}
}

//void microStructureHdl::updateGlobalVoxelContainer() {
//#pragma omp parallel
//	{
//		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
//				omp_get_thread_num());
//
//		//man darf hier parallel ausführen, aber ist das performant???
//		for (auto id : workload) {
//			if (id <= Settings::NumberOfGrains) {
//				m_grains[id]->copySubgrainsToGlobalContainer(m_container);
//			}
//		}
//
//	}
//}

void microStructureHdl::saveTexture() {

	string filename = string("Texture") + string(".txt");

	FILE* output = fopen(filename.c_str(), "wt");

	for (auto it : m_grains) {
		for (auto sub : it->get_SubGrains()) {

			myQuaternion* ori = sub->get_Orientation();
			double* bunge = ori->Quaternion2Euler();
			fprintf(output, "%lf \t%lf \t%lf\t %lf \n", bunge[0], bunge[1],
					bunge[2], sub->get_Volume());
			delete[] bunge;
		}
	}

	fclose(output);
}

void microStructureHdl::plotGrains() {
	vector<Grains*>::iterator it;
	for (it = ++m_grains.begin(); it != m_grains.end(); it++) {
		(*it)->plotLocalContainer();
	}
}


void microStructureHdl::testprng( unsigned int n, double mu, double sigma ) {
	//get some output from Marsaglia's SHR3 powered Ziggurat for standardnormal distributed
	/*for ( int tid = 0; tid < 128; tid++ ) {
		m_seqRND->initSHR3( (uint32_t) pow(2.0, 31) - tid - 1 );
		m_seqRND->initR4Uni( (uint32_t) pow(2.0, 31) - tid - 1 );
		m_seqRND->r4_nor_setup();
		for ( unsigned int i = 0; i < n; i++ ) {
			//rescale
			cout << setprecision(6) << (m_seqRND->r4_nor() * sigma) + mu << endl;
		}
	}*/


	for ( int tid = 0; tid < 128; tid++ ) {
		randomClass agen;
		agen.initMT( (uint32_t) pow(2.0, 31) - tid - 1 );
		for ( unsigned int i = 0; i < n; i++ ) {
			cout << setprecision(6) << agen.MersenneTwister() << endl;
		}
	}
}

