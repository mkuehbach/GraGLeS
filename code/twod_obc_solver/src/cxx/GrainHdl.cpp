/*
 GraGLeS 2D A grain growth simulation utilizing level set approaches
 Copyright (C) 2015  Christian Miessen, Nikola Velinov

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
#include "GrainHdl.h"
#include "Settings.h"
#include "GrahamScan.h"
#include "Spoint.h"
#include "RTree.h"
#include "Utilities.h"
//#include "stdafx.h"
#include <sys/time.h>
#include <numa.h>
#include <unistd.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include "rapidxml.hpp"
#include "rapidxml_print.hpp"
#include "IterativeGrainScheduler.h"
#include "SquaresGrainScheduler.h"
#include <fstream>
#include <string>
#include <string.h>
#include "Structs.h"
#ifdef USE_MKL
#include "mkl.h"
#endif

using namespace rapidxml;

grainhdl::grainhdl() :
	m_ThreadPoolCount(0), loop(0), IDField(NULL) {
}
grainhdl::~grainhdl() {
	delete mymath;
	delete TwinBoundary;
	if (Settings::LatticeType == E_HEXAGONAL)
		delete m_misOriHdl;
}

void grainhdl::setSimulationParameter() {
	initEnvironment();
	mymath = new mathMethods();

	if (Settings::LatticeType == E_HEXAGONAL)
		m_misOriHdl = new MisorientationHdl;

	if (Settings::MicrostructureGenMode	== E_READ_VOXELIZED_MICROSTRUCTURE) {
		read_header_from_nexusfile();
	}
	else {
		cerr << "Unsupported microstructureGenMode!" << "\n";
	}

	ds = (Settings::ConstantSectorRadius * 2 + Settings::InterpolatingSectorRadius * 2) / ((double) realDomainSize);
	TwinBoundary = new Quaternion(60 * PI / 180, 1.0, 1.0, 1.0, true);

	// this fancy choose of the timestep ensures that the interface of spherical grain with the fastest mobility
	// and radius 5 moves maximum the distance of 1 gridpoint per timestep
	switch (Settings::ConvolutionMode)
	{
		case E_LAPLACE: {
			dt = 0.8 / ((double) SQR(realDomainSize)) * Settings::NumberOfPointsPerGrain / 2. / 5.;
			TimeSlope = 0.8482;
			break;
		}
		case E_LAPLACE_RITCHARDSON: {
			dt = 0.8 / ((double) SQR(realDomainSize)) * Settings::NumberOfPointsPerGrain / 2. / 5.;
			TimeSlope = 0.8202;
			break;
		}
		case E_GAUSSIAN: {
			//CM::one grid point per integration step defined such that grain with radius 5 migrates at most 1 px per integration step
			dt = 0.8 / ((double) SQR(realDomainSize)) * Settings::NumberOfPointsPerGrain / 2. / 5.; 
			TimeSlope = 0.8359;
			break;
		}
		default:
		{
			break;
		}
	}
	h = 1.0 / double(realDomainSize);

	//Recalculate the setting parameters for the sector radiuses
	Settings::ConstantSectorRadius *= h; //TODO::use public non-static Settings!
	Settings::InterpolatingSectorRadius *= h;

	delta = Settings::DomainBorderSize * 1 / double(realDomainSize);
	grid_blowup = Settings::DomainBorderSize;
	BoundaryGrainTube = grid_blowup;
	ngridpoints = realDomainSize + (2 * grid_blowup);
	boundary = new LSbox(0, 0, 0, 0, this);
	grains.resize(Settings::NumberOfParticles + 1); //TODO::use actual number of grains!

	if ( Settings::MicrostructureGenMode == E_READ_VOXELIZED_MICROSTRUCTURE ) {
		read_microstructure_from_nexusfile();
		cout << "Read microstructure input files successfully" << endl;
	}
	else {
		cerr << "Unable to read microstructure!" << "\n";
	}

	GrainJunction::handler = this;
	find_correctTimestepSize();
	//construct_boundary();

	cout << endl << "******* PROGRAM OPTIONS: *******" << endl << endl;
	cout << "Number of Grains: " << ngrains << endl;
	cout << "simulated Timesteps: " << Settings::NumberOfTimesteps << endl;
	cout << "DELTA TUBE: " << delta << endl;
	cout << "Timestepwidth " << dt << endl;
	cout << "Number of Gridpoints: " << ngridpoints << endl << endl;

	cout << endl << "******* start simulation: *******" << endl << endl;
}

void grainhdl::find_correctTimestepSize() {
	double my_max = 0;
	double my_min = 1000000;
	if (Settings::UseMagneticField) {
		for (int i = 1; i < ngrains; i++) {
			if (grains[i] != NULL) {
				if (grains[i]->get_magneticEnergy() < my_min)
					my_min = grains[i]->get_magneticEnergy();
				if (grains[i]->get_magneticEnergy() > my_max)
					my_max = grains[i]->get_magneticEnergy();
			}
		}

		if (ngrains == 1)
			my_min = 0.0;
		m_Energy_deltaMAX = (my_max - my_min);
	}
	double vmax = (m_Energy_deltaMAX * dt);
	double m_dt_Correction = 0.5/(vmax / h) ;
	;
	if (m_dt_Correction > 1.0)
		m_dt_Correction = 1.0;
	dt *= m_dt_Correction;

}


void grainhdl::read_header_from_nexusfile()
{
	cout << __func__ << "\n";
	HdfFiveSeqHdl h5r = HdfFiveSeqHdl( Settings::AdditionalFilename );
	string grpnm = "/entry1/ms";
	string dsnm = "";

	//key quantities relevant to define the discretization of the simulated VE
	dsnm = grpnm + "/extent";
	vector<unsigned int> nxy = {0, 0};
	if ( h5r.nexus_read(dsnm, nxy) != MYHDF5_SUCCESS ) { return; }
	realDomainSize = nxy[0];

	dsnm = grpnm + "/number_of_subgrains";
	unsigned int n_subgr = 0;
	if ( h5r.nexus_read(dsnm, n_subgr) != MYHDF5_SUCCESS ) { return; }
	Settings::NumberOfParticles = n_subgr;
	ngrains = Settings::NumberOfParticles;
	currentNrGrains = ngrains;
	//Settings::NumberOfPointsPerGrain = (int) (realDomainSize / sqrt(ngrains));
	dsnm = grpnm + "/average_subgrain_discretization";
	vector<unsigned int> dxy;
	if ( h5r.nexus_read( dsnm, dxy) != MYHDF5_SUCCESS ) { return; }
	Settings::NumberOfPointsPerGrain = dxy.at(0);
	cout << "NumberOfPointsPerGrain " << Settings::NumberOfPointsPerGrain << "\n";
}


void grainhdl::read_microstructure_from_nexusfile()
{
	cout << __func__ << "\n";
	HdfFiveSeqHdl h5r = HdfFiveSeqHdl( Settings::AdditionalFilename );
	string grpnm = "/entry1/ms";
	string dsnm = "";

	//key quantities relevant to define the discretization of the simulated VE
	dsnm = grpnm + "/average_subgrain_discretization";
	vector<unsigned int> dxy = {0, 0};
	if ( h5r.nexus_read(dsnm, dxy) != MYHDF5_SUCCESS ) { return; }
	dsnm = grpnm + "/extent";
	vector<unsigned int> nxy = {0, 0};
	if ( h5r.nexus_read(dsnm, nxy) != MYHDF5_SUCCESS ) { return; }
	dsnm = grpnm + "/number_of_subgrains";
	unsigned int n_subgr = 0;
	if ( h5r.nexus_read(dsnm, n_subgr) != MYHDF5_SUCCESS ) { return; }

	//grain-specific quantities
	dsnm = grpnm + "/grain_size";
	vector<double> sz;
	if ( h5r.nexus_read(dsnm, sz) != MYHDF5_SUCCESS ) { return; }
	vector<int> counts = vector<int>(n_subgr + 1, 0);
	for( unsigned int i = 1; i <= n_subgr; i++ ) {
		counts[i] = (int) sz[i-1];
	}
	sz = vector<double>();

	int id_offset = 1;
	dsnm = grpnm + "/grain_identifier";
	vector<int> id_tmp;
	if ( h5r.nexus_read(dsnm, id_tmp) != MYHDF5_SUCCESS ) { return; }
	vector<int> ID = vector<int>(n_subgr + 1, -1);
	//grains with zero count or other problems are marked with -1
	for( unsigned int i = 1; i <= n_subgr; i++ ) {
		if ( counts[i] > 0 ) {
			ID[i] = id_tmp[i-1] - (id_offset - 1);
		}
	}
	id_tmp = vector<int>();

	dsnm = grpnm + "/grain_aabb";
	vector<int> aabb_tmp;
	if ( h5r.nexus_read(dsnm, aabb_tmp) != MYHDF5_SUCCESS ) { return; }
	vector<vector<SPoint>> vertices;
	vertices.resize(n_subgr + 1);
	for( unsigned int i = 1; i <= n_subgr; i++ ) {
		int xmi = aabb_tmp[4*(i-1)+0];
		int xmx = aabb_tmp[4*(i-1)+1];
		int ymi = aabb_tmp[4*(i-1)+2];
		int ymx = aabb_tmp[4*(i-1)+3];
		vertices[i].push_back(SPoint(xmi, ymi, 0, 0));
		vertices[i].push_back(SPoint(xmi, ymx, 0, 0));
		vertices[i].push_back(SPoint(xmx, ymi, 0, 0));
		vertices[i].push_back(SPoint(xmx, ymx, 0, 0));
	}
	aabb_tmp = vector<int>();

	/*
	dsnm = grpnm + "/grain_barycentre_naive";
	vector<double> xy_tmp;
	if ( h5r.nexus_read(dsnm, xy_tmp) != MYHDF5_SUCCESS ) { return; }
	vector<double> xy = vector<double>( 2*(n_subgr + 1), -1);
	for( unsigned int i = 1; i <= n_subgr; i++ ) {
		for ( unsigned int j = 0; j < 2; j++ ) {
			xy[(2*i)+j] = (int) xy_tmp[2*(i-1)+j];
		}
	}
	xy_tmp = vector<double>(); //TODO::int ????
	*/

	dsnm = grpnm + "/grain_orientation";
	vector<double> q4_tmp;
	if ( h5r.nexus_read(dsnm, q4_tmp) != MYHDF5_SUCCESS ) { return; }
	vector<double> quat = vector<double>(4*(n_subgr + 1), 0.); //scalar, vector???
	for ( unsigned int i = 1; i <= n_subgr; i++ ) {
		for( unsigned int j = 0; j < 4; j++ ) {
			quat[(4*i)+j] = q4_tmp[4*(i-1)+j];
		}
	}
	q4_tmp = vector<double>();

	pair<double, double> see_mimx = pair<double, double>(1.e20, 0.);
	vector<double> stored_elastic_energy = vector<double>( n_subgr + 1, 0.);
	if (Settings::UseStoredElasticEnergy) {
		dsnm = grpnm + "/grain_dislocation_density";
		vector<double> see_tmp;
		if ( h5r.nexus_read(dsnm, see_tmp) != MYHDF5_SUCCESS ) { return; }
		for( unsigned int i = 1; i <= n_subgr; i++ ) {
			stored_elastic_energy[i] = see_tmp[i-1];
			if ( see_tmp[i-1] <= see_mimx.first ) {
				see_mimx.first = see_tmp[i-1];
			}
			if ( see_tmp[i-1] >= see_mimx.second ) {
				see_mimx.second = see_tmp[i-1];
			}
		}
		see_tmp = vector<double>();
	}

	double real_time_scaling = 1. / Settings::DislocEnPerM * Settings::HAGB_Energy / Settings::Physical_Domain_Size;
	//1/((N/m^2)*m^2)*(J/m^2)/m = (1/N)*Nm/m = 1
	// find maximum Velocity driven exclusively by Stored Elastic Energy
	// adapt timestep with
	if (Settings::UseStoredElasticEnergy) {
		see_mimx.first = see_mimx.first * real_time_scaling;
		see_mimx.second = see_mimx.second * real_time_scaling;
		if (ngrains == 1) {
			see_mimx.first = 0.;
		}
		m_Energy_deltaMAX = Settings::DislocEnPerM * (see_mimx.second - see_mimx.first) / Settings::HAGB_Energy * Settings::Physical_Domain_Size;
	}

	dsnm = grpnm + "/subgrain_structure";
	vector<unsigned int> tmp;
	if( h5r.nexus_read( dsnm, tmp ) != MYHDF5_SUCCESS ) { return; }

	IDField = new DimensionalBuffer<int> (0, 0, ngridpoints, ngridpoints);
	cout << "ngridpoints " << ngridpoints << "\n";
	cout << "grid_blowup " << grid_blowup << "\n";
	for (int j = 0; j < ngridpoints; j++) {
		for (int i = 0; i < ngridpoints; i++) {
			IDField->setValueAt(j, i, 0);
		}
	}

	for (int j = grid_blowup; j < (ngridpoints - grid_blowup); j++) {
		int yoff = (j - grid_blowup) * nxy[0];
		for (int i = grid_blowup; i < (ngridpoints - grid_blowup); i++) {
			int ixy = (i - grid_blowup) + yoff;
			if ( ixy >= tmp.size() ) {
				cout << "j, i " << j << "\t\t" << i << "\t\t" << ixy << "\n";
			}
			int box_id = (int) tmp.at(ixy);
			box_id = box_id - (id_offset - 1);
			if (box_id < 0) {
				box_id = 0;
			}
			IDField->setValueAt(j, i, box_id);
		}
	}
	tmp = vector<unsigned int>();

	cout << "Build boxes..." << "\n";

	buildBoxVectors(ID, vertices, quat, stored_elastic_energy);
}


void grainhdl::distanceInitialisation() {
	for (int ii = 0; ii < omp_get_max_threads(); ii++) {
		vector<unsigned int>& workload =
				m_grainScheduler->getThreadWorkload(ii);
		cout << "workload thread " << ii << "= " << workload.size()
				<< " number of grains to process " << endl;
		//		for	(auto id : workload) {
		//			cout << id  << "  " ;
		//		}
		cout << endl;
	}
#pragma omp parallel

	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles)
		if (grains[id] != NULL)
		grains[id]->calculateDistanceFunction();
	}
}
if (IDField != NULL)
delete IDField;

}

void grainhdl::convolution(double& plan_overhead) {
	double timer;
	timeval t;

	plan_overhead = 0;
	gettimeofday(&t, NULL);
	timer = t.tv_sec + t.tv_usec / 1000000;
	createConvolutionPlans();
	gettimeofday(&t, NULL);
	plan_overhead += t.tv_sec + t.tv_usec / 1000000 - timer;
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles)
		if (grains[id] != NULL)
		grains[id]->executeConvolution(
				m_ThreadMemPool[omp_get_thread_num()]);
	}
}
gettimeofday(&t, NULL);
timer = t.tv_sec + t.tv_usec / 1000000;
destroyConvolutionPlans();
gettimeofday(&t, NULL);
plan_overhead += t.tv_sec + t.tv_usec / 1000000 - timer;
}
void grainhdl::createConvolutionPlans() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles)
		if (grains[id] != NULL)
		grains[id]->preallocateMemory(
				m_ThreadMemPool[omp_get_thread_num()]);
	}
}
#ifdef USE_FFTW
for(unsigned int i=0; i<Settings::MaximumNumberOfThreads; i++)
{
	vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(i);
	for(auto & id : workload)
	if(grains[id] != NULL)
	grains[id]->makeFFTPlans(m_ThreadMemPool[i]);
}
#endif
}

void grainhdl::destroyConvolutionPlans() {
#ifdef USE_FFTW
	for(unsigned int i=1; i<grains.size(); i++)
	{
		if(grains[i] != NULL)
		grains[i]->destroyFFTWs();
	}
#endif
}

void grainhdl::comparison_box() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles)
		if (grains[id] != NULL) {
			grains[id]->executeComparison();
			grains[id]->executeSetComparison();
		}
	}
}
}

void grainhdl::level_set() {
#pragma omp parallel

	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles) {
			if (grains[id] == NULL)
			continue;
			if (grains[id]->grainExists() == false) {
				delete grains[id];
				grains[id] = NULL;
			} else
			grains[id]->extractContour();
		}
	}
}
}

void grainhdl::redistancing() {

#pragma omp parallel
	{

		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles) {
			if (grains[id] == NULL)
			continue;
			grains[id]->computeVolumeAndEnergy();
			grains[id]->executeRedistancing();
			//#pragma omp atomic
			//				currentNrGrains += 1;
		}

	}
}
}


vector<unsigned int> grainhdl::get_nexus_grain_identifier()
{
	vector<unsigned int> retval;
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists()) {
			retval.push_back( grains[i]->getID() );
		}
	}
	return retval;
}


vector<double> grainhdl::get_nexus_grain_size()
{
	vector<double> retval;
	double A = Settings::Physical_Domain_Size * Settings::Physical_Domain_Size;
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists()) {
			retval.push_back( grains[i]->getVolume() ); //*0.5 ??
		}
	}
	return retval;
}


vector<double> grainhdl::get_nexus_grain_stored_elastic_energy()
{
	vector<double> retval;
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists()) {
			retval.push_back( grains[i]->get_StoredElasticEnergy() );
		}
	}
	return retval;
}


vector<unsigned char> grainhdl::get_nexus_grain_edge_contact()
{
	vector<unsigned char> retval;
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists()) {
			retval.push_back( grains[i]->has_edge_contact() );
		}
	}
	return retval;
}


vector<double> grainhdl::get_nexus_grain_orientation()
{
	vector<double> retval;
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists() ) {
			vector<double> quat = grains[i]->get_quaternion();
			retval.insert( retval.end(), quat.begin(), quat.end());
			/*
			for( size_t j = 0; j < bunge.size(); j++ ) {
				retval.push_back( bunge[j] );
			}
			*/
		}
	}
	return retval;
}


vector<double> grainhdl::get_nexus_grain_barycentre()
{
	vector<double> retval;
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists() ) {
			vector<double> xy = grains[i]->get_barycentre();
			retval.insert( retval.end(), xy.begin(), xy.end() );
			/*
			for( size_t j = 0; j < xy.size(); j++ ) {
				f64.push_back( xy[j] );
			}
			*/
		}
	}
	return retval;
}


void grainhdl::get_nexus_grain_boundary_vertices( vector<double> & vrts )
{
	nx_vrts_offsets = vector<size_t>( grains.size(), 0 );
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists() ) {
			vector<double> i_vrts;
			nx_vrts_offsets[i] = grains[i]->get_contour_vertices( i_vrts );
			vrts.insert( vrts.end(), i_vrts.begin(), i_vrts.end() );
			i_vrts = vector<double>();
		}
	}
}


void grainhdl::get_nexus_grain_boundary_xdmf_topology( vector<unsigned int> & inds )
{
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists() ) {
			vector<unsigned int> i_inds;
			unsigned int i_offset = 0;
			for ( int j = 0; j < i; j++ ) {
				i_offset += nx_vrts_offsets[j];
			}
			grains[i]->get_contour_xdmf_topology( i_offset, i_inds );
			inds.insert( inds.end(), i_inds.begin(), i_inds.end() );
			i_inds = vector<unsigned int>();
		}
	}
}


void grainhdl::get_nexus_grain_boundary_xdmf_grain_indices( vector<unsigned int> & grain_ids )
{
	//individual indices on pointer array grains
	//which is not necessarily for each grain the return value of grains[i]->getID() !
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists() ) {
			grain_ids.push_back( i );
		}
	}
}


void grainhdl::get_nexus_grain_boundary_info( vector<double> & ifo )
{
	for (int i = 1; i < grains.size(); i++) {
		if (grains[i] != NULL && grains[i]->grainExists() ) {
			vector<double> i_ifo;
			grains[i]->get_contour_xdmf_info( i_ifo );
			ifo.insert( ifo.end(), i_ifo.begin(), i_ifo.end() );
			i_ifo = vector<double>();
		}
	}
}


bool grainhdl::save_NeXus()
{
	HdfFiveSeqHdl h5w = HdfFiveSeqHdl( Settings::ResultsFileName );
	ioAttributes anno = ioAttributes();
	string grpnm = "";
	string dsnm = "";
	vector<unsigned int> u32;
	vector<double> f64;
	vector<unsigned char> u8;

	grpnm = "/entry1/ms/step" + to_string(loop);
	anno = ioAttributes();
	anno.add( "NX_class", string("NXms_snapshot") );
	if ( h5w.nexus_path_exists( grpnm ) == false ) {
		if ( h5w.nexus_write_group( grpnm, anno ) != MYHDF5_SUCCESS ) { return false; }
	}
	else {
		return true;
	}

	dsnm = grpnm + "/real_time";
	double rt = Realtime;
	anno = ioAttributes();
	anno.add( "unit", string("s ?????") );
	if ( h5w.nexus_write( dsnm, rt, anno ) != MYHDF5_SUCCESS ) { return false; }

	dsnm = grpnm + "/real_domain_size";
	int rd = realDomainSize;
	anno = ioAttributes();
	if ( h5w.nexus_write( dsnm, rd, anno ) != MYHDF5_SUCCESS ) { return false; }

	dsnm = grpnm + "/current_number_of_grains";
	unsigned int ngr = currentNrGrains;
	anno = ioAttributes();
	if ( h5w.nexus_write( dsnm, ngr, anno ) != MYHDF5_SUCCESS ) { return false; }



	dsnm = grpnm + "/grain_identifier";
	u32 = get_nexus_grain_identifier();
	anno = ioAttributes();
	if ( h5w.nexus_write(
		dsnm,
		io_info({u32.size()}, {u32.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		u32,
		anno ) != MYHDF5_SUCCESS ) { return false; }
	u32 = vector<unsigned int>();

	dsnm = grpnm + "/grain_size";
	f64 = get_nexus_grain_size();
	anno = ioAttributes();
	//anno.add( "unit", string("m") );
	if ( h5w.nexus_write(
		dsnm,
		io_info({f64.size()}, {f64.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		f64,
		anno ) != MYHDF5_SUCCESS ) { return false; }
	f64 = vector<double>();

	dsnm = grpnm + "/grain_stored_elastic_energy";
	f64 = get_nexus_grain_stored_elastic_energy();
	anno = ioAttributes();
	//anno.add( "unit", string("1/m^2") );
	if ( h5w.nexus_write(
		dsnm,
		io_info({f64.size()}, {f64.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		f64,
		anno ) != MYHDF5_SUCCESS ) { return false; }
	f64 = vector<double>();

	if ( loop == 1 ) {
		dsnm = grpnm + "/grain_orientation";
		f64 = get_nexus_grain_orientation();
		anno = ioAttributes();
		//anno.add( "unit", string("°") );
		if ( h5w.nexus_write(
			dsnm,
			io_info({f64.size() / 4, 4}, {f64.size() / 4, 4}, MYHDF5_COMPRESSION_GZIP, 0x01),
			f64,
			anno ) != MYHDF5_SUCCESS ) { return false; }
		f64 = vector<double>();
	}

	/*
	dsnm = grpnm + "/grain_barycentre";
	f64 = get_nexus_grain_barycentre();
	anno = ioAttributes();
	//anno.add( "unit", string("m") );
	if ( h5w.nexus_write(
		dsnm,
		io_info({f64.size() / 2, 2}, {f64.size() / 2, 2}, MYHDF5_COMPRESSION_GZIP, 0x01),
		f64,
		anno ) != MYHDF5_SUCCESS ) { return false; }
	f64 = vector<double>();
	*/

	dsnm = grpnm + "/grain_edge_contact";
	u8 = get_nexus_grain_edge_contact();
	anno = ioAttributes();
	if ( h5w.nexus_write(
		dsnm,
		io_info({u8.size()}, {u8.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
		u8,
		anno ) != MYHDF5_SUCCESS ) { return false; }
	u8 = vector<unsigned char>();

	if ( (loop - Settings::StartTime) % Settings::NetworkExport == 0) {
		vector<double> f64_vrts;
		//should be implemented with incrementally writing i.e.
		//first a probe_nexus_grain_boundary to find the bounds
		dsnm = grpnm + "/grain_boundary_vertices";
		get_nexus_grain_boundary_vertices( f64_vrts );
		anno = ioAttributes();
		if ( h5w.nexus_write(
			dsnm,
			io_info({f64_vrts.size() / 2, 2}, {f64_vrts.size() / 2, 2}, MYHDF5_COMPRESSION_GZIP, 0x01),
			f64_vrts,
			anno ) != MYHDF5_SUCCESS ) { return false; }
		f64_vrts = vector<double>();

		vector<unsigned int> u32_xdmf_topo;
		dsnm = grpnm + "/grain_boundary_xdmf_topology";
		get_nexus_grain_boundary_xdmf_topology( u32_xdmf_topo );
		anno = ioAttributes();
		if ( h5w.nexus_write(
			dsnm,
			io_info({u32_xdmf_topo.size()}, {u32_xdmf_topo.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
			u32_xdmf_topo,
			anno ) != MYHDF5_SUCCESS ) { return false; }
		u32_xdmf_topo = vector<unsigned int>();

		dsnm = grpnm + "/grain_boundary_xdmf_grain_id";
		get_nexus_grain_boundary_xdmf_grain_indices( u32 );
		anno = ioAttributes();
		if ( h5w.nexus_write(
			dsnm,
			io_info({u32.size()}, {u32.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
			u32,
			anno ) != MYHDF5_SUCCESS ) { return false; }
		u32 = vector<unsigned int>();

		vector<double> f64_ifo;
		dsnm = grpnm + "/grain_boundary_energy_times_mobility";
		get_nexus_grain_boundary_info( f64_ifo );
		anno = ioAttributes();
		//anno.add( "unit", string("(J/m^2)*(m^4/J/s)") );
		if ( h5w.nexus_write(
			dsnm,
			io_info({f64_ifo.size()}, {f64_ifo.size()}, MYHDF5_COMPRESSION_GZIP, 0x01),
			f64_ifo,
			anno ) != MYHDF5_SUCCESS ) { return false; }
		f64_ifo = vector<double>();
	}

	return true;
}


double parallelRest = 0;
double convo_time = 0;
double comparison_time = 0;
double levelset_time = 0;
double redistancing_time = 0;
double plan_overhead = 0;
void grainhdl::run_sim() {
	timeval time;
	double timer;
	double overhead;
	gettimeofday(&time, NULL);
	timer = time.tv_sec + time.tv_usec / 1000000.0;
	distanceInitialisation();
	get_biggestGrainVol();
	gettimeofday(&time, NULL);
	cout << "Time for Distancefunction Initialization: " << time.tv_sec
			+ time.tv_usec / 1000000.0 - timer << endl;

	Realtime = 0.0;
	find_neighbors();
	for (loop = Settings::StartTime; loop <= Settings::StartTime
			+ Settings::NumberOfTimesteps; loop++) {
		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		gridCoarsement();
		gettimeofday(&time, NULL);
		parallelRest += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		convolution(overhead);
		plan_overhead += overhead;
		gettimeofday(&time, NULL);
		convo_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		switchDistancebuffer();
		updateSecondOrderNeighbors();
		gettimeofday(&time, NULL);
		parallelRest += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		if (Settings::DecoupleGrains != 1) {
			gettimeofday(&time, NULL);
			timer = time.tv_sec + time.tv_usec / 1000000.0;
			comparison_box();
			gettimeofday(&time, NULL);
			comparison_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

			gettimeofday(&time, NULL);
			timer = time.tv_sec + time.tv_usec / 1000000.0;
			switchDistancebuffer();
			gettimeofday(&time, NULL);
			parallelRest += time.tv_sec + time.tv_usec / 1000000.0 - timer;
		}

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		level_set();
		gettimeofday(&time, NULL);
		levelset_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;

		gettimeofday(&time, NULL);
		timer = time.tv_sec + time.tv_usec / 1000000.0;
		redistancing();
		gettimeofday(&time, NULL);
		redistancing_time += time.tv_sec + time.tv_usec / 1000000.0 - timer;
		//		switchDistancebuffer();
		if (((loop - Settings::StartTime) % int(Settings::AnalysisTimestep))
				== 0 || loop == Settings::NumberOfTimesteps) {

			if ( save_NeXus() == true ) {
				cout << "Writing snapshot data for loop " << loop << " into NeXus file success." << "\n";
			}
			else {
				cerr << "Writing snapshot data for loop " << loop << " into NeXus file failed!" << "\n";
			}
			/*
			save_TextureFaces_Binary(); //MK::save_Texture(); MODF writing disabled to reduce number of files
			*/
			//if (Settings::MicrostructureGenMode == E_READ_VOXELIZED_MICROSTRUCTURE && loop == 0)
			//	save_Full_Microstructure_for_Restart();
			//! With activated ResearchMode the id to centroid assignment is plotted
		}
		Realtime += (dt * (Settings::Physical_Domain_Size
				* Settings::Physical_Domain_Size) / (TimeSlope
				* Settings::HAGB_Energy * Settings::HAGB_Mobility)); //##MK correction ok?

cout << "I am incrementing the real-time realTime/dt/PhysDomSize/realDomainSize/TimeSlope/Energy/Mobility/GridCoarsement = "
	<< Realtime << ";" << dt << ";" << Settings::Physical_Domain_Size << ";" << realDomainSize << ";" << TimeSlope << ";" << Settings::HAGB_Energy << ";" << Settings::HAGB_Mobility << "--" << (int) Settings::GridCoarsement << endl;

		get_biggestGrainVol();

		if (currentNrGrains < Settings::BreakupNumber) {
			//check which type of grain boundaries, in the old code with an totalenergy.back() >0.98*PI/2)
			cout << "Network has coarsed to less than 3% of the population. "
					<< "Remaining Grains: " << currentNrGrains
					<< " Break and save." << endl;
			/*
			save_GBContourPlot();
			save_GBJunctionPlot();
			save_TextureFaces_Binary(); //MK::save_Texture(); MODF writing disabled to reduce number of files
			*/
			break;
		}

	}
	// 	utils::CreateMakeGif();
	cout << "Simulation complete." << endl;
	cout << "Simulation Time: " << Realtime << endl;
	cout << "Detailed timings: " << endl;
	cout << "Convolution time: " << convo_time << endl;
	cout << "     Of which plan overhead is: " << plan_overhead << endl;
	cout << "Comparison time: " << comparison_time << endl;
	cout << "Redistancing time: " << redistancing_time << endl;
	cout << "Levelset time: " << levelset_time << endl;
	cout << "GridCoarse/SwitchBuffer/UpNeigh: " << parallelRest << endl;
	cout << "Sum parallel regions: " << convo_time + comparison_time
			+ levelset_time + parallelRest + redistancing_time << endl;

#ifdef USE_MKL
	int numMKLBuffer;
	mkl_mem_stat(&numMKLBuffer);
	cout << "Memory used by MKL : "<< numMKLBuffer << endl;
#endif

}


void grainhdl::updateSecondOrderNeighbors()
{
	#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(omp_get_thread_num());
		for	(auto id : workload) {
			if (grains[id] != NULL) {
				grains[id]->computeSecondOrderNeighbours();
			}
		}
	}
}


void grainhdl::find_neighbors()
{
	RTree<unsigned int, int, 2, float> tree;
	int min[2], max[2];
	for (unsigned int i = 1; i <= Settings::NumberOfParticles; i++) {
		if (grains[i] != NULL) {
			min[0] = grains[i]->getMinX();
			min[1] = grains[i]->getMinY();
			max[0] = grains[i]->getMaxX();
			max[1] = grains[i]->getMaxY();
			tree.Insert(min, max, i);
		}
	}

	#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(omp_get_thread_num());
		for	(auto id : workload) {
			if (grains[id] != NULL) {
				grains[id]->computeDirectNeighbours(tree);
			}
		}
	}
}


void grainhdl::switchDistancebuffer()
{
	#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(omp_get_thread_num());
		for	(auto id : workload) {
			if (grains[id] != NULL) {
				grains[id]->switchInNOut();
			}
		}
	}
}


void grainhdl::gridCoarsement()
{
	int newSize = pow((double) currentNrGrains, (1./3.)) * pow(PI * 4./3., 1./3.) / 2. * Settings::NumberOfPointsPerGrain;
	double population_reduction = (double) currentNrGrains / (double) ngrains;
	if (Settings::GridCoarsement == true && loop != 0 && population_reduction < Settings::GridCoarsementGradient ) {
		//int newSize = sqrt(currentNrGrains) * Settings::NumberOfPointsPerGrain;
		cout << "Coarsing the current grid in step " << loop << " to newSize " << newSize << "\n";
		#pragma omp parallel
		{
			vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(omp_get_thread_num());
			for (auto id : workload) {
				//if (id <= Settings::NumberOfParticles) {
				if (grains[id] != NULL) {
					grains[id]->resizeGrid(newSize);
				}
				//}
			}
		}

		realDomainSize = newSize;
		delta = Settings::DomainBorderSize * 1 / double(realDomainSize);
		ngridpoints = realDomainSize + 2 * grid_blowup;
		h = 1.0 / realDomainSize;
		//! DISCREPANCY: Compare to the application of dt in the convolution, time decreasing factor 0.8

		switch (Settings::ConvolutionMode) {
			case E_LAPLACE: {
				dt = 0.8 / ((double) SQR(realDomainSize)) * Settings::NumberOfPointsPerGrain / 2. / 5.;
				break;
			}
			case E_LAPLACE_RITCHARDSON: {
				dt = 0.8 / ((double) SQR(realDomainSize)) * Settings::NumberOfPointsPerGrain / 2. / 5.;
				break;
			}
			//CM::one grid point per integration step defined such that grain with radius 5 migrates at most 1 px per integration step
			case E_GAUSSIAN: {
				dt = 0.8 / ((double) SQR(realDomainSize)) * Settings::NumberOfPointsPerGrain / 2. / 5.;
				break;
			}
			default: {
				break;
			}
		}
		double m_dt_Correction = 0.5 / realDomainSize / m_Energy_deltaMAX / dt;
		if (m_dt_Correction > 1.0) {
			cerr << "Why is m_dt_Correction > 1.0 ?" << "\n";
		}
		dt *= m_dt_Correction;
		ngrains = currentNrGrains;

		#pragma omp parallel
		{
			vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(omp_get_thread_num());
			for (auto id : workload) {
				if (id <= Settings::NumberOfParticles) {
					if (grains[id] == NULL)
					continue;
					grains[id]->recalculateIDLocal();
				}
			}
		}
		#pragma omp parallel
		{
			vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(omp_get_thread_num());
			for (auto id : workload) {
				if (id <= Settings::NumberOfParticles) {
					if (grains[id] == NULL)
					continue;
					grains[id]->extractContour();
				}
			}
		}
	}
	else {
		switchDistancebuffer();
	}
}


void grainhdl::initEnvironment() {
	//Set up correct Maximum Number of threads
	if (Settings::ExecuteInParallel) {
		Settings::MaximumNumberOfThreads = omp_get_max_threads();

	} else {
		Settings::MaximumNumberOfThreads = 1;
		omp_set_num_threads(Settings::MaximumNumberOfThreads);
	}

	m_ThreadPoolCount = Settings::MaximumNumberOfThreads;
	m_ThreadMemPool.resize(m_ThreadPoolCount);

	//These lines might need to be moved if spatial distribution of grains is utilized
	//At best the grain scheduler should be configurable through the parameters file

	//m_grainScheduler = new IterativeGrainScheduler(Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);

	//choose grain scheduler:
	if (Settings::GrainScheduler == E_ITERATIVE) {
		m_grainScheduler = new IterativeGrainScheduler(
				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
	} else if (Settings::GrainScheduler == E_SQUARES) {
		m_grainScheduler = new SquaresGrainScheduler(
				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
	} else if (Settings::GrainScheduler == E_DEFAULT_SCHEDULER) {
		m_grainScheduler = new IterativeGrainScheduler(
				Settings::MaximumNumberOfThreads, Settings::NumberOfParticles);
	}
	initNUMABindings();
#pragma omp parallel
	{
		double max_size = Settings::NumberOfPointsPerGrain
				* Settings::NumberOfPointsPerGrain * 50;
		int power_of_two = 1 << (int) (ceil(log2(max_size)) + 0.5);
		//!int power_of_two = 1 << (int) (ceil(log2(2<<20)) + 0.5); //!27
		m_ThreadMemPool[omp_get_thread_num()].resize(power_of_two);
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

void grainhdl::initNUMABindings() {
	vector<NUMANode> nodes;
	nodes.reserve(16);
	numa_available();
	bitmask* mask = numa_get_run_node_mask();
	bitmask* cpus = numa_allocate_cpumask();
	for (unsigned int j = 0; j < mask->size; j++) {
		if (numa_bitmask_isbitset(mask, j)) {
			cout << "We are allowed to used node " << j << "\n";
			NUMANode node;
			memset(&node, 0xFF, sizeof(node));
			numa_node_to_cpus(j, cpus);
			node.num_cpus = my_numa_bitmask_weight(cpus);
			int cpuCounter = 0;
			for (unsigned int i = 0; i < cpus->size; i++) {

				if (numa_bitmask_isbitset(cpus, i) && numa_bitmask_isbitset(
						numa_all_cpus_ptr, i)) {
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
					cout << "Will bind thread " << omp_get_thread_num() << " to CPU " << nodes.at(i).numa_cpus[threadID];
					cpu_set_t set;
					CPU_ZERO(&set);
					CPU_SET(nodes.at(i).numa_cpus[threadID], &set);
					int res = sched_setaffinity(0, sizeof(set), &set);
					if ( res == 0 ) {
						cout << ", success";
					}
					else {
						cout << ", failed!";
					}
					cout << "\n";
				}
				break;
			}
			threadID -= nodes.at(i).num_cpus;
		}
	}
}

void grainhdl::buildBoxVectors(vector<vector<SPoint>>& contours) {
	m_grainScheduler->buildGrainWorkloads(contours, ngridpoints);
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles) {
			LSbox* grain = new LSbox(id, contours[id], this);
			grains[id] = grain;
		}
	}
}
}

void grainhdl::buildBoxVectors(vector<vector<SPoint>>& contours,
		vector<double>& q1, vector<double>& q2, vector<double>& q3,
		vector<double>& q4) {
	m_grainScheduler->buildGrainWorkloads(contours, ngridpoints);
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles) {
			//catch the grains with default q1 = -1000
			if (q1[id] == -1000) {
				grains[id] = NULL;
				continue;
			}
			LSbox* grain = new LSbox(id, contours[id], q1[id], q2[id],
					q3[id], q4[id], this);
			grains[id] = grain;
		}
	}
}
}

void grainhdl::buildBoxVectors(int* ID, vector<vector<SPoint>>& contours,
		Quaternion* Quaternionen, double* StoredElasticEnergy) {
	m_grainScheduler->buildGrainWorkloads(contours, ngridpoints);
	cout << "Construct LSbox objects" << endl;
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles) {
			if (ID[id] == -1) {
				grains[id] = NULL;
				continue;
			}
			LSbox* grain;
			if (Settings::UseStoredElasticEnergy) {
				grain = new LSbox(ID[id], contours[id], Quaternionen[id],
						StoredElasticEnergy[id], this);
			} else {
				grain = new LSbox(ID[id], contours[id], Quaternionen[id], 0,
						this);
			}
			grains[id] = grain;
		}
	}
}
}


void grainhdl::buildBoxVectors(vector<int> & ID, vector<vector<SPoint>> & contours,
	vector<double> & q, vector<double> & see )
{
	m_grainScheduler->buildGrainWorkloads(contours, ngridpoints);
	cout << "Construct LSbox objects" << endl;
	#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(omp_get_thread_num());
		for	(auto id : workload) {
			if (id <= Settings::NumberOfParticles) {
				if (ID[id] == -1) {
					grains[id] = NULL;
					continue;
				}
				LSbox* grain;
				double stored_elastic_energy = (Settings::UseStoredElasticEnergy) ? see[id] : 0.;
				vector<double> qt;
				for( unsigned int i = 0; i < 4; i++ ) { qt.push_back(q[(4*id)+i]); }

				grain = new LSbox(ID[id], contours[id], qt, stored_elastic_energy, this);
				grains[id] = grain;
			}
		}
	}
}


void grainhdl::set_h(double hn)
{
	h = hn;
}


void grainhdl::set_realDomainSize(int realDomainSizen)
{
	realDomainSize = realDomainSizen;
	ngridpoints = realDomainSize + 2 * grid_blowup;
}


void grainhdl::setResearchAdjustments()
{
	//convolutionCorrection = true;
	//! Calculate regression lines?
	//calcRegression = false;
	//! Calculate centroids?
	//calcCentroid = false;
	//! Is there a need to load curvature in the first timesteps?
	//loadCurvature = false;
	//! As the reference the observed values of the four-sided-grain case
	//! is used.
	//loadCurvatureLoop = ((Settings::NumberOfPointsPerGrain / 40.0)
	//		* (Settings::NumberOfPointsPerGrain / 40.0) * 200);
	return;
}


void grainhdl::get_biggestGrainVol() {
	currentNrGrains = 0;
	maxVol = 0.0;
	for (auto it = (++grains.begin()); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		if (maxVol < (*it)->getVolume())
			maxVol = (*it)->getVolume();
		currentNrGrains++;
	}
}
