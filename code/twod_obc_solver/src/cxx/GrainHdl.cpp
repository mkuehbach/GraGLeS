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

	if (Settings::MicrostructureGenMode == E_GENERATE_TESTCASE) {
		Settings::NumberOfParticles = read_ScenarioPoints();
	}
	if (Settings::MicrostructureGenMode == E_READ_VERTEX) {
		FILE * getNumber;
		getNumber = fopen(Settings::ReadFromFilename.c_str(), "r");
		fscanf(getNumber, "%d\n", &Settings::NumberOfParticles);
		fclose(getNumber);
	} else if (Settings::MicrostructureGenMode
			== E_READ_VOXELIZED_MICROSTRUCTURE) {
		//read_HeaderCPG();
		read_header_from_nexusfile();
	} else {
		ngrains = Settings::NumberOfParticles;
		currentNrGrains = ngrains;
		realDomainSize = sqrt(Settings::NumberOfParticles)
				* Settings::NumberOfPointsPerGrain; // half open container of VORO++
	}

	ds = (Settings::ConstantSectorRadius * 2
			+ Settings::InterpolatingSectorRadius * 2) / double(realDomainSize);
	MODF.resize(Settings::DiscreteSamplingRate);
	TwinBoundary = new Quaternion(60 * PI / 180, 1.0, 1.0, 1.0, true);
	fill(MODF.begin(), MODF.end(), 0);

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

	switch (Settings::MicrostructureGenMode) {
	case E_GENERATE_WITH_VORONOY: {
		if (Settings::UseTexture) {
			bunge = new double[3] { PI / 2, PI / 2, PI / 2 };
			deviation = Settings::MaxMisOrientation * PI / 180;
		} else {
			bunge = NULL;
			deviation = 0;
		}
		ST = NULL;

		VOROMicrostructure();
		// 			generateRandomEnergy();
		break;
	}
	case E_READ_VERTEX: {
		bunge = NULL;
		deviation = 0;
		ST = new double[ngrains * ngrains];
		std::fill_n(ST, ngrains * ngrains, 0);
		readMicrostructureFromVertex();

		break;
	}
	case E_READ_FROM_FILE: {
		cout << "Reading file..." << endl;
		bunge = NULL;
		deviation = 0;

		ST = NULL;
		readMicrostructure();
		break;
	}
	case E_READ_VOXELIZED_MICROSTRUCTURE: {
		cout << "Starting to read microstructure input files" << endl;
		//read_voxelized_microstructure();
		read_microstructure_from_nexusfile();
		cout << "Read microstructure input files successfully" << endl;
		break;
	}
	case E_INVALID_VAL: {
		cout << "Invalid Read mode. Correct the parameters.xml and try again!"
				<< endl;
		exit(2);
	}
		//!
		//! This case handles the processing of
		//! grain construction by means of a file input
		//! with 2D point information
		//!
	case E_GENERATE_TESTCASE: {
		if (Settings::UseTexture) {
			bunge = new double[3] { PI / 2, PI / 2, PI / 2 };
			deviation = 15 * PI / 180;
		} else {
			bunge = NULL;
			deviation = 0;
		}
		ST = NULL;
		VOROMicrostructure();
		// 			generateRandomEnergy();

		//! Read weights from a file
		if (Settings::MicrostructureGenMode == E_GENERATE_TESTCASE
				&& !Settings::IsIsotropicNetwork) {

			//! Read from file in vector

			ifstream mapStream("surfaceTension.dat");

			if (!mapStream) {
				cerr << "File can't be opened." << endl;
				exit(-1);
			}

			while (mapStream) {
				string line;
				getline(mapStream, line);
				istringstream iss(line);
				vector<double> values;
				copy(istream_iterator<double> (iss),
						istream_iterator<double> (),
						back_insert_iterator<vector<double> > (values));
				if (values.size() > 0)
					weightsMatrix.push_back(values);
			}

			cout << "Line: " << weightsMatrix.size() << endl;
			cout << "Column (in the first line): " << weightsMatrix[0].size()
					<< endl;
		}

		break;
	}
	}
	GrainJunction::handler = this;
	find_correctTimestepSize();
	// 	construct_boundary();
	//program options:
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
void grainhdl::read_HeaderCPG() {
	FILE * compressedGrainInfo;
	compressedGrainInfo = fopen(Settings::AdditionalFilename.c_str(), "r");
	if (compressedGrainInfo == NULL) {
		cout << "Could not read from specified file !";
		exit(2);
	}

	// Read header

	char c;
	string name;
	char buffer[100];
	int I_D;
	int DX;
	int DY;
	int NX;
	int NY;

	for (int i = 1; i < 8; i++) {
		if (i == 1) {
			do {
				c = fgetc(compressedGrainInfo);
				//cout << c;
			} while (c != '\n');

			continue;
		}
		if (i == 2) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &I_D);
			//cout << buffer << "\t" << I_D << endl;
		}
		if (i == 3) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &DX);
			//cout << buffer << "\t" << DX << endl;
		}
		if (i == 4) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &DY);
			//cout << buffer << "\t" << DY << endl;
		}
		if (i == 5) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &realDomainSize);
			//cout << buffer << "\t" << NX << endl;
		}
		if (i == 6) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &NY);
			//cout << buffer << "\t" << NY << endl;
		}
		if (i == 7) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer,
					&Settings::NumberOfParticles);
			//cout << buffer << "\t" << Settings::NumberOfParticles << endl;
		}
	}
	ngrains = Settings::NumberOfParticles;
	currentNrGrains = ngrains;
	//Settings::NumberOfPointsPerGrain = (int) (realDomainSize / sqrt(ngrains));
	// half open container of VORO++
	fclose(compressedGrainInfo);
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


void grainhdl::read_voxelized_microstructure() {
	FILE * compressedGrainInfo;
	compressedGrainInfo = fopen(Settings::AdditionalFilename.c_str(), "rt");
	if (compressedGrainInfo == NULL) {
		cout << "Could not read from specified file !";
		exit(2);
	}
	rewind(compressedGrainInfo);
	// Read header

	char c;
	string name;
	char buffer[100];
	int ID_offset;
	int DX;
	int DY;
	int NX;
	int NY;

	for (int i = 1; i <= 8; i++) {
		if (i == 1 || i == 8) {
			do {
				c = fgetc(compressedGrainInfo);
			} while (c != '\n');
			continue;
		}
		if (i == 2) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &ID_offset);
			//cout << buffer << "\t" << I_D << endl;
		}
		if (i == 3) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &DX);
			//cout << buffer << "\t" << DX << endl;
		}
		if (i == 4) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &DY);
			//cout << buffer << "\t" << DY << endl;
		}
		if (i == 5) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &NX);
			//cout << buffer << "\t" << NX << endl;
		}
		if (i == 6) {
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &NY);
			//cout << buffer << "\t" << NY << endl;
		}
		if (i == 7) {
			int grains = 0;
			fscanf(compressedGrainInfo, "%s \t %i\n", buffer, &grains);
			//cout << buffer << "\t" << grains << endl;
		}
	}

	double bunge[3];
	double *x, *y;
	double* StoredElasticEnergy;
	int* ID;
	int xmin, xmax, ymin, ymax;
	int *counts;
	vector < vector < SPoint >> vertices;
	vertices.resize(ngrains + 1);
	Quaternion* Quaternionen = new Quaternion[ngrains + 1];
	ID = new int[ngrains + 1];
	if (Settings::UseStoredElasticEnergy) {
		StoredElasticEnergy = new double[ngrains + 1];
	}
	double bufferStored;
	counts = new int[ngrains + 1];
	x = new double[ngrains + 1];
	y = new double[ngrains + 1];
	int id = 0;
	double StoredElasticEnergy_MIN = 1e20;
	double StoredElasticEnergy_MAX = 0;
	for (int nn = 1; nn <= ngrains; nn++) {
		//ID, x, y, bunge1, bunge2, bunge3, xmin, xmax, ymin, ymax, vol, stored energy
		fscanf(
				compressedGrainInfo,
				"%d\t %lf\t %lf\t %lf\t %lf\t %lf\t %d\t %d\t %d\t %d\t %d\t %lf \n",
				&id, &x[nn], &y[nn], &bunge[0], &bunge[1], &bunge[2], &xmin,
				&xmax, &ymin, &ymax, &counts[nn], &bufferStored);
		if (Settings::UseStoredElasticEnergy) {
			StoredElasticEnergy[nn] = bufferStored;
			if (abs(StoredElasticEnergy[nn]) < StoredElasticEnergy_MIN)
				StoredElasticEnergy_MIN = abs(StoredElasticEnergy[nn]);
			if (abs(StoredElasticEnergy[nn]) > StoredElasticEnergy_MAX)
				StoredElasticEnergy_MAX = abs(StoredElasticEnergy[nn]);
		}
		//printf("%d\t %lf\t %lf\t %lf\t %lf\t %lf\t %d\t %d\t %d\t %d\t %d \n",
		//		id, x[nn], y[nn], bunge[0], bunge[1], bunge[2], xmin,
		//		xmax, ymin, ymax, counts[nn]);
		if (counts[nn] == 0) {
			//skips zero pixel grains, with invalid box boundaries
			ID[nn] = -1;
			continue;
		}
		ID[nn] = id - (ID_offset - 1);
		Quaternionen[nn].euler2quaternion(bunge);
		vertices[nn].push_back(SPoint(xmin, ymin, 0, 0));
		vertices[nn].push_back(SPoint(xmin, ymax, 0, 0));
		vertices[nn].push_back(SPoint(xmax, ymin, 0, 0));
		vertices[nn].push_back(SPoint(xmax, ymax, 0, 0));
	}

	fclose(compressedGrainInfo);
	// find maximum Velocity driven exclusively by Stored Elastic Energy
	// adapt timestep with
	if (Settings::UseStoredElasticEnergy) {
		if (abs(boundary->get_StoredElasticEnergy()) / Settings::DislocEnPerM
				* Settings::HAGB_Energy / Settings::Physical_Domain_Size
				< StoredElasticEnergy_MIN)
			StoredElasticEnergy_MIN = abs(boundary->get_StoredElasticEnergy())
					/ Settings::DislocEnPerM * Settings::HAGB_Energy
					/ Settings::Physical_Domain_Size;
		if (abs(boundary->get_StoredElasticEnergy()) / Settings::DislocEnPerM
				* Settings::HAGB_Energy / Settings::Physical_Domain_Size
				> StoredElasticEnergy_MAX)
			StoredElasticEnergy_MAX = abs(boundary->get_StoredElasticEnergy())
					/ Settings::DislocEnPerM * Settings::HAGB_Energy
					/ Settings::Physical_Domain_Size;
		if (ngrains == 1)
			StoredElasticEnergy_MIN = 0.0;
		m_Energy_deltaMAX = Settings::DislocEnPerM * (StoredElasticEnergy_MAX
				- StoredElasticEnergy_MIN) / Settings::HAGB_Energy
				* Settings::Physical_Domain_Size;
	}
	//BINARY read-in

	FILE * voxelized_data;

	voxelized_data = fopen(Settings::ReadFromFilename.c_str(), "rb");
	if (voxelized_data == NULL) {
		cout << "Could not read from specified file !";
		exit(2);
	}

	IDField = new DimensionalBuffer<int> (0, 0, ngridpoints, ngridpoints);
	//	int min_ID = 100;

	for (int j = 0; j < ngridpoints; j++) {
		for (int i = 0; i < ngridpoints; i++) {
			if (i < grid_blowup || j < grid_blowup || i >= ngridpoints
					- grid_blowup || j >= ngridpoints - grid_blowup)
				IDField->setValueAt(j, i, 0);
			else {
				int box_id;
				fread(&box_id, sizeof(int), 1, voxelized_data);
				box_id = box_id - (ID_offset - 1);
				if (box_id < 0)
					box_id = 0;
				IDField->setValueAt(j, i, box_id);
				//				if (box_id < min_ID)
				//					min_ID = box_id;
			}

		}
	}
	//	fstream datei("Voxel.dat", ios::out);
	//	for (int i = 0; i < ngridpoints; i++) {
	//		for (int j = 0; j < ngridpoints; j++) {
	//			datei << IDField->getValueAt(i, j) << " | ";
	//		}
	//		datei << endl;
	//	}
	//	fclose(voxelized_data);

	buildBoxVectors(ID, vertices, Quaternionen, StoredElasticEnergy);

	delete[] ID;
	delete[] Quaternionen;
	delete[] x;
	delete[] y;
	delete[] counts;
	if (Settings::UseStoredElasticEnergy) {
		delete[] StoredElasticEnergy;
	}
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


void grainhdl::VOROMicrostructure() {

	stringstream filename, plotfiles;
	int cell_id;
	double x, y, z, rx, ry, rz;

	grains.resize(Settings::NumberOfParticles + 1);

	vector<LSbox*> local_grains;

	std::vector<LSbox*>::iterator itg;

	// stores the centroids of the cells ; access by (3*Id, 3*ID +1, 3*ID +2)
	part_pos = new double[3 * Settings::NumberOfParticles];

	bool randbedingung = false; // bei false ist der container halb offen?! d.h. gitterwert mit 1 werden keinem partikel zugeordnet
	if (randbedingung == false)
		realDomainSize -= 1;

	voronoicell_neighbor c;
	container con(0, 1, 0, 1, 0, 1, 5, 5, 5, randbedingung, randbedingung,
			randbedingung, 2);
	c_loop_all vl(con);

	//!
	//! Particles are added deliberately in the container according to the input file data.
	//!
	if (Settings::MicrostructureGenMode == E_GENERATE_TESTCASE) {
		FILE* pointSketch;
		pointSketch = fopen(Settings::ReadFromFilename.c_str(), "r");
		if (pointSketch == nullptr) {

			cout << "There does not exist a file";
			exit(1);
		}

		double pointX, pointY;

		for (int k = 0; k < ngrains; k++) {

			fscanf(pointSketch, "%lf\t%lf\n", &pointX, &pointY);
			con.put(k, pointX, pointY, 0);
		}

		fclose(pointSketch);
	}

	/**********************************************************/
	// Randomly add particles into the container
	if (Settings::MicrostructureGenMode != E_GENERATE_TESTCASE) {
		for (int i = 0; i < ngrains; i++) {
			x = rnd();
			y = rnd();
			z = 0;
			con.put(i, x, y, z);
		}
	}
	/**********************************************************/

	for (int i = 0; i < realDomainSize; i++)
		for (int j = 0; j < realDomainSize; j++) {
			x = double(i * h);
			y = double(j * h); // only point within the domain
			if (con.find_voronoi_cell(x, y, z, rx, ry, rz, cell_id)) {
				cell_id++;
				part_pos[3 * (cell_id - 1)] = rx;
				part_pos[3 * (cell_id - 1) + 1] = ry;
				part_pos[3 * (cell_id - 1) + 2] = rz;
			}

			else
				fprintf(stderr, "# find_voronoi_cell error for %g %g 0\n", x, y);

		}

	vector < vector < SPoint >> m_initialContours;
	m_initialContours.resize(Settings::NumberOfParticles + 1);

	if (vl.start())
		do {
			con.compute_cell(c, vl);
			int box_id = vl.pid() + 1;
			GrahamScan scanner(c, box_id, part_pos);
			scanner.generateCovnexHull(m_initialContours[box_id]);

		} while (vl.inc());

	buildBoxVectors( m_initialContours);

	delete[] part_pos;
}

void grainhdl::readMicrostructure() {
	FILE * levelset;
	levelset = fopen(Settings::ReadFromFilename.c_str(), "r");
	if (levelset == NULL) {
		cout << "Could not read from specified file !";
		exit(2);
	}
	int id;
	double xl, yl;

	int nvertices;
	//double* vertices = new double[1000];
	vector < vector<SPoint> > vertices;
	vertices.resize(ngrains + 1);
	vector<double> q1, q2, q3, q4;
	q1.resize(ngrains + 1);
	q2.resize(ngrains + 1);
	q3.resize(ngrains + 1);
	q4.resize(ngrains + 1);

	vertices.resize(ngrains + 1);

	for (int nn = 1; nn <= ngrains; nn++) {
		fscanf(levelset, "%d\t", &id);
		q1[nn] = -1000;
		if (id >= ngrains) {
			grains.resize(id + 1);
			q1.resize(id + 1);
			q2.resize(id + 1);
			q3.resize(id + 1);
			q4.resize(id + 1);
			Settings::NumberOfParticles = id;
			ngrains = id;
			currentNrGrains = id;
		}

		fscanf(levelset, "%d\t %lf\t %lf\t%lf\t%lf\n", &nvertices, &q1[nn],
				&q2[nn], &q3[nn], &q4[nn]);
		for (int j = 0; j < nvertices; j++) {
			fscanf(levelset, "%lf\t %lf\n", &xl, &yl);
			if (xl < 0)
				xl = 0;
			if (yl < 0)
				yl = 0;
			if (xl > 1)
				xl = 1;
			if (yl > 1)
				yl = 1;
			vertices[nn].push_back(SPoint(xl, yl, 0, 0));
		}
		fscanf(levelset, "\n");

	}
	fclose(levelset);

	buildBoxVectors(vertices, q1, q2, q3, q4);
}

void grainhdl::readMicrostructureFromVertex() {
	FILE * levelset;
	levelset = fopen(Settings::ReadFromFilename.c_str(), "r");

	long id;
	int nedges;
	double phi1, PHI, phi2, xr, yr, xl, yl;
	double* edges;

	fscanf(levelset, "%d\n", &ngrains);
	cout << "ngrains : " << ngrains << endl;

	grains.resize(ngrains + 1);

	for (int nn = 0; nn < ngrains; nn++) {

		fscanf(levelset, "%ld\t %d\t %lf\t %lf\t%lf\n", &id, &nedges, &phi1,
				&PHI, &phi2);
		edges = new double[nedges * 4];
		//		cout << id << " || " << nedges << " || " << phi1 << " || " << PHI
		//				<< " || " << phi2 << endl;

		for (int j = 0; j < nedges; j++) {
			fscanf(levelset, "%lf\t %lf\t %lf\t%lf\n", &xl, &yl, &xr, &yr);
			//			cout << xl << " ||\t " << yl << " ||\t " << xr << " ||\t " << yr
			//					<< " ||\t " << endl;
			int k = 4 * j;
			edges[k] = xl;
			edges[k + 1] = yl;
			edges[k + 2] = xr;
			edges[k + 3] = yr;
		}

		LSbox* newBox = new LSbox(id, nedges, edges, phi1, PHI, phi2, this);
		grains[id] = newBox;

		delete[] edges;
	}

	ST = new double[ngrains * ngrains]; //Create ST array and fill with zeros
	std::fill_n(ST, ngrains * ngrains, 0);

	for (int i = 0; i < ngrains; i++) {
		double buffer;
		fscanf(levelset, "%lf\t", &buffer);
		for (int j = 0; j < ngrains; j++) {
			while (j < i) {
				fscanf(levelset, "%lf\t", &buffer);
				j++;
			}
			fscanf(levelset, "%lf\t", &buffer);
			ST[j + (ngrains * i)] = (double) buffer;
			ST[i + (ngrains * j)] = ST[j + (ngrains * i)];
		}
		fscanf(levelset, "\n");
	}
	fclose(levelset);

	for (int i = 0; i < ngrains; i++) {
		for (int j = 0; j < ngrains; j++) {
			cout << ST[i + (ngrains * j)] << "  \t";
		}
		cout << endl;
	}
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


void grainhdl::save_Texture() {
	cerr << "grainhdl::save_Texture FIX ME !!" << "\n";
	/*
	ofstream myfile;
	FILE* enLenDis;
	stringstream filename;
	double totalLength = 0;
	double total_energy = 0.0;

	filename.str("");
	filename << "MODF_" << loop << ".txt";
	enLenDis = fopen(filename.str().c_str(), "w");

	double dh;
	if (Settings::LatticeType == 0)
		dh = 62.8 / (double) Settings::DiscreteSamplingRate;
	else
		dh = 92. / (double) Settings::DiscreteSamplingRate;
	double *euler = new double[3];
	vector<characteristics>::iterator it2;
	vector<LSbox*>::iterator it;

	std::fill(MODF.begin(), MODF.end(), 0.0);

	if (Settings::NeighbourTracking && loop > 0) {
		double timer = omp_get_wtime();
		stringstream filename;
		filename << "Texture" << "_" << loop << ".bin";
		FILE* binaryTexture = fopen(filename.str().c_str(), "wb");
		filename.str(std::string());
		filename.clear();
		filename << "Faces" << "_" << loop << ".bin";
		FILE* binaryFaces = fopen(filename.str().c_str(), "wb");

		for (int i = 1; i < grains.size(); i++) {
			if (grains[i] != NULL && grains[i]->grainExists()) {
				totalLength += grains[i]->getVolume() * 0.5;
				total_energy += grains[i]->getEnergy() * 0.5;
				fwrite( &grains[i]->collectTextureData(), sizeof(TextureData), 1, binaryTexture ); //##MK::better struct pointer
				vector<Face>* myfaces = grains[i]->get_Faces();
				for (auto it : *myfaces) {
					fwrite(&(it), sizeof(Face), 1, binaryFaces);
				}
				delete myfaces;
			}
		}
		fclose(binaryFaces);
		fclose(binaryTexture);
		double spent = omp_get_wtime();
		cout << "Texture/Faces binary I/O " << filename.str().c_str() << " took " << spent << " second" << endl;
	} else {
		filename.str("");
		filename << "Texture" << "_" << loop << ".ori";
		myfile.open(filename.str());
		for (it = ++grains.begin(); it != grains.end(); it++) {
			if (*it != NULL && (*it)->grainExists() == true
					&& (*it)->isMotionRegular() == true) {
				total_energy += (*it)->getEnergy();
				euler = (*it)->getOrientationQuat()->quaternion2EulerConst();

				totalLength += (*it)->getPerimeter() * 0.5;
				map<int, double>& localMODF = (*it)->getlocalMODF();
				int MODFsize = MODF.size();
				for (const auto& iterator : localMODF) {
					if (iterator.first < MODFsize) {
						MODF[iterator.first] += iterator.second / 2;
					}
				}
				//! If ResearchMode is activated additional data (e.g. the mean (euclidean) triple-junction-distance) is stored in the Texture files
				if (Settings::ResearchMode) {
					double mean_a = (*it)->getMeanA();
					//! Determine the mean number "m" of the direct neighbors' number of faces

					double mean_m = (*it)->getMeanM();

					myfile << (*it)->getID() << '\t'
							<< (*it)->getDirectNeighbourCount() << '\t'
							<< (*it)->intersectsBoundaryGrain() << '\t'
							<< (*it)->getVolume() << '\t' << (*it)->getMeanDa()
							/ Settings::AnalysisTimestep << '\t'
							<< (*it)->getPerimeter() << '\t'
							<< (*it)->getEnergy() << '\t';
					if (Settings::UseMagneticField == 1)
						myfile << (*it)->get_magneticEnergy() << '\t';
					else
						myfile << (*it)->get_StoredElasticEnergy() << '\t';
					myfile << mean_a << '\t' << mean_m << '\t' << euler[0]
							<< '\t' << euler[1] << '\t' << euler[2];

				} else {
					myfile << (*it)->getID() << '\t'
							<< (*it)->getDirectNeighbourCount() << '\t'
							<< (*it)->intersectsBoundaryGrain() << '\t'
							<< (*it)->getVolume() << '\t' << (*it)->getMeanDa()
							/ Settings::AnalysisTimestep << '\t'
							<< (*it)->getPerimeter() << '\t'
							<< (*it)->getEnergy() << '\t';
					if (Settings::UseMagneticField == 1)
						myfile << (*it)->get_magneticEnergy() << '\t';
					else
						myfile << (*it)->get_StoredElasticEnergy() << '\t';
					myfile << euler[0] << '\t' << euler[1] << '\t' << euler[2];
				}
				//				if (Settings::NeighbourTracking) {
				//					vector<int> ids = (*it)->getDirectNeighbourIDs();
				//					vector<double> lengths = (*it)->getGBLengths();
				//					if (ids.size() != lengths.size())
				//						cout << "length differences in ID and Length array! "
				//								<< endl;
				//					myfile << '\t'
				//							<< (*it)->getMinX()
				//									+ (((*it)->getMaxX() - (*it)->getMinX()) / 2)
				//											* h;
				//					myfile << '\t'
				//							<< (*it)->getMinY()
				//									+ (((*it)->getMaxY() - (*it)->getMinY()) / 2)
				//											* h;
				//					for (int i = 0; i < (*it)->getDirectNeighbourCount(); i++) {
				//						myfile << '\t' << ids[i] << '\t' << lengths[i];
				//					}
				//				}

				myfile << endl;

				double sum = 0.0;
				for (unsigned int i = 0; i < Settings::DiscreteSamplingRate; i++) {
					sum += MODF[i];
				}
				for (unsigned int i = 0; i < Settings::DiscreteSamplingRate; i++) {
					if (MODF[i] > 0)
						MODF[i] = MODF[i] / sum;
				}

			}
		}
		delete[] euler;
		myfile.close();
	}
	if (!Settings::IsIsotropicNetwork) {

		for (unsigned int i = 0; i < Settings::DiscreteSamplingRate; i++) {
			fprintf(enLenDis, "%lf\t%lf\n", (float) (dh * (i + 1)),
					(float) MODF[i]);
			//printf("%lf\t%lf\n", (float) (dh * (i + 1)), (float) MODF[i]);
		}

	}
	totalenergy.push_back(0.5 * total_energy);
	nr_grains.push_back(currentNrGrains);
	time.push_back(Realtime);
	cout << "Timestep " << loop << " complete:" << endl;
	cout << "Number of grains remaining in the Network :" << currentNrGrains
			<< endl;
	cout << "Amount of free Energy in the Network :" << 0.5 * total_energy
			<< endl;
	cout << "Total GB Length in Network :" << totalLength << endl << endl
			<< endl;
	myfile.close();
	fclose(enLenDis);
	*/
}


void grainhdl::save_TextureFaces_Binary() {
	cerr << "grainhdl::save_TextureFaces_Binary FIXME !!" << "\n";
	/*
	double timer = omp_get_wtime();
	vector<LSbox*>::iterator it;

	if (Settings::NeighbourTracking && loop > 0) {
		stringstream filename;
		filename << "Texture" << "_" << loop << ".bin";
		FILE* binaryTexture = fopen(filename.str().c_str(), "wb");
		filename.str(std::string());
		filename.clear();
		filename << "Faces" << "_" << loop << ".bin";
		FILE* binaryFaces = fopen(filename.str().c_str(), "wb");

		for (int i = 1; i < grains.size(); i++) {
			if ( grains[i] != NULL && grains[i]->grainExists() ) {
				fwrite(&grains[i]->collectTextureData(), sizeof(TextureData), 1, binaryTexture);
				vector<Face>* myfaces = grains[i]->get_Faces();
				for (auto it : *myfaces) {
					fwrite(&(it), sizeof(Face), 1, binaryFaces);
				}
				delete myfaces;
			}
		}
		fclose(binaryFaces);
		fclose(binaryTexture);
	}
	cout << "GBTextureAndFaces binary I/O Texture_" << loop << ".bin and Faces_" << loop << ".bin in " << setprecision(6) << (double) (omp_get_wtime() - timer) << endl;
	*/
}


void grainhdl::save_RealtimeLog()
{
	double timer = omp_get_wtime();
	double total_energy = 0.0;
	vector<LSbox*>::iterator it;

	if (Settings::NeighbourTracking && loop > 0) {
		for (int i = 1; i < grains.size(); i++) {
			if ( grains[i] != NULL && grains[i]->grainExists() ) { //most likely
				total_energy += grains[i]->getEnergy() * 0.5;
			}
		}
	}

	totalenergy.push_back(0.5 * total_energy);
	nr_grains.push_back(currentNrGrains);
	time.push_back(Realtime);

	cout << "Timestep " << loop << " complete;" << (int) currentNrGrains << " grains remaining" << endl; //, I/O time; " << setprecision(6) << (double) (omp_get_wtime() - timer) << endl;
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
		} else {
			tweakIDLocal();
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
			if ( ((loop - Settings::StartTime) % ( int(Settings::NetworkExport * Settings::AnalysisTimestep) ) == 0) && loop != 0 ) {
				//##MK::for debug purposes still do plotting of network
				//save_NetworkPlot();
				/*
				save_GBCurvApprx();
				save_GBContourPlot();
				save_GBJunctionPlot();
				*/
			}

			//if ( loop >= 4000 )

			if ( save_NeXus() == true ) {
				cout << "Writing snapshot data for loop " << loop << " into NeXus file success." << "\n";
			}
			else {
				cerr << "Writing snapshot data for loop " << loop << " into NeXus file failed!" << "\n";
			}
			/*
			save_TextureFaces_Binary(); //MK::save_Texture(); MODF writing disabled to reduce number of files
			*/
			save_RealtimeLog();
			save_NrGrainsStats();
			//if (Settings::MicrostructureGenMode == E_READ_VOXELIZED_MICROSTRUCTURE && loop == 0)
			//	save_Full_Microstructure_for_Restart();
			//! With activated ResearchMode the id to centroid assignment is plotted
			if (Settings::ResearchMode && calcCentroid)
				save_id();
			if (loop == Settings::NumberOfTimesteps)
				save_Full_Microstructure_for_Restart();
			//save_Memory_Print();
		}
		Realtime += (dt * (Settings::Physical_Domain_Size
				* Settings::Physical_Domain_Size) / (TimeSlope
				* Settings::HAGB_Energy * Settings::HAGB_Mobility)); //##MK correction ok?

cout << "I am incrementing the real-time realTime/dt/PhysDomSize/realDomainSize/TimeSlope/Energy/Mobility/GridCoarsement = "
	<< Realtime << ";" << dt << ";" << Settings::Physical_Domain_Size << ";" << realDomainSize << ";" << TimeSlope << ";" << Settings::HAGB_Energy << ";" << Settings::HAGB_Mobility << "--" << (int) Settings::GridCoarsement << endl;

		get_biggestGrainVol();
		if(Settings::DecoupleGrains && ngrains == 1 && totalenergy.back() >0.98*PI/2)
			break;
		if (currentNrGrains < Settings::BreakupNumber) {
			cout << "Network has coarsed to less than 3% of the population. "
					<< "Remaining Grains: " << currentNrGrains
					<< " Break and save." << endl;
			//save_NetworkPlot();
			/*
			save_GBContourPlot();
			save_GBJunctionPlot();
			save_TextureFaces_Binary(); //MK::save_Texture(); MODF writing disabled to reduce number of files
			save_GBCurvApprx();
			*/
			save_RealtimeLog();

			save_NrGrainsStats();
			save_Full_Microstructure_for_Restart();
			//! With activated ResearchMode the id to centroid assignment is plotted
			//                        if(Settings::ResearchMode && calcCentroid)
			//                                                        save_id();
			//                                                                                if (loop == Settings::NumTimesteps)
			//                                                                                                             saveMicrostructure();
			//
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

void grainhdl::tweakIDLocal() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles) {
			if (grains[id] == NULL)
			continue;
			grains[id]->setIDLocal(boundary->getID());
		}
	}
}
}

void grainhdl::save_Full_Microstructure_for_Restart() {
	stringstream param_xml_name;
	param_xml_name << "PARAMETERS_SIM_NONAME_TIMESTEP_" << loop << "_GRAINS_"
			<< currentNrGrains << ".xml";
	stringstream vertex_dump_name;
	vertex_dump_name << "NETWORK_NONAME_TIMESTEP_" << loop << "_GRAINS_"
			<< currentNrGrains << ".dat";

	createParamsForSim(param_xml_name.str().c_str(),
			vertex_dump_name.str().c_str());

	ofstream output;
	output.open(vertex_dump_name.str());
	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL || (*it)->grainExists() == false)
			continue;
		if (loop != 0 || Settings::MicrostructureGenMode
				== E_READ_VOXELIZED_MICROSTRUCTURE)
			(*it)->plot_full_grain(loop, false, &output, true);
		else
			(*it)->plot_full_grain(loop, false, &output, false);
	}
	output.close();

}
void grainhdl::createParamsForSim(const char* param_filename,
		const char* vertex_dump_filename) {
	xml_document<> doc_tree;

	xml_node<>* declaration = doc_tree.allocate_node(node_declaration);
	declaration->append_attribute(doc_tree.allocate_attribute("version", "1.0"));
	declaration->append_attribute(
			doc_tree.allocate_attribute("encoding", "utf-8"));
	doc_tree.append_node(declaration);

	doc_tree.append_node(
			Settings::generateXMLParametersNode(&doc_tree,
					vertex_dump_filename, loop, currentNrGrains));
	ofstream output;
	output.open(param_filename);
	output << doc_tree;
	output.close();

}

void grainhdl::save_NrGrainsStats() {
	// 	(*my_weights).plot_weightmap(ngridpoints, ID, ST, zeroBox);
	ofstream myfile;
	if (loop == 0)
		myfile.open("NrGrains&EnergyStatistics.txt");
	else
		myfile.open("NrGrains&EnergyStatistics.txt", ios::app);

	if ( time.size() > 0 && nr_grains.size() > 0 && totalenergy.back() > 0 ) {
		myfile << time.back() << "\t";
		myfile << nr_grains.back() << "\t";
		myfile << totalenergy.back() << "\t";
		myfile << realDomainSize << endl;
	}

	myfile.close();

	// 	if (SAVEIMAGE)utils::PNGtoGIF("test.mp4");
	//cout << "number of distanzmatrices: "<< domains.size() << endl;

}

void grainhdl::updateSecondOrderNeighbors() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles) {
			if (grains[id] == NULL)
			continue;
			grains[id]->computeSecondOrderNeighbours();
		}
	}
}
}

void grainhdl::find_neighbors() {
	RTree<unsigned int, int, 2, float> tree;
	int min[2], max[2];
	for (unsigned int i = 1; i <= Settings::NumberOfParticles; i++) {
		if (grains[i] == NULL)
			continue;
		min[0] = grains[i]->getMinX();
		min[1] = grains[i]->getMinY();
		max[0] = grains[i]->getMaxX();
		max[1] = grains[i]->getMaxY();
		tree.Insert(min, max, i);
	}

#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles) {
			if (grains[id] == NULL)
			continue;
			grains[id]->computeDirectNeighbours(tree);
		}
	}
}
}

void grainhdl::saveSpecialContourEnergies(int id) {
	if (grains[id] == NULL)
		return;
	grains[id]->plot_box_contour(loop, true);
}

void grainhdl::save_NetworkPlot() {
	double timer = omp_get_wtime();
	ofstream output;
	stringstream filename;
	filename << "Network_Timestep_" << loop << ".gnu";
	output.open(filename.str());

	std::vector<LSbox*>::iterator it;
	//for (it = ++grains.begin(); it != grains.end(); it++) {
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++) {
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(i);
		for	(auto id : workload) {
			if (grains[id] == NULL)
				continue;
			grains[id]->plot_box_contour(loop, true, &output, true, i);
		}
	}
	output.close();

	double spent = omp_get_wtime() - timer;
	cout << "NetworkPlot ascii I/O " << filename.str().c_str() << " in " << spent << " second" << endl;
}

void grainhdl::save_GBCurvApprx() {
	double tic = omp_get_wtime();
	//store estimation of curvature distribution along contour after redistancing
	//this is a lower bound estimate of the capillary driving forces acting by fitting circle to set of three points

	//first parallel processing of results for disjoint grains
	#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload( omp_get_thread_num() );
		for	(auto id : workload) {
			if (id <= Settings::NumberOfParticles)
				if (grains[id] != NULL)
					grains[id]->approximateCurvature( this->get_grid_blowup() );
		}
	} //implicit barrier

	double toc = omp_get_wtime();
	cout << "GBCurvatureApproximation computation in parallel in " << (toc-tic) << endl;


	//secondly storing results to file sequentially
	tic = omp_get_wtime();
	ofstream output;
	stringstream filename;
	filename << "GBCurvatureApproximation_" << loop << ".bin";
	FILE* binaryCurv = fopen(filename.str().c_str(), "wb");

	std::vector<LSbox*>::iterator it;
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++) {
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(i);
		for	(auto id : workload) {
			if ( grains[id] == NULL )
				continue;
			if ( grains[id]->grainExists() == false )
				continue;

			grains[id]->writeGBCurvatureApprx( binaryCurv );
		}
	}

	fclose(binaryCurv);
	toc = omp_get_wtime();
	cout << "GBCurvatureApproximation binary I/O " << filename.str().c_str() << " in " << (toc-tic) << endl;
}

void grainhdl::save_GBContourPlot() {
	double timer = omp_get_wtime();
	//store geometry of the boundary after redistancing, this is enables to compute a lower bound estimate of the
	//capillary driving force by fitting circle to set of three points
	ofstream output;
	stringstream filename;
	filename << "GBContourPoints_" << loop << ".bin";
	FILE* binaryContours = fopen(filename.str().c_str(), "wb");

	std::vector<LSbox*>::iterator it;
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++) {
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(i);
		for	(auto id : workload) {
			if ( grains[id] == NULL )
				continue;
			if ( grains[id]->grainExists() == false )
				continue;

			grains[id]->writeGBContourPoints( binaryContours, this->get_grid_blowup() );
		}
	}

	fclose(binaryContours);
	double spent = omp_get_wtime() - timer;
	cout << "GBContour binary I/O " << filename.str().c_str() << " in " << spent << endl;
}


void grainhdl::save_GBJunctionPlot() {
	double timer = omp_get_wtime();
	//store geometry of the junctions after redistancing
	ofstream output;
	stringstream filename;
	filename << "GBJunctionPoints_" << loop << ".bin";
	FILE* binaryJunctions = fopen(filename.str().c_str(), "wb");

	std::vector<LSbox*>::iterator it;
	for (unsigned int i = 0; i < Settings::MaximumNumberOfThreads; i++) {
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(i);
		for	(auto id : workload) {
			if ( grains[id] == NULL )
				continue;
			if ( grains[id]->grainExists() == false )
				continue;

			grains[id]->writeGBJunctionPoints( binaryJunctions, this->get_grid_blowup() );
		}
	}

	fclose(binaryJunctions);
	double spent = omp_get_wtime() - timer;
	cout << "GBJunction binary I/O " << filename.str().c_str() << " in " << spent << endl;
}



void grainhdl::save_id() {
	ofstream output;
	stringstream filename;
	filename << "Id_At_Network_Timestep_" << loop << ".gnu";
	output.open(filename.str());

	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		SPoint centroid = (*it)->getCentroid();
		output << centroid.x << "\t" << centroid.y << "\t" << (*it)->getID()
				<< endl;
	}
	output.close();

	//! Plot junction points
	stringstream filename2;
	filename2 << "Junctions_At_Network_Timestep_" << loop << ".gnu";
	output.open(filename2.str());

	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		else
			(*it)->plot_grain_junctions(loop, &output);
	}
	output.close();

}

void grainhdl::save_regLine() {
	ofstream output;
	stringstream filename;
	filename << "RegLine_At_Network_Timestep_" << loop << ".gnu";
	output.open(filename.str());

	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		for (unsigned int j = 0; j < (*it)->getRegressionPoints().size(); j++) {

			output << (*it)->getRegressionPoints()[j].x << "\t"
					<< (*it)->getRegressionPoints()[j].y << "\t" << endl;
			if ((j - 2) % 3 == 0)
				output << "\n" << endl;
		}

		//!output << "\n" << endl;
	}
	output.close();
}

void grainhdl::saveAllContourLines() {
	stringstream filename;
	filename << "NetworkAtTime_" << loop << ".gnu";
	ofstream dateiname;
	dateiname.open(filename.str());
	std::vector<LSbox*>::iterator it;
	for (it = ++grains.begin(); it != grains.end(); it++) {
		if (*it == NULL)
			continue;
		(*it)->plot_box_contour(loop, false);
	}
	dateiname.close();
}
void grainhdl::removeGrain(int id) {
	grains[id] = NULL;
}

void grainhdl::switchDistancebuffer() {
#pragma omp parallel
	{
		vector<unsigned int>& workload = m_grainScheduler->getThreadWorkload(
				omp_get_thread_num());
for	(auto id : workload) {
		if (id <= Settings::NumberOfParticles) {
			if (grains[id] == NULL)
			continue;
			grains[id]->switchInNOut();
		}
	}
}
}

void grainhdl::gridCoarsement() {
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

void grainhdl::clear_mem() {
	if (ST != NULL) {
		delete[] ST;
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


void grainhdl::set_h(double hn) {
	h = hn;
}
void grainhdl::set_realDomainSize(int realDomainSizen) {
	realDomainSize = realDomainSizen;
	ngridpoints = realDomainSize + 2 * grid_blowup;
}
/**
 * This function analyzes the input file for MicrostructureGenMode 4.
 * The amount of lines of the input file is determined. This number
 * indicates the number of points specified in the file, i.e the number of
 * grains. This means one pair of x-y-coordinates each in every line. A
 * tabulator is used as the separator.
 *
 * @return the amount of points in the input file
 */
int grainhdl::read_ScenarioPoints() {

	int counter = 0;
	string line;
	ifstream reader(Settings::ReadFromFilename.c_str());
	while (std::getline(reader, line)) {
		counter++;
	}
	return counter;
}

/**
 * This piece of code cares for the proper program behavior as
 * intended by the current research project.
 */
void grainhdl::setResearchAdjustments(E_RESEARCH_PROJECT pro) {

	//! Some minor control variables to enhance speed if necessary,
	//! It is important to check their properties before starting a
	//! large simulation project.
	//! Activate calibration factor in convolution corrector step
	convolutionCorrection = true;
	//! Calculate regression lines?
	calcRegression = false;
	//! Calculate centroids?
	calcCentroid = false;
	//! Is there a need to load curvature in the first timesteps?
	loadCurvature = false;
	//! As the reference the observed values of the four-sided-grain case
	//! is used.
	loadCurvatureLoop = ((Settings::NumberOfPointsPerGrain / 40.0)
			* (Settings::NumberOfPointsPerGrain / 40.0) * 200);
	//! See the ENUM declaration to get the matching of integer and project
	if (pro == 0)
		project = pro;

	//! This project is characterized by constant energies for every
	//! grain boundary in the network. The few exceptions are the
	//! junctions sections on the domain boundaries which need to be
	//! fixed to their position (weight = 0) because of small deviation
	//! from symmetry in the algorithm.
	//!
	//! Application: Investigation of triple junction and grain boundary
	//! regime by means of Lambda.
	//!
	//! Changes are made in: junction.cpp
	//!
	//! Works with: MicrostructureGenMode == E_READ_VERTEX
	if (pro == 1)
		project = pro;

	//! This project is characterized by constant energies for every
	//! grain boundary in the network.
	//!
	//! Application: Investigation of triple junction and grain boundary
	//! regime by means of Lambda.
	//!
	//! Changes are made in: junction.cpp
	//!
	//! Works with: MicrostructureGenMode == E_GENERATE_WITH_VORONOY,
	//! 			MicrostructureGenMode == E_READ_FROM_FILE
	if (pro == 2)
		project = pro;
}

void grainhdl::setResearchAdjustments() {
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
void grainhdl::save_Memory_Print() {
	stringstream name;
	name << "MemoryPrint_T" << loop << ".txt";
	ofstream file;
	file.open(name.str());
	file << ngridpoints << endl;
	for (int i = 1; i < Settings::NumberOfParticles; i++) {
		if (grains[i] != NULL) {
			grains[i]->outputMemoryUsage(file);
		}
	}
	file.close();
}
