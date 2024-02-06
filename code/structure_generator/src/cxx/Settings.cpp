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


//*********************************************************************************************************************************************************************************************
#include "Settings.h"
#include "rapidxml.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>
#include <stdexcept>
#include <string>
#include <math.h>

using namespace std;
using namespace rapidxml;

//*********************************************************************************************************************************************************************************************

//Initialization of the setting variables *****************************************************************************************************************************************************

unsigned int Settings::SimID = 0; //default naming
unsigned int Settings::NumberOfGrains = 0;
unsigned int Settings::NumberOfSubgrains = 0;
unsigned int Settings::NumberOfGridpoints = 0;				//MK::if 0, program detects size automatically, otherwise sets to user-defined size
unsigned long Settings::NumberOfPointsPerSubGrain = 0;
unsigned long Settings::MaximumNumberOfThreads = 1;
double Settings::Settings::SubgrainOriScatter = 0.0;
double Settings::StoredElasticEnergyMax = 0.0;
double Settings::StoredElasticEnergyMin = 0.0;
double Settings::StoredElasticScatterGrain = 0.0;
double Settings::StoredElasticScatterSubgrain = 0.0;
double Settings::DefGrainRelDimensionX = 1.0; //by definition
double Settings::DefGrainRelDimensionY = 1.0;
double Settings::DefGrainRelDimensionZ = 1.0;
double Settings::RandomnessX = 0.0;
double Settings::RandomnessY = 0.0;
double Settings::RandomnessZ = 0.0;

bool Settings::VoronoiPeriodic = false;
bool Settings::ExecuteInParallel = false;
bool Settings::PlotIPF2DSection = false;
bool Settings::StatusHealthy = true;
double Settings::PlotWindowXMin = 0.0;
double Settings::PlotWindowXMax = 1.0;
double Settings::PlotWindowYMin = 0.0;
double Settings::PlotWindowYMax = 1.0;
double Settings::PlotWindowZMin = 0.0;
double Settings::PlotWindowZMax = 1.0;

bool Settings::BreakPerX = false; //by default domains are periodic
bool Settings::BreakPerY = false;
bool Settings::BreakPerZ = false;

bool Settings::IO_DAMASK = false;
bool Settings::IO_HDF5 = false;


//bool Settings::UseOrientationSpace = false;

E_CRYSTAL_STRUCTURE Settings::CrystalStructure = E_DEFAULT_STRUCTURE;
E_MICROGENMODE Settings::MicroGenMode = E_DEFAULT_GENERATOR;
E_SAMPLING Settings::TextureSampling = E_DEFAULT_SAMPLING;
E_TEXTURE Settings::TextureGEN = E_DEFAULT_TEXTURE;
E_PLOT_DIMENSION Settings::PlotDimension = E_DEFAULT_DIMENSION;
E_GRAIN_AGGREGATE Settings::GrainAggregation = E_POLYCRYSTAL;
E_GRAIN_SHAPE Settings::GrainShape = E_GLOBULITIC;

string Settings::ReadFromFilename;
string Settings::AdditionalFilename;

//Definition of the needed functions **********************************************************************************************************************************************************

void Settings::readXML(string filename) {

	//Find the wished .xml file ***********************************************************************************************************************************************************

	if (0 == filename.compare(""))
		filename = string("parameters.xml");
	ifstream file(filename);
	if (file.fail()) {
		throw runtime_error(string("Unable to locate file ") + filename);
	}
	stringstream contents;
	contents << file.rdbuf();
	string xmlDocument(contents.str());

	xml_document<> tree;
	tree.parse<0>(&xmlDocument[0]);

	xml_node<>* rootNode = tree.first_node();
	if (0 != strcmp(rootNode->name(), "Parameters")) {
		throw runtime_error("undefined parameters file!");
	}

	//Read all the required parameters ****************************************************************************************************************************************************
	if (0 != rootNode->first_node("SimulationID")) {
		SimID = std::stoul(
				rootNode->first_node("SimulationID")->value());
	}

	if (0 != rootNode->first_node("PlotDimension")) {
		PlotDimension = (E_PLOT_DIMENSION) std::stoi(
				rootNode->first_node("PlotDimension")->value());

		if (PlotDimension >= E_DEFAULT_DIMENSION)
			PlotDimension = E_DEFAULT_DIMENSION;
	}

	if (0 != rootNode->first_node("NumberOfGrains")) {
		NumberOfGrains = std::stoul(
				rootNode->first_node("NumberOfGrains")->value());
	}

	if (0 != rootNode->first_node("NumberOfSubgrains")) {
		NumberOfSubgrains = std::stoul(
				rootNode->first_node("NumberOfSubgrains")->value());
	}

	if (0 != rootNode->first_node("NumberOfPointsPerSubGrain")) { //resolution
		NumberOfPointsPerSubGrain = std::stoul(
				rootNode->first_node("NumberOfPointsPerSubGrain")->value());
		if ( NumberOfPointsPerSubGrain <= MINIMUM_DISCRETIZATION )
			NumberOfPointsPerSubGrain = MINIMUM_DISCRETIZATION;
	}

	if (0 != rootNode->first_node("NumberOfGridpoints")) {
		if ( std::stoul(rootNode->first_node("NumberOfGridpoints")->value()) > AUTODETECTION_OF_GRIDSIZE ) //user desires explicit setting of the domain size
			NumberOfGridpoints = std::stoul(
			rootNode->first_node("NumberOfGridpoints")->value());
	}

	
	if ( NumberOfGridpoints != AUTODETECTION_OF_GRIDSIZE ) { //adjust resolution per subgrain accordingly if user set hard size of the domain
		if ( PlotDimension == E_3D ) {
			NumberOfPointsPerSubGrain = pow( (((double) (NumberOfGridpoints*NumberOfGridpoints*NumberOfGridpoints) / (double) (NumberOfGrains*NumberOfSubgrains)) / (4.0/3.0 * PI)), (1.0/3.0) );
		}
		if ( PlotDimension == E_2D ) {
			NumberOfPointsPerSubGrain = pow( (((double) (NumberOfGridpoints*NumberOfGridpoints*NumberOfGridpoints) / (double) (NumberOfGrains*NumberOfSubgrains)) / (4.0/3.0 * PI)), (1.0/3.0) );
		}

		if ( NumberOfPointsPerSubGrain <= MINIMUM_DISCRETIZATION ) {
			cout << "Insufficient resolution per subgrain with this domain resolution!" << endl;
			Settings::StatusHealthy = false;
			return;
		}
	}

	//physical quantities and range limitors*****************************************************************************************************************************
	if (0 != rootNode->first_node("StoredElasticEnergyMin")) {
		StoredElasticEnergyMin = std::stod(
				rootNode->first_node("StoredElasticEnergyMin")->value());
		if ( StoredElasticEnergyMin <= 0.0 ) 
			StoredElasticEnergyMin = 0.0;
	}
	if (0 != rootNode->first_node("StoredElasticEnergyMax")) {
		StoredElasticEnergyMax = std::stod(
				rootNode->first_node("StoredElasticEnergyMax")->value());
		if ( StoredElasticEnergyMax <= 0.0 ) 
			StoredElasticEnergyMax = 0.0;
		if ( StoredElasticEnergyMax <= StoredElasticEnergyMin ) 
			StoredElasticEnergyMax = StoredElasticEnergyMin;
	}
	if (0 != rootNode->first_node("StoredElasticScatterGrain")) {
		StoredElasticScatterGrain = std::stod(
				rootNode->first_node("StoredElasticScatterGrain")->value());
		if ( StoredElasticScatterGrain <= 0.0 )
			StoredElasticScatterGrain = 0.0;
	}
	if (0 != rootNode->first_node("StoredElasticScatterSubgrain")) {
		StoredElasticScatterSubgrain = std::stod(
				rootNode->first_node("StoredElasticScatterSubgrain")->value());
		if ( StoredElasticScatterSubgrain <= 0.0 ) 
			StoredElasticScatterSubgrain = 0.0;
	}

	if (0 != rootNode->first_node("DefGrainRelDimensionY")) {
		DefGrainRelDimensionY = std::stod(
				rootNode->first_node("DefGrainRelDimensionY")->value());
		if ( DefGrainRelDimensionY >= 1.0 )
			DefGrainRelDimensionY = 1.0;
	}

	if (0 != rootNode->first_node("DefGrainRelDimensionZ")) {
		DefGrainRelDimensionZ = std::stod(
				rootNode->first_node("DefGrainRelDimensionZ")->value());
		if ( DefGrainRelDimensionZ >= 1.0 )
			DefGrainRelDimensionZ = 1.0;
	}

	if (0 != rootNode->first_node("RandomnessX")) {
		RandomnessX = std::stod(
				rootNode->first_node("RandomnessX")->value());
		if ( RandomnessX >= 0.5 )
			RandomnessX = 0.5;
	}
	if (0 != rootNode->first_node("RandomnessY")) {
		RandomnessY = std::stod(
				rootNode->first_node("RandomnessY")->value());
		if ( RandomnessY >= 0.5 )
			RandomnessY = 0.5;
	}
	if (0 != rootNode->first_node("RandomnessZ")) {
		RandomnessZ = std::stod(
				rootNode->first_node("RandomnessZ")->value());
		if ( RandomnessZ >= 0.5 )
			RandomnessZ = 0.5;
	}

	//*****************************************************************************************************************************

	if (0 != rootNode->first_node("SubgrainOriScatter")) {
		SubgrainOriScatter = std::stod(
				rootNode->first_node("SubgrainOriScatter")->value());
		if ( SubgrainOriScatter <= 0.5 )
			SubgrainOriScatter = 0.5;
	}

	if (0 != rootNode->first_node("VoronoiPeriodic")) {
		VoronoiPeriodic = (bool) std::stoul(rootNode->first_node("VoronoiPeriodic")->value());
	}


	if (0 != rootNode->first_node("ExecuteInParallel")) {
		ExecuteInParallel = (bool) std::stoul(rootNode->first_node("ExecuteInParallel")->value());
	}


	//if (0 != rootNode->first_node("UseOrientationSpace")) {
	//	UseOrientationSpace = (bool) std::stoul(
	//			rootNode->first_node("UseOrientationSpace")->value());
	//}

	//*****************************************************************************************************************************

	if (0 != rootNode->first_node("CrystalStructure")) {
		CrystalStructure = (E_CRYSTAL_STRUCTURE) std::stoi(
				rootNode->first_node("CrystalStructure")->value());

		if (CrystalStructure >= E_DEFAULT_STRUCTURE)
			CrystalStructure = E_DEFAULT_STRUCTURE;
	}

	if (0 != rootNode->first_node("MicroGenMode")) {
		MicroGenMode = (E_MICROGENMODE) std::stoi(
				rootNode->first_node("MicroGenMode")->value());

		if (MicroGenMode >= E_DEFAULT_GENERATOR)
			MicroGenMode = E_DEFAULT_GENERATOR;
	}

	if (0 != rootNode->first_node("TextureSampling")) {
		TextureSampling = (E_SAMPLING) std::stoi(
			rootNode->first_node("TextureSampling")->value());
		
		if (TextureSampling > E_PICK_RANDOMLY)
			TextureSampling = E_DEFAULT_SAMPLING;
	}

	if (0 != rootNode->first_node("TextureGEN")) {
		TextureGEN = (E_TEXTURE) std::stoi(
				rootNode->first_node("TextureGEN")->value());

		if (TextureGEN > E_USE_PREFERENCEORI)
			TextureGEN = E_DEFAULT_TEXTURE;
	}

	if (0 != rootNode->first_node("GrainAggregateType")) {
		GrainAggregation = (E_GRAIN_AGGREGATE) std::stoi( 
			rootNode->first_node("GrainAggregateType")->value());

		if ( GrainAggregation != E_SINGLECRYSTAL && GrainAggregation != E_BICRYSTAL && GrainAggregation != E_POLYCRYSTAL )
			GrainAggregation = E_POLYCRYSTAL;
		if ( (GrainAggregation == E_POLYCRYSTAL && NumberOfGrains < 1) || (GrainAggregation == E_BICRYSTAL && NumberOfGrains != 2) || (GrainAggregation == E_SINGLECRYSTAL && NumberOfGrains != 1) ) {
			cout << "Inconsistent setting of NumberOfGrains and GrainAggregation!" << endl; Settings::StatusHealthy = false; return;
		}
	}

	if (0 != rootNode->first_node("GrainShape")) {
		GrainShape = (E_GRAIN_SHAPE) std::stoi(
				rootNode->first_node("GrainShape")->value());

		if ( GrainShape >= E_DEFAULT_SHAPE)
			GrainShape = E_GLOBULITIC;
	}

	if (0 != rootNode->first_node("PlotIPF2DSection")) {
		PlotIPF2DSection = (bool) std::stoul(
				rootNode->first_node("PlotIPF2DSection")->value());
	}
	if (0 != rootNode->first_node("PlotWindowXMin")) {
			PlotWindowXMin = std::stod(
					rootNode->first_node("PlotWindowXMin")->value());
	}
	if (0 != rootNode->first_node("PlotWindowXMax")) {
				PlotWindowXMax = std::stod(
						rootNode->first_node("PlotWindowXMax")->value());
	}
	if (0 != rootNode->first_node("PlotWindowYMin")) {
				PlotWindowYMin = std::stod(
						rootNode->first_node("PlotWindowYMin")->value());
	}
	if (0 != rootNode->first_node("PlotWindowYMax")) {
				PlotWindowYMax = std::stod(
						rootNode->first_node("PlotWindowYMax")->value());
	}
	if (0 != rootNode->first_node("PlotWindowZMin")) {
				PlotWindowYMax = std::stod(
						rootNode->first_node("PlotWindowZMin")->value());
	}
	if (0 != rootNode->first_node("PlotWindowZMax")) {
				PlotWindowYMax = std::stod(
						rootNode->first_node("PlotWindowZMax")->value());
	}
	if ( PlotWindowXMin <= 0.0 ) PlotWindowXMin = 0.0;
	if ( PlotWindowXMax >= 1.0 ) PlotWindowXMax = 1.0;
	if ( PlotWindowYMin <= 0.0 ) PlotWindowYMin = 0.0;
	if ( PlotWindowYMax >= 1.0 ) PlotWindowYMax = 1.0;
	if ( PlotWindowZMin <= 0.0 ) PlotWindowZMin = 0.0;
	if ( PlotWindowZMax >= 1.0 ) PlotWindowZMax = 1.0;

	if (0 != rootNode->first_node("BreakPeriodicityX")) {
				BreakPerX = (bool) std::stoul(
					rootNode->first_node("BreakPeriodicityX")->value());
	}
	if (0 != rootNode->first_node("BreakPeriodicityY")) {
				BreakPerY = (bool) std::stoul(
					rootNode->first_node("BreakPeriodicityY")->value());
	}
	if (0 != rootNode->first_node("BreakPeriodicityZ")) {
				BreakPerZ = (bool) std::stoul(
					rootNode->first_node("BreakPeriodicityZ")->value());
	}

	if ( (BreakPerX == true || BreakPerY == true || BreakPerZ == true) && (PlotDimension != E_3D) ) {
		cout << "Breaking periodic boundary conditions currently only implemented in 3D!" << endl;
		Settings::StatusHealthy = false; return;
	}

	if ( BreakPerX == true || BreakPerY == true || BreakPerZ == true ) {
		if ( Settings::NumberOfGridpoints < (2 + 2 + 1) ) {
			cout << "Adding 2-voxel thick layers around the domain not possible for breaking boundary conditions because domain is too small!" << endl;
			Settings::StatusHealthy = false; return;
		}
	}

	if (0 != rootNode->first_node("WriteDAMASKOutput")) {
				IO_DAMASK = (bool) std::stoul(
					rootNode->first_node("WriteDAMASKOutput")->value());
	}
	if (0 != rootNode->first_node("WriteHDF5Output")) {
				IO_HDF5 = (bool) std::stoul(
					rootNode->first_node("WriteHDF5Output")->value());
	}


	//*****************************************************************************************************************************

	if (0 != rootNode->first_node("ReadFromFilename")) {
		ReadFromFilename = rootNode->first_node("ReadFromFilename")->value();
	}

	if (0 != rootNode->first_node("AdditionalFilename")) {
		AdditionalFilename =
				rootNode->first_node("AdditionalFilename")->value();
	}

}

