/*
	IMM_MicrostructureGenerator
	A program to instantiate two-staged Poisson-Voronoi tessellation microstructures of 
	parent grains and their sub-grains with adjustable properties such as orientation, and dislocation density
	Copyright (C) 2016
	Christian Miessen (data structures), Markus K�hbach (physical metallurgy functionalities, PRNGs), 
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

unsigned int Settings::SimulationId = 0; //default naming
unsigned int Settings::NumberOfGrains = 0;
unsigned int Settings::NumberOfSubgrains = 0;
unsigned long Settings::NumberOfPointsPerSubGrain = 0;
unsigned int Settings::NumberOfGridpoints = 0;
unsigned long Settings::MaximumNumberOfThreads = 1;
double Settings::PhysEdgeLenMatPoint = 0.;
double Settings::SubgrainOriScatter = 0.;
double Settings::StoredElasticEnergyMax = 0.;
double Settings::StoredElasticEnergyMin = 0.;
double Settings::StoredElasticScatterGrain = 0.;
double Settings::StoredElasticScatterSubgrain = 0.;
double Settings::DefGrainRelDimensionX = 1.; //by definition
double Settings::DefGrainRelDimensionY = 1.;
double Settings::DefGrainRelDimensionZ = 1.;
double Settings::RandomnessX = 0.;
double Settings::RandomnessY = 0.;
double Settings::RandomnessZ = 0.;
bool Settings::VoronoiPeriodic = false;
bool Settings::ExecuteInParallel = false;
bool Settings::StatusHealthy = true;
bool Settings::BreakPerX = false; //by default domains are periodic
bool Settings::BreakPerY = false;
bool Settings::BreakPerZ = false;
E_CRYSTAL_STRUCTURE Settings::CrystalStructure = E_FCC;
E_MICROGENMODE Settings::MicroGenMode = E_VORONOI;
E_SAMPLING Settings::TextureSampling = E_PICK_RANDOMLY;
E_TEXTURE Settings::TextureGen = E_USE_PREFERENCEORI;
E_PLOT_DIMENSION Settings::Dimensionality = E_3D;
E_GRAIN_AGGREGATE Settings::GrainAggregation = E_POLYCRYSTAL;
E_GRAIN_SHAPE Settings::GrainShape = E_GLOBULITIC;
string Settings::ConfigFileName = "";
string Settings::ReadFromFilename = "";
string Settings::AdditionalFilename = "";
string Settings::ResultsFileName = "";

//Definition of the needed functions **********************************************************************************************************************************************************

void Settings::ReadXmlConfig(string filename)
{
	cout << __func__ << "\n";
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
	string key = "";
	//Read all the required parameters ****************************************************************************************************************************************************
	key = "SimulationID";
	if (0 != rootNode->first_node(key.c_str())) {
		SimulationId = std::stoul(rootNode->first_node(key.c_str())->value());
	}
	key = "Dimensionality";
	if (0 != rootNode->first_node(key.c_str())) {
		Dimensionality = (E_PLOT_DIMENSION) std::stoi(rootNode->first_node(key.c_str())->value());
		if ( Dimensionality != E_2D && Dimensionality != E_3D ) {
			cerr << key << " needs to be 2 or 3 !" << "\n";
			StatusHealthy = false; return;
		}
	}
	key = "NumberOfGrains";
	if (0 != rootNode->first_node(key.c_str())) {
		NumberOfGrains = std::stoul(rootNode->first_node(key.c_str())->value());
		if (NumberOfGrains < 1) {
			cerr << "WARNING::Auto-resetting " << key << " to 1" << "\n";
			NumberOfGrains = 1;
		}
	}
	key = "NumberOfSubgrains";
	if (0 != rootNode->first_node(key.c_str())) {
		NumberOfSubgrains = std::stoul(rootNode->first_node(key.c_str())->value());
		if (NumberOfSubgrains < 1) {
			cerr << "WARNING::Auto-resetting " << key << " to 1" << "\n";
			NumberOfSubgrains = 1;
		}
	}
	key = "NumberOfPointsPerSubGrain";
	if (0 != rootNode->first_node(key.c_str())) { //average radius argument
		NumberOfPointsPerSubGrain = std::stoul(rootNode->first_node(key.c_str())->value());
		if ( NumberOfPointsPerSubGrain <= MINIMUM_DISCRETIZATION ) {
			cerr << "WARNING::Auto-resetting " << key << " to " << MINIMUM_DISCRETIZATION << "\n";
			NumberOfPointsPerSubGrain = MINIMUM_DISCRETIZATION;
		}
	}
	key = "PhysicalEdgeLengthPerMaterialPoint";
	if (0 != rootNode->first_node(key.c_str())) {
		PhysEdgeLenMatPoint = std::stod(rootNode->first_node(key.c_str())->value());
	}
	//NumberOfGridpoints will be defined later

	//physical quantities and range limitors*****************************************************************************************************************************
	key = "CrystalStructure";
	if (0 != rootNode->first_node(key.c_str())) {
		CrystalStructure = (E_CRYSTAL_STRUCTURE) std::stoi(rootNode->first_node(key.c_str())->value());
		if ( CrystalStructure != E_FCC && CrystalStructure != E_BCC && CrystalStructure != E_HCP ) {
			cerr << "Unsupported " << key << "\n";
			StatusHealthy = false; return;
		}
	}
	key = "StoredElasticEnergyMin";
	if (0 != rootNode->first_node(key.c_str())) {
		StoredElasticEnergyMin = std::stod(rootNode->first_node(key.c_str())->value());
		if ( StoredElasticEnergyMin <= __DBL_EPSILON__ ) { 
			StoredElasticEnergyMin = 0.;
		}
	}
	key = "StoredElasticEnergyMax";
	if (0 != rootNode->first_node(key.c_str())) {
		StoredElasticEnergyMax = std::stod(rootNode->first_node(key.c_str())->value());
		if ( StoredElasticEnergyMax <= __DBL_EPSILON__ ) {
			StoredElasticEnergyMax = 0.;
		}
		if ( StoredElasticEnergyMax <= StoredElasticEnergyMin ) {
			StoredElasticEnergyMax = StoredElasticEnergyMin;
		}
	}
	key = "StoredElasticEnergyScatterGrain";
	if (0 != rootNode->first_node(key.c_str())) {
		StoredElasticScatterGrain = std::stod(rootNode->first_node(key.c_str())->value());
		if ( StoredElasticScatterGrain <= __DBL_EPSILON__ ) {
			StoredElasticScatterGrain = 0.;
		}
	}
	key = "StoredElasticEnergyScatterSubgrain";
	if (0 != rootNode->first_node(key.c_str())) {
		StoredElasticScatterSubgrain = std::stod(rootNode->first_node(key.c_str())->value());
		if ( StoredElasticScatterSubgrain <= __DBL_EPSILON__ ) {
			StoredElasticScatterSubgrain = 0.;
		}
	}
	key = "OrientationScatterSubgrain";
	if (0 != rootNode->first_node(key.c_str())) {
		SubgrainOriScatter = std::stod(rootNode->first_node(key.c_str())->value());
		if ( SubgrainOriScatter <= MINIMUM_SUBGRAIN_ORISCATTER ) {
			SubgrainOriScatter = MINIMUM_SUBGRAIN_ORISCATTER;
		}
	}
	key = "DefGrainRelDimensionY";
	if (0 != rootNode->first_node(key.c_str())) {
		DefGrainRelDimensionY = std::stod(rootNode->first_node(key.c_str())->value());
		if ( DefGrainRelDimensionY <= MINIMUM_ASPECT_RATIO ) {
			DefGrainRelDimensionY = MINIMUM_ASPECT_RATIO;
		}
		if ( DefGrainRelDimensionY >= (1. - __DBL_EPSILON__) ) {
			DefGrainRelDimensionY = 1.;
		}
	}
	key = "DefGrainRelDimensionZ";
	if (0 != rootNode->first_node(key.c_str())) {
		DefGrainRelDimensionZ = std::stod(rootNode->first_node(key.c_str())->value());
		if ( DefGrainRelDimensionZ <= MINIMUM_ASPECT_RATIO ) {
			DefGrainRelDimensionZ = MINIMUM_ASPECT_RATIO;
		}
		if ( DefGrainRelDimensionZ >= (1. - __DBL_EPSILON__) ) {
			DefGrainRelDimensionZ = 1.;
		}
	}
	key = "RandomnessX";
	if (0 != rootNode->first_node(key.c_str())) {
		RandomnessX = std::stod(rootNode->first_node(key.c_str())->value());
		if ( RandomnessX <= __DBL_EPSILON__ ) {
			RandomnessX = 0.;
		}
		if ( RandomnessX >= (0.5 - __DBL_EPSILON__) ) {
			RandomnessX = 0.5;
		}
	}
	key = "RandomnessY";
	if (0 != rootNode->first_node(key.c_str())) {
		RandomnessY = std::stod(rootNode->first_node(key.c_str())->value());
		if ( RandomnessY <= __DBL_EPSILON__ ) {
			RandomnessY = 0.;
		}
		if ( RandomnessY >= (0.5 - __DBL_EPSILON__) ) {
			RandomnessY = 0.5;
		}
	}
	key = "RandomnessZ";
	if (0 != rootNode->first_node(key.c_str())) {
		RandomnessZ = std::stod(rootNode->first_node(key.c_str())->value());
		if ( RandomnessZ <= __DBL_EPSILON__ ){
			RandomnessZ = 0.;
		}
		if ( RandomnessZ >= (0.5 - __DBL_EPSILON__) ) {
			RandomnessZ = 0.5;
		}
	}
	key = "UseMultithreading";
	if (0 != rootNode->first_node(key.c_str())) {
		ExecuteInParallel = (bool) std::stoul(rootNode->first_node(key.c_str())->value());
	}
	MicroGenMode = E_VORONOI;
	TextureSampling = E_PICK_RANDOMLY;
	TextureGen = E_USE_PREFERENCEORI;
	key = "GrainShape";
	if (0 != rootNode->first_node(key.c_str())) {
		GrainShape = (E_GRAIN_SHAPE) std::stoi(rootNode->first_node(key.c_str())->value());
		if ( GrainShape != E_GLOBULITIC && GrainShape != E_FLAT ) {
			cerr << "Unsupported GrainShape!" << "\n";
			StatusHealthy = false; return;
		}
	}
	key = "BreakPeriodicityX";
	if (0 != rootNode->first_node(key.c_str())) {
		BreakPerX = (bool) std::stoul(rootNode->first_node(key.c_str())->value());
	}
	key = "BreakPeriodicityY";
	if (0 != rootNode->first_node(key.c_str())) {
		BreakPerY = (bool) std::stoul(rootNode->first_node(key.c_str())->value());
	}
	key = "BreakPeriodicityZ";
	if (0 != rootNode->first_node(key.c_str())) {
		BreakPerZ = (bool) std::stoul(rootNode->first_node(key.c_str())->value());
	}
	if ( (BreakPerX == true || BreakPerY == true || BreakPerZ == true) && (Dimensionality != E_3D) ) {
		cerr << "Breaking periodic boundary conditions currently only implemented in 3D !" << "\n";
		StatusHealthy = false; return;
	}
	if ( BreakPerX == true || BreakPerY == true || BreakPerZ == true ) {
		if ( Settings::NumberOfGridpoints < (2 + 2 + 1) ) {
			cerr << "Adding 2-voxel thick layers around the domain not possible for breaking boundary conditions because domain is too small !" << "\n";
			StatusHealthy = false; return;
		}
	}
	//*****************************************************************************************************************************
	key = "ReadFromFilename";
	if (0 != rootNode->first_node(key.c_str())) {
		ReadFromFilename = rootNode->first_node(key.c_str())->value();
	}
	key = "AdditionalFilename";
	if (0 != rootNode->first_node(key.c_str())) {
		AdditionalFilename = rootNode->first_node(key.c_str())->value();
	}
}
