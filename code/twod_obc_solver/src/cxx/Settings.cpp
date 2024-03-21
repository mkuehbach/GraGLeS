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
#include "Settings.h"
#include "rapidxml.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring>

using namespace std;
using namespace rapidxml;

//Initializing the static setting variables
E_CONVOLUTION_MODE Settings::ConvolutionMode = E_GAUSSIAN;
E_MICROSTRUCTURE_GEN_MODE Settings::MicrostructureGenMode = E_READ_VOXELIZED_MICROSTRUCTURE;
string Settings::InputStructureFilename = "";
string Settings::ConfigFileName = "";
string Settings::ResultsFileName = "";

unsigned long Settings::StartTime = 0;
unsigned int Settings::SimulationId = 0;
unsigned int Settings::NumberOfParticles = 0;
unsigned long Settings::NumberOfPointsPerGrain = 0;
unsigned long Settings::BreakupNumber = 0;
unsigned long Settings::NumberOfTimesteps = 0;
unsigned long Settings::GrainExport = 0;
unsigned long Settings::NetworkExport = 0;
unsigned long Settings::DomainBorderSize = 0;
unsigned long Settings::GrainScheduler = E_SQUARES;
double Settings::GridCoarsementGradient = 0.95;

//physics
unsigned short Settings::LatticeType = 0;
double Settings::HAGB_Energy = 1.;
double Settings::HAGB_Mobility = 1.;
bool Settings::IdentifyTwins = false;
double Settings::Physical_Domain_Size = 0.;
double Settings::DislocEnPerM = 0.;
//double Settings::TriplePointDrag = 0.;
Magnetic Settings::MagneticParams = Magnetic();
double Settings::UserDefTimeSlope = 0.8; //empirical was 0.8359;
bool Settings::IsIsotropicNetwork = false;
bool Settings::GridCoarsement = false;
int Settings::MaxNumberOfOpenMpThreads = 1;
double Settings::ConstantSectorRadius = 0.0;
double Settings::InterpolatingSectorRadius = 0.0;
bool Settings::UseStoredElasticEnergy = false;
bool Settings::UseMagneticField = false;


bool Settings::initializeParameters()
{
	cout << __func__ << "\n";
	ifstream file(Settings::ConfigFileName);
	if (file.fail()) {
		cerr << "Unable to read " << Settings::ConfigFileName << "\n"; return false;
	}

	stringstream contents;
	contents << file.rdbuf();
	string xmlDocument(contents.str());

	xml_document<> tree;
	try {
		tree.parse<0> (&xmlDocument[0]);
	}
	catch (parse_error& Error) {
		cerr << Settings::ConfigFileName << " is not a valid XML!" << "\n"; return false;
	}
	xml_node<>* rootNode = tree.first_node();
	if (0 != strcmp(rootNode->name(), "Parameters")) {
		cerr << "Root node is not named 'Parameters'!" << "\n"; return false;
	}

	string key = "";
	/*
	if (0 != rootNode->first_node(key.c_str())) {
		StartTime = std::stoul(rootNode->first_node(key.c_str())->value());
	}
	if (0 != rootNode->first_node("NumberOfParticles")) {
		NumberOfParticles = std::stoul(
				rootNode->first_node("NumberOfParticles")->value());
	}
	if (0 != rootNode->first_node("NumberOfPointsPerGrain")) {
		NumberOfPointsPerGrain = std::stoul(
				rootNode->first_node("NumberOfPointsPerGrain")->value());
	}
	*/
	key = "GrainExport";
	if (0 != rootNode->first_node(key.c_str())) {
		NetworkExport = std::stoul(rootNode->first_node(key.c_str())->value());
	}

	key = "NumberOfTimesteps";
	if (0 != rootNode->first_node(key.c_str())) {
		NumberOfTimesteps = std::stoul(rootNode->first_node(key.c_str())->value());
	}
	key = "BreakupNumber";
	if (0 != rootNode->first_node(key.c_str())) {
		BreakupNumber = std::stoul(rootNode->first_node(key.c_str())->value());
	}
	key = "DomainBorderSize";
	if (0 != rootNode->first_node(key.c_str())) {
		DomainBorderSize = std::stoul(rootNode->first_node(key.c_str())->value());
	}
	key = "InputStructureFilename";
	if (0 != rootNode->first_node(key.c_str())) {
		InputStructureFilename = rootNode->first_node(key.c_str())->value();
	}
	key = "SpaceGroup";
	if (0 != rootNode->first_node(key.c_str())) {
		unsigned long space_group = std::stoul(rootNode->first_node(key.c_str())->value());
		switch(space_group)
		{
			case 225: { LatticeType = E_FCC; break; }
			case 229: { LatticeType = E_BCC; break; }
			case 194: { LatticeType = E_HCP; break; }
			default: { LatticeType = 0; break; }
		}
	}
	key = "HagbEnergy";
	if (0 != rootNode->first_node(key.c_str())) {
		HAGB_Energy = std::stod(rootNode->first_node(key.c_str())->value());
	}
	key = "HagbMobility";
	if (0 != rootNode->first_node(key.c_str())) {
		HAGB_Mobility = std::stod(rootNode->first_node(key.c_str())->value());
	}
	key = "StrainEnergyDislocations";
	if (0 != rootNode->first_node(key.c_str())) {
		DislocEnPerM = std::stod(rootNode->first_node(key.c_str())->value());
	}
	key = "PhysicalEdgeLengthVe";
	if (0 != rootNode->first_node(key.c_str())) {
		Physical_Domain_Size = std::stod(rootNode->first_node(key.c_str())->value());
	}
	key = "IdentifyTwins";
	if (0 != rootNode->first_node(key.c_str())) {
		IdentifyTwins = (bool) std::stoul(rootNode->first_node(key.c_str())->value());
	}
	key = "UseMagneticField";
	if (0 != rootNode->first_node(key.c_str())) {
		UseMagneticField = (bool) std::stoul(rootNode->first_node(key.c_str())->value());
	}
	/*
	if (0 != rootNode->first_node("MagneticParams")) {
		MagneticParams = rootNode->first_node("MagneticParams")->value();
	}
	*/
	key = "IsIsotropicNetwork";
	if (0 != rootNode->first_node(key.c_str())) {
		IsIsotropicNetwork = (bool) std::stoul(rootNode->first_node(key.c_str())->value());
	}
	key = "MaxNumberOfOpenMpThreads";
	if (0 != rootNode->first_node(key.c_str())) {
		MaxNumberOfOpenMpThreads = std::stoul(rootNode->first_node(key.c_str())->value());
	}
	if (0 != rootNode->first_node("GridCoarsement")) {
		GridCoarsement = (bool) std::stoul(rootNode->first_node("GridCoarsement")->value());
	}
	key = "GridCoarsementGradient";
	if (0 != rootNode->first_node(key.c_str())) {
		GridCoarsementGradient = std::stod(rootNode->first_node(key.c_str())->value());
	}
	ConvolutionMode = E_GAUSSIAN;
	key = "ConstantSectorRadius";
	if (0 != rootNode->first_node(key.c_str())) {
		ConstantSectorRadius = std::stod(rootNode->first_node(key.c_str())->value());
	}
	key = "IpolSectorRadius";
	if (0 != rootNode->first_node(key.c_str())) {
		InterpolatingSectorRadius = std::stod(rootNode->first_node(key.c_str())->value());
	}
	key = "UseStoredElasticEnergy";
	if (0 != rootNode->first_node(key.c_str())) {
		UseStoredElasticEnergy = (bool) std::stoul(rootNode->first_node(key.c_str())->value());
	}
	if (0 != rootNode->first_node("UserDefTimeSlope")) {
		UserDefTimeSlope = std::stod(rootNode->first_node("UserDefTimeSlope")->value());
	}
	file.close();


	return true;
}
