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
unsigned long Settings::StartTime = 0;
unsigned int Settings::SimulationId = 0;
unsigned int Settings::NumberOfParticles = 0;
unsigned long Settings::NumberOfPointsPerGrain = 0;
unsigned long Settings::BreakupNumber = 0;
unsigned long Settings::NumberOfTimesteps = 0;
unsigned long Settings::AnalysisTimestep = 0;
unsigned long Settings::NetworkExport = 0;
unsigned long Settings::DiscreteSamplingRate = 0;
unsigned long Settings::DomainBorderSize = 0;
unsigned long Settings::MaximumNumberOfThreads = 0;
unsigned long Settings::GrainScheduler = 0;
double Settings::GridCoarsementGradient = 0.95;
E_CONVOLUTION_MODE Settings::ConvolutionMode = E_INVALID_VALUE;
E_MICROSTRUCTURE_GEN_MODE Settings::MicrostructureGenMode = E_INVALID_VAL;
E_RESEARCH_PROJECT Settings::ResearchProject = E_DEFAULT;
string Settings::ReadFromFilename = "";
string Settings::AdditionalFilename = "";
string Settings::ConfigFileName = "";
string Settings::ResultsFileName = "";
unsigned long Settings::LatticeType = 0;
double Settings::HAGB_Energy = 1.0;
double Settings::HAGB_Mobility = 1.0;
double Settings::Physical_Domain_Size = 0.0;
double Settings::DislocEnPerM = 0.0;
double Settings::TriplePointDrag = 0.0;
unsigned long Settings::UseMobilityModel = 0;
bool Settings::IdentifyTwins = 0;
bool Settings::IsIsotropicNetwork = false;
bool Settings::UseTexture = false;
double Settings::MaxMisOrientation = 0.0;
bool Settings::ExecuteInParallel = false;
bool Settings::GridCoarsement = false;
bool Settings::DecoupleGrains = false;
//! Enables the generation and capture of a wider spectrum of analysing data
bool Settings::ResearchMode = false;
double Settings::ConstantSectorRadius = 0.0;
double Settings::InterpolatingSectorRadius = 0.0;
unsigned long Settings::NeighbourTracking = 0;
bool Settings::UseStoredElasticEnergy = false;
bool Settings::UseMagneticField = false;
Magnetic Settings::MagneticParams = Magnetic();
double Settings::UserDefTimeSlope = 0.8; //empirical was 0.8359;


bool Settings::initializeParameters(string filename)
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
	try {
		tree.parse<0> (&xmlDocument[0]);
	} catch (parse_error& Error) {
		cout << filename.c_str() << "is not a valid XML!" << endl;
		cout << "Exception is " << Error.what() << endl;
	}
	xml_node<>* rootNode = tree.first_node();
	if (0 != strcmp(rootNode->name(), "Parameters")) {
		cout << "Root node is not 'Parameters'! Will now halt!" << endl;
		exit(2);
	}

	//Now read all parameters if present
	if (0 != rootNode->first_node("StartTime")) {
		StartTime = std::stoul(rootNode->first_node("StartTime")->value());
	}
	if (0 != rootNode->first_node("NumberOfParticles")) {
		NumberOfParticles = std::stoul(
				rootNode->first_node("NumberOfParticles")->value());
	}
	/*
	if (0 != rootNode->first_node("NumberOfPointsPerGrain")) {
		NumberOfPointsPerGrain = std::stoul(
				rootNode->first_node("NumberOfPointsPerGrain")->value());
	}
	*/
	if (0 != rootNode->first_node("AnalysisTimestep")) {
		AnalysisTimestep = std::stoul(
				rootNode->first_node("AnalysisTimestep")->value());
	}
	if (0 != rootNode->first_node("PlotInterval")) {
		NetworkExport
				= std::stoul(rootNode->first_node("PlotInterval")->value());
	}
	if (0 != rootNode->first_node("NumberOfTimesteps")) {
		NumberOfTimesteps = std::stoul(
				rootNode->first_node("NumberOfTimesteps")->value());
	}
	if (0 != rootNode->first_node("BreakupNumber")) {
		BreakupNumber = std::stoul(
				rootNode->first_node("BreakupNumber")->value());
	}
	if (0 != rootNode->first_node("DiscreteSamplingRate")) {
		DiscreteSamplingRate = std::stoul(
				rootNode->first_node("DiscreteSamplingRate")->value());
	}
	if (0 != rootNode->first_node("DomainBorderSize")) {
		DomainBorderSize = std::stoul(
				rootNode->first_node("DomainBorderSize")->value());
	}
	if (0 != rootNode->first_node("MicrostructureGenMode")) {
		MicrostructureGenMode = (E_MICROSTRUCTURE_GEN_MODE) std::stoi(
				rootNode->first_node("MicrostructureGenMode")->value());
		if (MicrostructureGenMode >= E_INVALID_VAL)
			MicrostructureGenMode = E_INVALID_VAL;
	}
	if (0 != rootNode->first_node("ReadFromFilename")) {
		ReadFromFilename = rootNode->first_node("ReadFromFilename")->value();
	}
	if (0 != rootNode->first_node("AdditionalFilename")) {
		AdditionalFilename
				= rootNode->first_node("AdditionalFilename")->value();
	}
	if (0 != rootNode->first_node("LatticeType")) {
		LatticeType = (E_LATTICE_TYPE) std::stoi(
				rootNode->first_node("LatticeType")->value());
		if (LatticeType >= E_INVALID_LATTICE)
			LatticeType = E_INVALID_LATTICE;
	}
	if (0 != rootNode->first_node("HAGB_Energy")) {
		HAGB_Energy = std::stod(rootNode->first_node("HAGB_Energy")->value());
	}
	if (0 != rootNode->first_node("HAGB_Mobility")) {
		HAGB_Mobility = std::stod(
				rootNode->first_node("HAGB_Mobility")->value());
	}
	if (0 != rootNode->first_node("DislocEnPerM")) {
		DislocEnPerM = std::stod(rootNode->first_node("DislocEnPerM")->value());
	}
	if (0 != rootNode->first_node("Physical_Domain_Size")) {
		Physical_Domain_Size = std::stod(
				rootNode->first_node("Physical_Domain_Size")->value());
	}
	if (0 != rootNode->first_node("TriplePointDrag")) {
		TriplePointDrag = std::stod(
				rootNode->first_node("TriplePointDrag")->value());
	}
	if (0 != rootNode->first_node("UseMobilityModel")) {
		UseMobilityModel = std::stoul(
				rootNode->first_node("UseMobilityModel")->value());
	}
	if (0 != rootNode->first_node("IdentifyTwins")) {
		IdentifyTwins = (bool) std::stoul(
				rootNode->first_node("IdentifyTwins")->value());
	}
	if (0 != rootNode->first_node("UseMagneticField")) {
		UseMagneticField = (bool) std::stoul(
				rootNode->first_node("UseMagneticField")->value());
	}
	if (0 != rootNode->first_node("MagneticParams")) {
		MagneticParams = rootNode->first_node("MagneticParams")->value();
	}
	if (0 != rootNode->first_node("IsIsotropicNetwork")) {
		IsIsotropicNetwork = (bool) std::stoul(
				rootNode->first_node("IsIsotropicNetwork")->value());
	}
	if (0 != rootNode->first_node("UseTexture")) {
		UseTexture = (bool) std::stoul(
				rootNode->first_node("UseTexture")->value());
	}
	if (0 != rootNode->first_node("MaxMisOrientation")) {
		MaxMisOrientation = std::stoul(
				rootNode->first_node("MaxMisOrientation")->value());
	}
	if (0 != rootNode->first_node("ExecuteInParallel")) {
		ExecuteInParallel = (bool) std::stoul(
				rootNode->first_node("ExecuteInParallel")->value());
	}
	if (0 != rootNode->first_node("MaximumNumberOfThreads")) {
		MaximumNumberOfThreads = std::stoul(
				rootNode->first_node("MaximumNumberOfThreads")->value());

	}
	if (0 != rootNode->first_node("GridCoarsement")) {
		GridCoarsement = std::stoul(
				rootNode->first_node("GridCoarsement")->value());
	}
	if (0 != rootNode->first_node("GridCoarsementGradient")) {
		GridCoarsementGradient = std::stod(
				rootNode->first_node("GridCoarsementGradient")->value());
	}
	if (0 != rootNode->first_node("ConvolutionMode")) {
		ConvolutionMode = (E_CONVOLUTION_MODE) std::stoi(
				rootNode->first_node("ConvolutionMode")->value());
		if (ConvolutionMode >= E_INVALID_VALUE)
			ConvolutionMode = E_INVALID_VALUE;
	}
	//!
	if (0 != rootNode->first_node("ResearchMode")) {
		ResearchMode
				= std::stoul(rootNode->first_node("ResearchMode")->value());
	}
	//! Research modification
	if (0 != rootNode->first_node("ResearchProject")) {
		ResearchProject = (E_RESEARCH_PROJECT) std::stoi(
				rootNode->first_node("ResearchProject")->value());
		if (ResearchProject >= E_DEFAULT)
			ResearchProject = E_DEFAULT;
	}
	if (0 != rootNode->first_node("ConstantSectorRadius")) {
		ConstantSectorRadius = std::stod(rootNode->first_node("ConstantSectorRadius")->value());
	}
	if (0 != rootNode->first_node("InterpolatingSectorRadius")) {
		InterpolatingSectorRadius = std::stod(rootNode->first_node("InterpolatingSectorRadius")->value());
	}
	if (0 != rootNode->first_node("NeighbourTracking")) {
		NeighbourTracking = std::stoul(rootNode->first_node("NeighbourTracking")->value());
	}

	if (0 != rootNode->first_node("UseStoredElasticEnergy")) {
		UseStoredElasticEnergy = (bool) std::stoul(rootNode->first_node("UseStoredElasticEnergy")->value());
	}
	if (0 != rootNode->first_node("GrainScheduler")) {
		GrainScheduler = (E_GRAIN_SCHEDULER) std::stoi(rootNode->first_node("GrainScheduler")->value());
		if (GrainScheduler >= E_DEFAULT_SCHEDULER)
			GrainScheduler = E_DEFAULT_SCHEDULER;
	}
	if (0 != rootNode->first_node("UserDefTimeSlope")) {
		UserDefTimeSlope = std::stod(rootNode->first_node("UserDefTimeSlope")->value());
	}
	file.close();


	return true;
}
