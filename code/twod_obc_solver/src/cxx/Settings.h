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
#ifndef __SETTINGS_H__
#define __SETTINGS_H__

#include <string>
#include "rapidxml.hpp"

/*!
 * \enum E_MICROSTRUCTURE_GEN_MODE
 * \brief Enumeration used to control how the microstructure in the simulation will be
 * generated.
 */	
enum E_MICROSTRUCTURE_GEN_MODE {
	E_READ_VOXELIZED_MICROSTRUCTURE
};
/*!
 * \enum E_CONVOLUTION_MODE
 * \brief Enumeration used to control the convolution kernel.
 */
enum E_CONVOLUTION_MODE {
	E_LAPLACE, E_LAPLACE_RITCHARDSON, E_GAUSSIAN, E_INVALID_VALUE
};

#define E_FCC 		225		//Face-centered cubic
#define E_BCC		229		//Body-centered cubic
#define	E_HCP		194		//Hexagonal

/*!
 * \enum E_RESEARCH_PROJECT
 * \brief Enumeration used to control the research project execution.
 */
enum E_RESEARCH_PROJECT {
	E_DEFAULT_PROJECT
};

enum E_GRAIN_SCHEDULER {
	E_ITERATIVE,
	E_SQUARES,
	E_DEFAULT_SCHEDULER
};
/*!
 * \class Settings
 * \brief Class that holds all global settings that are simulation specific.
 */

struct Magnetic
{
	double VacuumPermeability;
	double MagneticVector_x;
	double MagneticVector_y;
	double MagneticVector_z;
	double deltaMagSys;
	double MagneticForceField;
	//double C_Value;
	//double A_Value;
	//values by C. Mießen
	//https://github.com/GraGLeS/GraGLeS2D/blob/master/params/MagneticField.xml
	Magnetic() :
		VacuumPermeability(12.566e-7),
		MagneticVector_x(0.), MagneticVector_y(0.848048096), MagneticVector_z(0.529919264),
		deltaMagSys(1.18e-5), MagneticForceField(1.35e7) {};
		//C_Value(8.), A_Value(3.0) {};
};


class Settings
{
public:
	static unsigned long StartTime;
	static unsigned int SimulationId;
	static unsigned int NumberOfParticles;
	static unsigned long NumberOfPointsPerGrain;
	static unsigned long NumberOfTimesteps;
	static unsigned long BreakupNumber;
	static unsigned long NetworkExport;
	static unsigned long DomainBorderSize;
	static unsigned long GrainScheduler;
	static E_MICROSTRUCTURE_GEN_MODE MicrostructureGenMode;
	static E_RESEARCH_PROJECT ResearchProject;
	static std::string InputStructureFilename;
	static std::string ConfigFileName;
	static std::string ResultsFileName;
	static unsigned short LatticeType;
	static double HAGB_Mobility;
	static double HAGB_Energy;
	static double Physical_Domain_Size;
	static double TriplePointDrag;
	static double DislocEnPerM;
	static bool IdentifyTwins;
	static bool IsIsotropicNetwork;
	static double MaxMisOrientation;
	static bool ExecuteInParallel;
	static bool GridCoarsement;
	static double GridCoarsementGradient;
	static unsigned long MaximumNumberOfThreads;
	static E_CONVOLUTION_MODE ConvolutionMode;
	static double ConstantSectorRadius;
	static double InterpolatingSectorRadius;
	static bool UseStoredElasticEnergy;
	static bool UseMagneticField;
	static bool DecoupleGrains;
	static Magnetic MagneticParams;
	static double UserDefTimeSlope;

	static bool initializeParameters();
};

#endif	//__SETTINGS_H__
