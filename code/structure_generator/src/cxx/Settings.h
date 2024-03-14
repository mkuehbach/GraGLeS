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
#ifndef SETTINGS_H_
#define SETTINGS_H_
using namespace std;
#include <string>

#define MINIMUM_DISCRETIZATION			8
#define AUTODETECTION_OF_GRIDSIZE		0
#define PI								3.14159265358979323846

//*********************************************************************************************************************************************************************************************

enum E_SAMPLING {
	E_PICK_RANDOMLY,
	E_DEFAULT_SAMPLING
};

enum E_TEXTURE {
	E_DEFAULT_TEXTURE,
	E_USE_PREFERENCEORI
};

enum E_MICROGENMODE {
	E_VORONOI,
	E_PREDEFINEDSTRUCTURE,
	E_DEFAULT_GENERATOR
};

enum E_CRYSTAL_STRUCTURE { //Corresponds to the CRYSTAL STRUCTURE of the sample the user wants to work with
	E_HCP, //Hexagonal compact
	E_FCC, //Face-centered cubic
	E_BCC, //Body-centered cubic
	E_DEFAULT_STRUCTURE
};

enum E_PLOT_DIMENSION {
	E_DEFAULT_DIMENSION,
	E_1D,
	E_2D,
	E_3D
};

enum E_GRAIN_AGGREGATE {
	E_SINGLECRYSTAL,
	E_BICRYSTAL,
	E_POLYCRYSTAL,
};

enum E_GRAIN_SHAPE { //grains as a Poisson-Voronoi process, subgrains Poisson Voronoi
	E_GLOBULITIC,
	E_FLAT,
	E_DEFAULT_SHAPE
};

//*********************************************************************************************************************************************************************************************

class Settings {
public:

	static unsigned int SimID;
	static unsigned int NumberOfGrains; //Number of wished grains in the sample
	static unsigned int NumberOfSubgrains;
	static unsigned int NumberOfGridpoints;
	static unsigned long NumberOfPointsPerSubGrain; //Number of interpolating points for each grain
	static unsigned long MaximumNumberOfThreads;
	//choose Texture Gen Mode
	static double SubgrainOriScatter;
	static double StoredElasticEnergyMax; //Energy of a given dislocation in the sample
	static double StoredElasticEnergyMin;
	static double StoredElasticScatterGrain;
	static double StoredElasticScatterSubgrain;
	static double DefGrainRelDimensionX;
	static double DefGrainRelDimensionY;
	static double DefGrainRelDimensionZ;
	static double RandomnessX;
	static double RandomnessY;
	static double RandomnessZ;

	static bool VoronoiPeriodic;
	static bool ExecuteInParallel; //Activate the execution in parallel
	static bool PlotIPF2DSection;
	static bool StatusHealthy;	//global indication of successful and consistent parameter file
	static double PlotWindowXMin;
	static double PlotWindowXMax;
	static double PlotWindowYMin;
	static double PlotWindowYMax;
	static double PlotWindowZMin;
	static double PlotWindowZMax;

	static bool BreakPerX;
	static bool BreakPerY;
	static bool BreakPerZ;
	static bool IO_DAMASK;
	static bool IO_HDF5;

	//static bool UseOrientationSpace;
	static E_MICROGENMODE MicroGenMode;
	static E_CRYSTAL_STRUCTURE CrystalStructure;
	static E_SAMPLING TextureSampling;
	static E_TEXTURE TextureGEN;
	static E_PLOT_DIMENSION PlotDimension;
	static E_GRAIN_AGGREGATE GrainAggregation;
	static E_GRAIN_SHAPE GrainShape;

	static string ReadFromFilename;
	static string AdditionalFilename;
	static string ResultsFileName;

	static void readXML(string filename = "");

};

#endif /* SETTINGS_H_ */
