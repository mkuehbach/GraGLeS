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
#define MINIMUM_VE_SIZE					10
#define AUTODETECTION_OF_GRIDSIZE		0
#define MINIMUM_ASPECT_RATIO			0.001
#define MINIMUM_SUBGRAIN_ORISCATTER		0.5  //degree
#define MINIMUM_SEE						1.e10  //1/m^2

#define MYMAX_NUMA_NODES				64

#define PI								3.14159265358979323846
#define SQR(x)							((x)*(x))
#define CUBE(x)							((x)*(x)*(x))
#define DEG2RAD(deg)					((deg)/(180.)*(PI))
#define RAD2DEG(rad)					((rad)/(PI)*(180.))

//*********************************************************************************************************************************************************************************************

enum E_SAMPLING {
	E_PICK_RANDOMLY
};

enum E_TEXTURE {
	E_USE_PREFERENCEORI
};

enum E_MICROGENMODE {
	E_VORONOI
};

enum E_CRYSTAL_STRUCTURE { //Corresponds to the CRYSTAL STRUCTURE of the sample the user wants to work with
	E_FCC, //Face-centered cubic, assuming space group 225
	E_BCC, //Body-centered cubic, assuming space group 229
	E_HCP  //Hexagonal, assuming space group 194
};

enum E_PLOT_DIMENSION {
	E_0D,
	E_1D,
	E_2D,
	E_3D
};

enum E_GRAIN_AGGREGATE {
	E_SINGLECRYSTAL,
	E_BICRYSTAL,
	E_POLYCRYSTAL
};

enum E_GRAIN_SHAPE { //grains as a Poisson-Voronoi process, subgrains Poisson Voronoi
	E_GLOBULITIC,
	E_FLAT
};

//*********************************************************************************************************************************************************************************************

class Settings {
public:
	static unsigned int SimulationId;
	static unsigned int NumberOfGrains; //Number of wished grains in the sample
	static unsigned int NumberOfSubgrains;
	static unsigned long NumberOfPointsPerSubGrain; //Number of interpolating points for each grain
	static unsigned int NumberOfGridpoints;
	static unsigned long MaximumNumberOfThreads;
	//choose Texture Gen Mode
	static double PhysEdgeLenMatPoint;
	static double SubgrainOriScatter;
	static double StoredElasticEnergyMin;
	static double StoredElasticEnergyMax; //Energy of a given dislocation in the sample
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
	static bool StatusHealthy;	//global indication of successful and consistent parameter file
	static bool BreakPerX;
	static bool BreakPerY;
	static bool BreakPerZ;

	//static bool UseOrientationSpace;
	static E_MICROGENMODE MicroGenMode;
	static E_CRYSTAL_STRUCTURE CrystalStructure;
	static E_SAMPLING TextureSampling;
	static E_TEXTURE TextureGen;
	static E_PLOT_DIMENSION Dimensionality;
	static E_GRAIN_AGGREGATE GrainAggregation;
	static E_GRAIN_SHAPE GrainShape;

	static string ConfigFileName;
	static string ReadFromFilename;
	static string AdditionalFilename;
	static string ResultsFileName;

	static void ReadXmlConfig(string filename = "");

};

#endif /* SETTINGS_H_ */
