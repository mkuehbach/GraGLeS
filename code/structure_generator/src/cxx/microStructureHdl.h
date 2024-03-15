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

#ifndef MICROSTRUCTUREHDL_H_
#define MICROSTRUCTUREHDL_H_

#include <string>
#include "dimensionalBuffer.h"
#include "Eigen/Dense"
#include "myQuaternion.h"
#include "utilities.h"
#include "Settings.h"
#include "../../../utils/src/cxx/PARAPROBE_HDF5Core.h"

using namespace std;
using namespace Eigen;
class Grains;
class SubGrain;
class IterativeGrainScheduler;
class mathMethods;
class myQuaternion;
struct myPreferenceOri;
class randomClass;

#define WHY_TWENTY_SEVEN				27  //TODO::why 27 ?
#define WHY_NINE						9  //TODO::why 9, likely empirical allocation demand see C. Mießen PhD for details


struct RGB
{
	unsigned char R;
	unsigned char G;
	unsigned char B;
	unsigned char ALPHA;
	RGB() : R(0), G(0), B(0), ALPHA(255) {} //white as the default
};


struct BinDiagData {
	double id;							//sub-grain id global view
	double x;							//center of generating tessellation cell
	double y;
	double z;

	double xmi;						//bounding boxes sub-grain
	double xmx;
	double ymi;
	double ymx;

	double zmi;
	double zmx;
	double parentid;				//parent grain id inheriting

	double bunge1;					//Bunge ZXZ orientation grain
	double bunge2;
	double bunge3;
	double volume;					//volume of the grain (area in 2D)

	double stored;					//stored elastic energy of grain
	double dis2parent_meas;			//disorientation angle measured to parent 
	double dis2parent_targ;			//disorientation angle desired on average to parent
	double parent_bunge1;			//parent ori Bunge ZXZ

	double parent_bunge2;
	double parent_bunge3;
	double parent_stored;			//parent stored elastic energy
	double parent_stored_grsc;		//scatter in stored elastic energy for gr

	double parent_stored_sgrsc;		//...for sub-grain
	double parent_size_sc;			//size scaler
	BinDiagData() : id(-1.), x(-1.), y(-1.), z(-1.), xmi(-1.), xmx(-1.), ymi(-1.), ymx(-1.), zmi(-1.), zmx(-1.), parentid(-1.), bunge1(-1.), bunge2(-1.), bunge3(-1.),
		volume(-1.), stored(-1.), dis2parent_meas(-1.), dis2parent_targ(-1.), parent_bunge1(-1.), parent_bunge2(-1.), parent_bunge3(-1.), parent_stored(-1.),
		parent_stored_grsc(-1.), parent_stored_sgrsc(-1.), parent_size_sc(-1.) {}

	void ClearingForSubgrains( void ) { //make sure that again default negative values are shown to identify errors in the dataset as not physical quantity studied here makes sense to be negative
		id = -1.;
		x = -1.;
		y = -1.;
		z = -1.;
		xmi = -1.;
		xmx = -1.;
		ymi = -1.;
		ymx = -1.;
		zmi = -1.;
		zmx = -1.;
		bunge1 = -1.;
		bunge2 = -1.;
		bunge3 = -1.;
		volume = -1.;
		stored = -1.;
		//MK::other data do not need clearing!
	}
	void ClearingComplete( void ) {
		id = -1.;
		x = -1.;
		y = -1.;
		z = -1.;
		xmi = -1.;
		xmx = -1.;
		ymi = -1.;
		ymx = -1.;
		zmi = -1.;
		zmx = -1.;
		parentid = -1.;
		bunge1 = -1.;
		bunge2 = -1.;
		bunge3 = -1.;
		volume = -1.;
		stored = -1.;
		dis2parent_meas = -1.;
		dis2parent_targ = -1.;
		parent_bunge1 = -1.;
		parent_bunge2 = -1.;
		parent_bunge3 = -1.;
		parent_stored = -1.;
		parent_stored_grsc = -1.;
		parent_stored_sgrsc = -1.;
		parent_size_sc = -1.;
	}
};


class timeLogger {

public:
	timeLogger();
	~timeLogger();

	void logev(const string title, double time);
	unsigned int get_entries( void );
	vector <double> times;
	vector <string> titles;

private:
	unsigned int entries;
};



typedef struct {
	int		gid;			//my ID
	int		prid;			//my parent's original ID
	int		x;
	int		y;

	int		z;
	int		xmi;
	int		xmx;
	int		ymi;

	int		ymx;
	int		zmi;
	int		zmx;
	
	double	bunge1;
	double	bunge2;
	double	bunge3;
	double	size;

	double	see;			//my SEE
	double	psee;			//parent's SEE
	double	d2pr_meas;
	
	double	d2pr_targ;
	double	psee_grsc;		//parent's grain SEE scatter
	double	psee_sgrsc;		//parent's subgrain SEE scatter
	double	pszsc;			//parent's size scaler

	double	pbunge1;
	double	pbunge2;
	double	pbunge3;
} SubgrainMetadata;


typedef struct {
	int simid;
	int NGrains;
	int firstgid;
	int DX;
	int DY;
	int DZ;
	int NX;
	int NY;
	int NZ;
} SimMetadata;


/*!
 * \class microStructureHdl
 * \brief Class encapsulating a Level Set Box.
 *
 * LSbox class contains a container of type unsigned it:
 * for each voxel the container stores a grain ID <br>
 *
 * The class is used to organize all routines performed on the hosted grains <br>
 *
 *
 *
 */
class microStructureHdl {

public:
	microStructureHdl();
	virtual ~microStructureHdl();
	void ReadAdditionalInputFiles();
	void ReadPreferenceOrientationFromFile();
	void GeneratePolycrystallineStructureOfGrains();
	void GenerateSubgrainStructureInEachGrain();
	void GeneratorVoro();
	void ReadDiscreteOrientationSet();
	myPreferenceOri FindNextPreferenceOrientation( myQuaternion ori);
	void InitializeGrains(vector<vector<Eigen::Vector3d>> hulls, vector<double> grainVolume);
	void FindNeighbors();
	void ConstructSubgrains();
	void DistributeGrainOriAndSee();
	void DistributeSubgrainOrientations();
	void DistributeSubgrainSee();
	unsigned int CountNumberOfSubgrains();
	void RehashGrainIds();
	void BreakPeriodicity();
	void SaveNeXus();
	void ReportProfile();
	void ReportConfig();
	void InitEnvironment();
	void InitNumaBinding();
	void CopyContainer();

	//getter setter
	unsigned long GetFirstId(){ return first_id; }
	void SetFirstId(bool add, unsigned long val ) {
		if ( add == true ) first_id += val;
	}
	bool healthy;

private:
	DimensionalBuffer<unsigned int>* m_container;
	vector<Grains*> m_grains;
	double m_h;
	unsigned long first_id;
	IterativeGrainScheduler* m_grainScheduler;
	vector<myQuaternion>* m_OrientationSpace;
	vector<Vector3d>* m_ParticlesPositions;
	vector<myPreferenceOri>* m_PreferenceOrientations;
	randomClass* m_seqRND;
	timeLogger myprofiler;
};

#endif /* MICROSTRUCTUREHDL_H_ */
