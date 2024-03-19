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

#ifndef GRAINHDL_h
#define GRAINHDL_h

#include "TwodObcSolver.h"
#include "ExpandingVector.h"
#include "Spoint.h"
#include "Settings.h"
#include <omp.h>
#include "Misorientation.h"
#include "GrainScheduler.h"
#include "DimensionalBuffer.h"
#include "../../../utils/src/cxx/PARAPROBE_HDF5Core.h"

#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])
// #define min(x,y) (x<y?x:y)
// #define max(x,y) (x>y?x:y)

#define NX_IDENTIFIER		0	//id
#define NX_SIZE				1	//volume
//double perimeter;
//double GBEnergy;
#define NX_SEE				4	//BulkEnergy
#define NX_ORI				5	//phi1, PHI, phi2
#define NX_BARY				6	//x, y
#define NX_NBOR_CNT			7	//NeighbourCount
#define NX_EDGE_CONTACT		8	//intersectsBoundaryGrain



using namespace voro;
using namespace std;

class LSbox;
class mathMethods;

class Quaternion;

/*!
 * \class grainhdl
 * \brief Class that manages the grain growth simulation.
 */
class grainhdl {
protected:
	int ngrains;
	double dt;
	double h;
	double m_Energy_deltaMAX;
	int realDomainSize;
	int ngridpoints;
	int grid_blowup;

	int Mode;
	GrainScheduler* m_grainScheduler;

public:
	unsigned int currentNrGrains;
	mathMethods* mymath;
	unsigned int loop;
	//! control variable for research mode
	bool loadCurvature;
	unsigned int loadCurvatureLoop;
	bool convolutionCorrection;
	E_RESEARCH_PROJECT project;
	bool constantE;

	//! A 2D vector which stores weights.
	vector<vector<double> > weightsMatrix;

	double ds;
	double delta;
	double *bunge;
	double deviation;
	double BoundaryGrainTube;
	double Realtime;
	double TimeSlope;
	double maxVol;
	MisorientationHdl* m_misOriHdl;
	Quaternion *TwinBoundary;
	DimensionalBuffer<int>* IDField;

	vector<LSbox*> grains;
	LSbox* boundary;

	grainhdl();
	~grainhdl();

	void setResearchAdjustments();
	void setSimulationParameter();
	void read_header_from_nexusfile();
	void read_microstructure_from_nexusfile();
	void find_neighbors();

	void distanceInitialisation();
	void convolution(double& plan_overhead);
	void createConvolutionPlans();
	void destroyConvolutionPlans();
	void comparison_box();

	void updateSecondOrderNeighbors();
	void level_set();
	void redistancing();

	virtual void run_sim();

	vector<unsigned int> get_nexus_grain_identifier();
	vector<double> get_nexus_grain_size();
	vector<double> get_nexus_grain_stored_elastic_energy();
	vector<unsigned char> get_nexus_grain_edge_contact();
	vector<double> get_nexus_grain_orientation();
	vector<double> get_nexus_grain_barycentre();
	void get_nexus_grain_boundary_vertices(vector<double> & vrts);
	vector<size_t> nx_vrts_offsets;
	void get_nexus_grain_boundary_xdmf_topology(vector<unsigned int> & inds );
	void get_nexus_grain_boundary_xdmf_grain_indices( vector<unsigned int> & grain_ids );
	void get_nexus_grain_boundary_info(vector<double> & ifo );
	bool save_NeXus();

	void gridCoarsement();
	void switchDistancebuffer();

	void set_h(double hn);
	void set_realDomainSize(int realDomainSizen);
	//! Used if points are set manually
	void get_biggestGrainVol();
	void find_correctTimestepSize();
	inline LSbox* getGrainByID(unsigned int ID) {
		if (ID == 0)
			return boundary;
		else if (ID > 0 && ID < grains.size())
			return grains[ID];
		else
			return NULL;
	}

	inline long get_ngrains() {
		return ngrains;
	}
	inline int get_realDomainSize() {
		return realDomainSize;
	}
	inline int get_ngridpoints() {
		return ngridpoints;
	}
	inline double get_h() {
		return h;
	}
	inline int get_grid_blowup() {
		return grid_blowup;
	}
	inline int get_loop() {
		return loop;
	}
	inline double get_dt() {
		return dt;
	}
	inline double getBoundaryGrainTube() {
		return BoundaryGrainTube;
	}
	inline double get_ds() {
		return ds;
	}
	inline double get_maxVol() {
		return maxVol;
	}

protected:
	void initEnvironment();
	void initNUMABindings();
	void buildBoxVectors(vector<int> & ID, vector<vector<SPoint>> & contours, vector<double> & q, vector<double> & see );

	int m_ThreadPoolCount;
	vector<ExpandingVector<char> > m_ThreadMemPool;
};
#endif
