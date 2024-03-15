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

#include "Grains.h"
#include "microStructureHdl.h"
#include "../../../thirdparty/mandatory/voroxx/voro/src/voro++.hh"
#include "utilities.h"
#include "Settings.h"
#include "SubGrain.h"
#include <stdexcept>
#include <iostream>
#include "random.h"

using namespace std;
using namespace voro;
using namespace Eigen;

Grains::Grains()
{
}


Grains::Grains(int id, vector<Vector3d>& hull, DimensionalBuffer<unsigned int>* container,
	double volume) : m_Volume(volume), m_ID(id), m_oldID(id)
{
	m_orientation = NULL;
	int xmin = Settings::NumberOfGridpoints, xmax = 0, ymin =
			Settings::NumberOfGridpoints, ymax = 0, zmin =
			Settings::NumberOfGridpoints, zmax = 0;
	double h = 1. / Settings::NumberOfGridpoints;
	for (unsigned int j = 0; j < hull.size(); j++) {
//		cout << hull[j][0] / h << "  " << hull[j][1] / h << "  "
//				<< hull[j][2] / h << endl;

		double test = hull[j][0] / h;
		if (test < xmin)
			xmin = int(test);
		if (test > xmax)
			xmax = int(test);
		test = hull[j][1] / h;
		if (test < ymin)
			ymin = int(test);
		if (test > ymax)
			ymax = int(test);
		test = hull[j][2] / h;
		if (test < zmin)
			zmin = int(test);
		if (test > zmax)
			zmax = int(test);
	}
	if (xmax < (Settings::NumberOfGridpoints - 2))
		xmax = xmax + 1;
	if (ymax < (Settings::NumberOfGridpoints - 2))
		ymax = ymax + 1;
	if (zmax < (Settings::NumberOfGridpoints - 2))
		zmax = zmax + 1;
	if (Settings::Dimensionality == E_2D) {
		zmin = 0;
		zmax = 1;
	}

	m_orientation = NULL;
	m_PrefOri = NULL;

	m_localContainer = new DimensionalBuffer<unsigned int>(ymin, xmin, zmin,
			ymax, xmax, zmax);
}


void Grains::copyContainerToGrain(DimensionalBuffer<unsigned int>* container)
{
	copyVoxelData(container);
}


Grains::Grains(myQuaternion ori, double SEE) : m_SEE(SEE) {
	m_orientation = NULL;
	set_Orientation(ori);
}


Grains::~Grains()
{
	delete m_localContainer;
	if (m_orientation != NULL) {
		delete m_orientation;
	}
	if(m_PrefOri != NULL ) {
		delete m_PrefOri;
	}
	for ( vector<SubGrain*>::iterator itS = m_SubGrains.begin(); itS != m_SubGrains.end(); itS++ ) {
		delete *itS;
	}
}

// Erstellung des neuen lokalen Containers = Umbox --> um Subkörner zu erhalten

void Grains::set_Orientation(myQuaternion ori) {
	if (m_orientation == NULL)
		m_orientation = new myQuaternion();
	*m_orientation = ori;
}

void Grains::computeDirectNeighbours(
		const RTree<unsigned int, int, 3, float>& tree) {
	int min[3], max[3];
	min[0] = getMinX();
	min[1] = getMinY();
	min[2] = getMinZ();
	max[0] = getMaxX();
	max[1] = getMaxY();
	max[2] = getMaxZ();
	vector<unsigned int> intersectingGrains;
	tree.Search(min, max, intersectingGrains);
	for (unsigned int k = 0; k < intersectingGrains.size(); k++) {
		if (m_ID != intersectingGrains[k]) {
			m_NeighborCandidates.push_back(intersectingGrains[k]);
		}
	}
}

void Grains::copyVoxelData(DimensionalBuffer<unsigned int> *container) {
	for (int k = getMinZ(); k < getMaxZ(); k++)
		for (int j = getMinY(); j < getMaxY(); j++)
			for (int i = getMinX(); i < getMaxX(); i++) {
				unsigned int value = container->getValueAt(j, i, k);
				m_localContainer->setValueAt(j, i, k, value);
			}
}


int Grains::copySubgrainsToGlobalContainer(
		DimensionalBuffer<unsigned int> *container, int offset) {
	unsigned int newOffset = 0;
	if (m_SubGrains.size() > 1) {
		newOffset = rehashAllSubGrains(offset);
		for (int k = getMinZ(); k < getMaxZ(); k++)
			for (int j = getMinY(); j < getMaxY(); j++) {
				for (int i = getMinX(); i < getMaxX(); i++)
					if (m_localContainer->getValueAt(j, i, k) > 0)
						container->setValueAt(j, i, k,
								(unsigned int) offset
										+ m_localContainer->getValueAt(j, i,
												k));
			}
	} else {
		newOffset = offset + 1;
		for (int k = getMinZ(); k < getMaxZ(); k++)
			for (int j = getMinY(); j < getMaxY(); j++) {
				for (int i = getMinX(); i < getMaxX(); i++)
					if (m_localContainer->getValueAt(j, i, k) == m_ID)
						container->setValueAt(j, i, k, (unsigned int) newOffset);
			}
		m_ID = newOffset; //but the m_oldID is not overwritten to allow identifying the parent grain
	}
	return newOffset;
}

/*
void Grains::correctSubgrainVolWhenBreakingPer() {
	if ( Settings::BreakPerX == true ) {
		if (m_SubGrains.size() > 1) {
			if ( getMinX() >= 0 ) { //boundary contact with AIR_GRAIN_XM
				for (int k = getMinZ(); k < getMaxZ(); k++)
					for (int j = getMinY(); j < getMaxY(); j++) {
						//if (m_localContainer->getValueAt(j, 0, k) > 0) // erode subgrain volume by this voxel as it will belong to AIR_GRAIN_XM
					}
			}
			if ( getMinX() >= 1 ) {
				for (int k = getMinZ(); k < getMaxZ(); k++)
					for (int j = getMinY(); j < getMaxY(); j++) {
						//if (m_localContainer->getValueAt(j, 0, k) > 0) // substract a voxel from that subgrain
					}
			}
			//AIR_GRAIN_XP
			//###.....
	} else {
		newOffset = offset + 1;
		for (int k = getMinZ(); k < getMaxZ(); k++)
			for (int j = getMinY(); j < getMaxY(); j++) {
				for (int i = getMinX(); i < getMaxX(); i++)
					if (m_localContainer->getValueAt(j, i, k) == m_ID)
						container->setValueAt(j, i, k, (unsigned int) newOffset);
			}
		m_ID = newOffset; //but the m_oldID is not overwritten to allow identifying the parent grain
	}
	return newOffset;
}
*/


int Grains::rehashAllSubGrains(int offset) {
	vector<SubGrain*>::iterator itG;
	if (Settings::NumberOfSubgrains != 0 && m_SubGrains.size() > 1)
		for (itG = ++m_SubGrains.begin(); itG != m_SubGrains.end(); itG++) {
			if (*itG != NULL)
				(*itG)->setID(offset + (*itG)->get_ID());
		}
	return offset + m_SubGrains.size();
}


void Grains::generateSubGrainOri(randomClass& r) {
	
	double deviation = m_PrefOri->subgrainsScatterOri;
	if (Settings::NumberOfSubgrains != 0 && m_SubGrains.size() > 1) {
		for (vector<SubGrain*>::iterator itG = ++m_SubGrains.begin(); itG != m_SubGrains.end(); itG++) {
			if ( (*itG) != NULL ) {
				myQuaternion* newori = m_orientation->DisorientedFromReference(deviation, r);
				(*itG)->set_Orientation(*newori); //##MK::passes pointer to allocated heap object, requires only to copy properties via the operator=
			}
		}
	}
}

void Grains::generateSubGrainSEE(randomClass& r) {
	vector<SubGrain*>::iterator itG;
	double NrGrains = 0;
	double see_mu = m_SEE; //internal scatter about mean for the grain not m_PrefOri->SEE; because then there would be significantly less variety
	double see_sig_sgr = m_PrefOri->SEESubgrainScatter;
	if (Settings::NumberOfSubgrains != 0 && m_SubGrains.size() > 1)
		for (itG = ++m_SubGrains.begin(); itG != m_SubGrains.end(); itG++) {
			//MK::SEE of sub-grain becomes normally-distributed with width scatter about mean of the grain
			double see = r.r4_nor( see_mu, see_sig_sgr );
			(*itG)->set_SEE(see);
		}
}


void Grains::SubGrainConstructor(randomClass& r)
{
	if (Settings::NumberOfSubgrains == 0) {
		return;
	}

	bool is_periodic = false;
	if ( Settings::NumberOfGrains == 1 && Settings::VoronoiPeriodic == true ) {
		is_periodic = true;
	}

	int i = 1;
	double f = m_Volume / double(1. / Settings::NumberOfGrains) * pow( (1. / this->get_RelSizeScalingFromPrefori() ), 3.);
	int subgrains = Settings::NumberOfSubgrains * f;
	if (subgrains < 2) {
		cout << "grain " << m_ID << " gets " << 0 << " subgrains" << endl;
		return;
	}
	int blocks = (int) (pow((subgrains / 8), (1 / 3.)) + 1);
	if (blocks < 1)
		blocks = 1;
	voro::container con(getMinY(), getMaxY(), getMinX(), getMaxX(), getMinZ(),
			getMaxZ(), blocks, blocks, blocks, is_periodic, is_periodic, is_periodic, 8);

	// 5,5,5 = the number of grid blocks in each of the three coordinate directions
	// 2 = the initial memory allocation for each block

	/**********************************************************/
	// Randomly add particles into the container
	while (i <= subgrains) {
		double x = r.MersenneTwister() * (getMaxX() - getMinX()) + getMinX();
		double y = r.MersenneTwister() * (getMaxY() - getMinY()) + getMinY();
		double z = r.MersenneTwister() * (getMaxZ() - getMinZ()) + getMinZ();

		if ( int(x) == getMaxX() || int(y) == getMaxY() || int(z) == getMaxZ() ) {
			#pragma omp critical
			{
				cerr << "grain " << m_ID << " hitting bounds!" << "\n";
			}
		}

		if (Settings::Dimensionality == E_2D) {
			z = getMinZ();
		}
		if (m_localContainer->getValueAt(int(y), int(x), int(z)) == m_ID) {
			con.put(i, y, x, z); //Jedem Korn i werden die Koordinaten x,y,z zugeordnet = Koordinaten des Partikels (1 Korn = 1 Partikel)
			i++;
		}
	}

	#pragma omp critical
	{
		cout << "grain " << m_ID << " gets " << --i << " subgrains" << endl;
	}

	voro::voronoicell_neighbor c;
	vector<double> cellCoordinates;
	for (int k = getMinZ(); k < getMaxZ(); k++) {
		for (int j = getMinY(); j < getMaxY(); j++) {
			for (int i = getMinX(); i < getMaxX(); i++) {
				int cell_ID = 0;
				double rx, ry, rz;
				if (m_localContainer->getValueAt(j, i, k) == m_ID) {
					if (con.find_voronoi_cell(j, i, k, rx, ry, rz, cell_ID)) {
						m_localContainer->setValueAt(j, i, k, (unsigned int) cell_ID);
					}
				}
				else {
					m_localContainer->setValueAt(j, i, k, 0);
				}
			}
		}
	}
	c_loop_all vl(con);
	m_SubGrains.resize(subgrains + 1);
	if (vl.start()) {
		do {
			double rx, ry, rz;
			con.compute_cell(c, vl);
			//new: get the grain_id
			int box_id = vl.pid();
			vl.pos(ry, rx, rz);
			c.vertices(ry, rx, rz, cellCoordinates);
			SubGrain* newSubgrain = new SubGrain(box_id, cellCoordinates,
					c.volume(), this);
			m_SubGrains[box_id] = newSubgrain;
		} while (vl.inc());
	}
	else {
		throw runtime_error("Voronoi container error at start() method!");
	}
	m_SubGrains[0] = NULL;
}
