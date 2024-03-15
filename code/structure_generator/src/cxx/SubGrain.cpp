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

#include "SubGrain.h"
#include "Grains.h"
#include "microStructureHdl.h"
#include "RTree.h"
#include "../../../thirdparty/mandatory/voroxx/voro/src/voro++.hh"
#include "utilities.h"
#include "Settings.h"
#include "Eigen/Dense"

using namespace std;
using namespace voro;
using namespace Eigen;

SubGrain::SubGrain() {
	m_orientation = NULL;
	m_owner = NULL;
}

SubGrain::SubGrain(myQuaternion ori, double SEE, double volume) :
		m_SEE(SEE), m_Volume(volume) {
	m_orientation = NULL;
	set_Orientation(ori);

	//##MK::m_owner = NULL; ?
}

SubGrain::SubGrain(int _id, vector<double> cellcoordinates, double volume,
		Grains* owner) :
		m_ID(_id), m_Volume(volume), m_owner(owner) {
	m_orientation = NULL;
	double xmin = Settings::NumberOfGridpoints, xmax = 0., ymin =
			Settings::NumberOfGridpoints, ymax = 0., zmin =
			Settings::NumberOfGridpoints, zmax = 0.;
	for (unsigned int i = 0; i < cellcoordinates.size() / 3; i++) {
		double candidate = cellcoordinates.at(3 * i + 1);
		if (candidate < xmin)
			xmin = candidate;
		if (candidate > xmax)
			xmax = candidate;
		candidate = cellcoordinates.at(3 * i);
		if (candidate < ymin)
			ymin = candidate;
		if (candidate > ymax)
			ymax = candidate;
		candidate = cellcoordinates.at(3 * i + 2);
		if (candidate < zmin)
			zmin = candidate;
		if (candidate > zmax)
			zmax = candidate;
//		cout << cellcoordinates.at(3 * i + 1) << "  "
//				<< cellcoordinates.at(3 * i) << "  "
//				<< cellcoordinates.at(3 * i + 2) << " \n";
	}
	m_xmin = int(xmin);
	m_ymin = int(ymin);
	m_zmin = int(zmin);

	m_xmax = int(xmax);
	m_ymax = int(ymax);
	m_zmax = int(zmax);

	if (m_xmax < (Settings::NumberOfGridpoints - 2.))
		m_xmax++;

	if (m_ymax < (Settings::NumberOfGridpoints - 2.))
		m_ymax++;

	if (m_zmax < (Settings::NumberOfGridpoints - 2.))
		m_zmax++;
	m_SEE = 0.;

	int x_min = Settings::NumberOfGridpoints, x_max = 0, y_min =
			Settings::NumberOfGridpoints, y_max = 0, z_min =
			Settings::NumberOfGridpoints, z_max = 0;
	m_Volume =0;
	for (int k = m_zmin; k < m_zmax; k++)
		for (int j = m_ymin; j < m_ymax; j++)
			for (int i = m_xmin; i < m_xmax; i++) {
				if (m_owner->m_localContainer->isPointInside(j, i, k))
					if (m_ID
							== m_owner->m_localContainer->getValueAt(j, i, k)) {
						if (i < x_min)
							x_min = i;
						if (i > x_max)
							x_max = i;
						if (j < y_min)
							y_min = j;
						if (j > y_max)
							y_max = j;
						if (k < z_min)
							z_min = k;
						if (k > z_max)
							z_max = k;
						m_Volume = m_Volume+1;
					}
			}
	m_xmin = x_min;
	m_ymin = y_min;
	m_zmin = z_min;
	m_xmax = x_max;
	m_ymax = y_max;
	m_zmax = z_max;
}


void SubGrain::setID(int id) {
	m_ID = id;
}


SubGrain::~SubGrain() {
	if ( m_orientation != NULL )
		delete m_orientation;
	//MK::do not delete owner only raw-pointer backreference!
}

