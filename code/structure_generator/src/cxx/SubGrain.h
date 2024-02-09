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

#ifndef SUBGRAIN_H_
#define SUBGRAIN_H_

#include <string>
#include "dimensionalBuffer.h"
#include "myQuaternion.h"
#include "RTree.h"
class container;
class Grains;

using namespace std;

class SubGrain {
public:
	//Constructors
	SubGrain();
	SubGrain(myQuaternion ori, double SEE, double volume);
	SubGrain(int id, vector<double> cellcoordinates, double volume, Grains* owner);

	//Destructor
	virtual ~SubGrain();

	//Set Functions
	void setID(int id);
	inline void set_Orientation(myQuaternion ori) {
		if (m_orientation == NULL)
			m_orientation = new myQuaternion(); //##MK::potential memory leak, because specificallyDisorient has already allocated heap object
		*m_orientation = ori;
	}
	inline myQuaternion* get_Orientation() {
		return m_orientation;
	}
	inline double get_Volume() {
		return m_Volume;
	}
	inline int get_ID() {
		return m_ID;
	}
	inline double get_SEE() {
		return m_SEE;
	}
	inline void set_SEE(double see) {
		m_SEE = see;
	}
	inline int getMinX() const {
		return m_xmin;
	}
	inline int getMaxX() const {
		return m_xmax;
	}
	inline int getMinY() const {
		return m_ymin;
	}
	inline int getMaxY() const {
		return m_ymax;
	}
	inline int getMinZ() const {
		return m_zmin;
	}
	inline int getMaxZ() const {
		return m_zmax;
	}


private:
	myQuaternion* m_orientation;
	double m_SEE;
	double m_Volume;
	int m_ID;
	int m_xmin;
	int m_xmax;
	int m_ymin;
	int m_ymax;
	int m_zmin;
	int m_zmax;
	Grains* m_owner; //MK::do not delete only backreference

};

#endif /* SUBGRAIN_H_ */
