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
#include "LsBox.h"
#include "Settings.h"
#include "DimensionalBufferReal.h"
#include "pooledDimensionalBufferReal.h"
#include "ContourSector.h"
#include "GrahamScan.h"
#include "Utilities.h"
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include"Structs.h"

#define PERIODIC(x, f) (((x)+f)%f)


LSbox::LSbox(int id, double phi1, double PHI, double phi2, grainhdl* owner) :
		m_ID(id), m_exists(true), m_grainHandler(owner), m_grainBoundary(this), m_isMotionRegular(
				true), m_intersectsBoundaryGrain(false), m_volume(0), m_energy(
				0), m_perimeter(0) {
	m_orientationQuat = new Quaternion();
	double euler[3] = { phi1, PHI, phi2 };
	m_orientationQuat->euler2quaternion(euler);
	m_StoredElasticEnergy = 0.;
	m_StoredElasticEnergy *= Settings::DislocEnPerM / Settings::HAGB_Energy * Settings::Physical_Domain_Size; // normierung

	m_inputDistance = new DimensionalBufferReal(0, 0, 0, 0);
	m_outputDistance = new DimensionalBufferReal(0, 0, 0, 0);
	if (Settings::UseMagneticField)
		calculateMagneticEnergy();
}

/*
LSbox::LSbox(int aID, vector<SPoint>& contour, grainhdl* owner) :
		m_ID(aID), m_exists(true), m_grainHandler(owner), m_grainBoundary(this), m_isMotionRegular(
				true), m_intersectsBoundaryGrain(false), m_volume(0), m_energy(
				0), m_perimeter(0) {
	m_StoredElasticEnergy = 0.0;
	int grid_blowup = owner->get_grid_blowup();
	if (Settings::UseMagneticField)
		calculateMagneticEnergy();
	double h = owner->get_h();
	// determine size of grain
	m_orientationQuat = new Quaternion();
#pragma omp critical
	{
		if (Settings::UseTexture) {
			double newOri[3];
			(*(m_grainHandler->mymath)).newOrientationFromReference(
					m_grainHandler->bunge, m_grainHandler->deviation, newOri);
			m_orientationQuat->euler2quaternion(newOri);
		} else
			m_orientationQuat->randomOriShoemakeQuat(m_grainHandler->mymath);
	}
	int xmax = 0;
	int xmin = m_grainHandler->get_ngridpoints();
	int ymax = 0;
	int ymin = xmin;

	m_grainBoundary.getRawBoundary() = contour;

	double x, y;
	for (unsigned int k = 0; k < m_grainBoundary.getRawBoundary().size(); k++) {
		y = m_grainBoundary.getRawBoundary()[k].y;
		x = m_grainBoundary.getRawBoundary()[k].x;
		if (y / h < ymin)
			ymin = y / h;
		if (y / h > ymax)
			ymax = y / h + 1;
		if (x / h < xmin)
			xmin = x / h;
		if (x / h > xmax)
			xmax = x / h + 1;
	}
	xmax += 2 * grid_blowup;
	ymax += 2 * grid_blowup;

	m_inputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);
	m_outputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);

	m_inputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());
	m_outputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());
	//	inputDistance->clearValues(0.0);
	//	outputDistance->clearValues(0.0);

	reizeIDLocalToDistanceBuffer();

	// 	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
}


//constructor for read from file
LSbox::LSbox(int id, const vector<SPoint>& vertices, double q1, double q2,
		double q3, double q4, grainhdl* owner) :
		m_ID(id), m_exists(true), m_grainHandler(owner), m_grainBoundary(this), m_isMotionRegular(
				true), m_intersectsBoundaryGrain(false), m_volume(0), m_energy(
				0), m_perimeter(0) {
	m_orientationQuat = new Quaternion(q1, q2, q3, q4);
	m_grainBoundary.getRawBoundary() = vertices;
	if (Settings::UseMagneticField)
		calculateMagneticEnergy();
	m_StoredElasticEnergy = 0.0;
	int grid_blowup = m_grainHandler->get_grid_blowup();
	double h = m_grainHandler->get_h();
	// determine size of grain
	int xmax = 0;
	int xmin = m_grainHandler->get_ngridpoints();
	int ymax = 0;
	int ymin = xmin;

	double y, x;
	for (int k = 0; k < vertices.size(); k++) {
		y = vertices[k].y;
		x = vertices[k].x;
		if (y / h < ymin)
			ymin = y / h;
		if (y / h > ymax)
			ymax = y / h;
		if (x / h < xmin)
			xmin = x / h;
		if (x / h > xmax)
			xmax = x / h;
	}
	xmax += 2 * grid_blowup;
	ymax += 2 * grid_blowup;

	if (ymax > m_grainHandler->get_ngridpoints())
		ymax = m_grainHandler->get_ngridpoints();
	if (xmax > m_grainHandler->get_ngridpoints())
		xmax = m_grainHandler->get_ngridpoints();
	//	cout << "constructed a box with size: "<< xmin << "  " << xmax << "  " << ymin << "  " << xmax << "  " << endl;
	m_inputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);
	m_outputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);

	m_inputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());
	m_outputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());
	//	inputDistance->clearValues(0.0);
	//	outputDistance->clearValues(0.0);

	reizeIDLocalToDistanceBuffer();

	// 	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
}

// Constructor for Voxelized InputData
LSbox::LSbox(int id, const vector<SPoint>& vertices, Quaternion ori,
		double StoredElasticEnergy, grainhdl* owner) :
		m_ID(id), m_exists(true), m_grainHandler(owner), m_grainBoundary(this), m_isMotionRegular(
				true), m_intersectsBoundaryGrain(false), m_volume(0), m_energy(
				0), m_perimeter(0), m_StoredElasticEnergy(StoredElasticEnergy) {
	m_orientationQuat = new Quaternion(ori.get_q0(), ori.get_q1(), ori.get_q2(),
			ori.get_q3());
	//	m_grainBoundary.getRawBoundary() = vertices;

	if (Settings::UseMagneticField)
		calculateMagneticEnergy();
	if (Settings::UseStoredElasticEnergy == true)
		m_StoredElasticEnergy *= Settings::DislocEnPerM / Settings::HAGB_Energy
				* Settings::Physical_Domain_Size; // normierung
	int grid_blowup = m_grainHandler->get_grid_blowup();
	// determine size of grain
	int xmax = 0;
	int xmin = m_grainHandler->get_ngridpoints();
	int ymax = 0;
	int ymin = xmin;

	double y, x;
	for (int k = 0; k < vertices.size(); k++) {
		y = vertices[k].y;
		x = vertices[k].x;
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
	}
	xmax += 2 * grid_blowup;
	ymax += 2 * grid_blowup;

	if (ymax > m_grainHandler->get_ngridpoints())
		ymax = m_grainHandler->get_ngridpoints();
	if (xmax > m_grainHandler->get_ngridpoints())
		xmax = m_grainHandler->get_ngridpoints();
	if (ymin < 0)
		ymin = 0;
	if (xmin < 0)
		xmin = 0;
	//	cout << "constructed a box with size: "<< xmin << "  " << xmax << "  " << ymin << "  " << xmax << "  " << endl;

	m_inputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);
	m_outputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);
	m_inputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());
	m_outputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());
	//	inputDistance->clearValues(0.0);
	//	outputDistance->clearValues(0.0);

	reizeIDLocalToDistanceBuffer();
	// 	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
}
*/


LSbox::LSbox(int id, const vector<SPoint>& vertices, vector<double> const & quaternion, 
		const double StoredElasticEnergy, grainhdl* owner) :
		m_ID(id), m_exists(true), m_grainHandler(owner), m_grainBoundary(this), m_isMotionRegular(
				true), m_intersectsBoundaryGrain(false), m_volume(0), m_energy(
				0), m_perimeter(0), m_StoredElasticEnergy(StoredElasticEnergy) {
	m_orientationQuat = new Quaternion(quaternion[0], quaternion[1], quaternion[2], quaternion[3]);
	if (Settings::UseMagneticField) {
		calculateMagneticEnergy();
	}
	if (Settings::UseStoredElasticEnergy == true) {
		//real time scaling
		m_StoredElasticEnergy *= Settings::DislocEnPerM / Settings::HAGB_Energy
				* Settings::Physical_Domain_Size;
	}
	int grid_blowup = m_grainHandler->get_grid_blowup();
	// determine size of grain
	int xmax = 0;
	int xmin = m_grainHandler->get_ngridpoints();
	int ymax = 0;
	int ymin = xmin;

	double y, x;
	for (int k = 0; k < vertices.size(); k++) {
		y = vertices[k].y;
		x = vertices[k].x;
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
	}
	xmax += 2 * grid_blowup;
	ymax += 2 * grid_blowup;

	if (ymax > m_grainHandler->get_ngridpoints())
		ymax = m_grainHandler->get_ngridpoints();
	if (xmax > m_grainHandler->get_ngridpoints())
		xmax = m_grainHandler->get_ngridpoints();
	if (ymin < 0)
		ymin = 0;
	if (xmin < 0)
		xmin = 0;
	//	cout << "constructed a box with size: "<< xmin << "  " << xmax << "  " << ymin << "  " << xmax << "  " << endl;

	m_inputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);
	m_outputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);
	m_inputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());
	m_outputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());
	//	inputDistance->clearValues(0.0);
	//	outputDistance->clearValues(0.0);

	reizeIDLocalToDistanceBuffer();
	// 	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;
}


/*
LSbox::LSbox(int id, int nedges, double* edges, double phi1, double PHI,
		double phi2, grainhdl* owner) :
		m_ID(id), m_exists(true), m_grainHandler(owner), m_grainBoundary(this), m_isMotionRegular(
				true), m_intersectsBoundaryGrain(false), m_volume(0), m_energy(
				0), m_perimeter(0) {
	if (id == 1) {
		if (Settings::UseMagneticField)
			calculateMagneticEnergy();
//		cout << "Volumenenergy Korn 1: " << m_magneticEnergy << endl;
	} else
		m_magneticEnergy = 0;
	m_orientationQuat = new Quaternion();
	double euler[3] = { phi1, PHI, phi2 };
	m_orientationQuat->euler2quaternion(euler);

	//! Add contour grain, the plus one takes account of the
	//! assumed datastructure in the following code segments:
	//! the first point is also the last one in the contourGrain-Array
	m_grainBoundary.getRawBoundary().resize(nedges + 1);

	int grid_blowup = owner->get_grid_blowup();
	double h = owner->get_h();
	// determine size of grain
	int xmax = 0;
	int xmin = m_grainHandler->get_ngridpoints();
	int ymax = 0;
	int ymin = xmin;

	double x1[2], x2[2];
	for (int k = 0; k < nedges; k++) {
		x1[0] = edges[(4 * k) + 1];
		x1[1] = edges[4 * k];
		x2[0] = edges[(4 * k) + 3];
		x2[1] = edges[(4 * k) + 2];

		//! Add contour grain points
		m_grainBoundary.getRawBoundary()[k].x = edges[4 * k];
		m_grainBoundary.getRawBoundary()[k].y = edges[(4 * k) + 1];
		//! Consider that the last point equals the first one
		if (k == nedges - 1) {
			m_grainBoundary.getRawBoundary()[k + 1].x = edges[(4 * k) + 2];
			m_grainBoundary.getRawBoundary()[k + 1].y = edges[(4 * k) + 3];
		}

		//	for convention:
		//	x[i][j]:
		//	i = Zeilenindex(y-direction)
		// 	j = Spaltenindex(x-direction)

		// check for "Zeilen" Minima/Maxima
		if (x1[0] / h < ymin)
			ymin = x1[0] / h;
		if (x2[0] / h < ymin)
			ymin = x2[0] / h;

		if (x1[0] / h > ymax)
			ymax = (x1[0] / h);
		if (x2[0] / h > ymax)
			ymax = (x2[0] / h);

		// check for "Spalten" Minima/Maxima
		if (x1[1] / h < xmin)
			xmin = x1[1] / h;
		if (x2[1] / h < xmin)
			xmin = x2[1] / h;

		if (x1[1] / h > xmax)
			xmax = (x1[1] / h);
		if (x2[1] / h > xmax)
			xmax = (x2[1] / h);
	}

	xmax += 2 * grid_blowup;
	ymax += 2 * grid_blowup;

	m_inputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);
	m_outputDistance = new DimensionalBufferReal(xmin, ymin, xmax, ymax);

	m_inputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());
	m_outputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());
	//	inputDistance->clearValues(0.0);
	//	outputDistance->clearValues(0.0);

	reizeIDLocalToDistanceBuffer();

	// 	cout << "made a new box: xmin="<<xmin<< " xmax="<<xmax <<" ymin="<<ymin << " ymax="<<ymax<<endl;

}
*/

LSbox::~LSbox() {
	if (m_orientationQuat != NULL)
		delete m_orientationQuat;
	delete m_inputDistance;
	delete m_outputDistance;
}

double LSbox::get_h() {
	return m_grainHandler->get_h();
}
void LSbox::calculateDistanceFunction() {

	SPoint to_test;

	if (Settings::MicrostructureGenMode == E_READ_VOXELIZED_MICROSTRUCTURE) {
		int min = m_grainHandler->get_grid_blowup();
		int max = m_grainHandler->get_ngridpoints() - min - 1;
		for (int i = m_inputDistance->getMinY(); i < m_inputDistance->getMaxY();
				i++) {
			for (int j = m_inputDistance->getMinX();
					j < m_inputDistance->getMaxX(); j++) {
				if (i < min || i > max || j < min || j > max)
					m_inputDistance->setValueAt(i, j, -m_grainHandler->get_h());
				else if (m_ID == m_grainHandler->IDField->getValueAt(i, j))
					m_inputDistance->setValueAt(i, j, m_grainHandler->get_h());
				else
					m_inputDistance->setValueAt(i, j, -m_grainHandler->get_h());
			}
		}
		executeRedistancing();
		//plot_box(true, 1, "init", true);
		//plot_box(true, 2, "init", true);
		marchingSquares(m_outputDistance);
		//plot_box(true, 1, "init", true);
	} else {
		int grid_blowup = m_grainHandler->get_grid_blowup();
		double h = m_grainHandler->get_h();
		int i = 0, j = 0;
		SPoint to_test;
		vector<SPoint> &contourGrain = m_grainBoundary.getRawBoundary();
		int contour_size = contourGrain.size();

		//! Added the visit of the last point in the outputDistance for both directions
		for (i = m_outputDistance->getMinY(); i < m_outputDistance->getMaxY();
				i++) {
			for (j = m_outputDistance->getMinX();
					j < m_outputDistance->getMaxX(); j++) {
				to_test.x = (j - grid_blowup) * h;
				to_test.y = (i - grid_blowup) * h;

				bool isInside = false;

				for (int k = 1, l = 0; k < contour_size; k++) {
					//! This PointInPolygon test proofed more valid
					//! for a larger amount of geometrical configurations
					if (((contourGrain[l].y > to_test.y)
							!= (contourGrain[k].y > to_test.y))
							&& (to_test.x
									< (contourGrain[k].x - contourGrain[l].x)
											* (to_test.y - contourGrain[l].y)
											/ (contourGrain[k].y
													- contourGrain[l].y)
											+ contourGrain[l].x)) {
						isInside = !isInside;
					}
					l = k;
				}

				double minDist = 1000000.0;
				for (int k = 1, l = 0; k < contour_size; k++) {
					SPoint u = contourGrain[k] - contourGrain[l];
					double lambda = (to_test - contourGrain[l]).dot(u);
					lambda /= u.dot(u);
					//! For a lamdba that equals 0 or 1 the point to point distance calculation is used
					double dist;
					if (lambda <= 0) {
						dist = (to_test - contourGrain[l]).len();
					} else if (lambda >= 1) {
						dist = (contourGrain[k] - to_test).len();
					} else {
						dist = (to_test - (contourGrain[l] + u * lambda)).len();
					}
					minDist = min(minDist, dist);
					l = k;
				}
				if (minDist > m_grainHandler->delta)
					minDist = m_grainHandler->delta;
				m_outputDistance->setValueAt(i, j,
						isInside ? minDist : -minDist);
			}
		}
		m_volume = computeVolume()
				/ (m_grainHandler->get_h() * m_grainHandler->get_h());
	}
}

// Convolution und Helperfunctions
/**************************************/
/**************************************/

void LSbox::executeConvolution(ExpandingVector<char>& mem_pool) {

	//plot_box(true, 1, "BeforeConvo", true);
	double h = m_grainHandler->get_h();
	if (grainExists() != true)
		return;
	//  set references for the convolution step
#ifdef USE_FFTW
	fftwp_complex *fftTemp = (fftwp_complex*) &mem_pool[0];
	convolutionGeneratorFFTW(fftTemp, m_forwardPlan, m_backwardsPlan);
#elif defined USE_MKL
	MKL_Complex16* fftTemp = (MKL_Complex16*) &mem_pool[0];
	convolutionGeneratorMKL(fftTemp);
#endif
	/*********************************************************************************/
	// Velocity Corrector Step:
	/*********************************************************************************/
	// hier soll energycorrection gerechnet werden.
	// in der domainCl steht die urspruenglich distanzfunktion, in dem arry die gefaltete
	//TEST CODE
	if (!Settings::IsIsotropicNetwork && m_grainHandler->loop != 0
			&& m_isMotionRegular == true) {
		constructBoundarySectors();
		vector<LSbox*>::iterator it;
		int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax;
		double weight;
		double val;
		vector<LSbox*> IDs;
		vector<LSbox*> IDsActive;

		//correction of the shrinkage of small grains with radius rCrit < rLimit
		double actualGrainRadius = sqrt(getVolume() / PI) / h;
		double rLimit = 1.0;
		if (Settings::ConvolutionMode == E_LAPLACE)
			rLimit = 20.0;
		else if (Settings::ConvolutionMode == E_LAPLACE_RITCHARDSON)
			rLimit = 25.0;
		else if (Settings::ConvolutionMode == E_GAUSSIAN)
			rLimit = 20.0;
		//! Linear fitting

		//! Quadratic fitting
		//! Square root fitting

		if (m_IDLocal.getMinX() < m_outputDistance->getMinX())
			intersec_xmin = m_outputDistance->getMinX();
		else
			intersec_xmin = m_IDLocal.getMinX();

		if (m_IDLocal.getMinY() < m_outputDistance->getMinY())
			intersec_ymin = m_outputDistance->getMinY();
		else
			intersec_ymin = m_IDLocal.getMinY();

		if (m_IDLocal.getMaxX() > m_outputDistance->getMaxX())
			intersec_xmax = m_outputDistance->getMaxX();
		else
			intersec_xmax = m_IDLocal.getMaxX();

		if (m_IDLocal.getMaxY() > m_outputDistance->getMaxY())
			intersec_ymax = m_outputDistance->getMaxY();
		else
			intersec_ymax = m_IDLocal.getMaxY();

		for (int i = intersec_ymin; i < intersec_ymax; i++) {
			for (int j = intersec_xmin; j < intersec_xmax; j++) {
				val = m_inputDistance->getValueAt(i, j);
				if (val > -m_grainHandler->delta) {
					double radiuscorrection = 1.0;
					if (m_grainHandler->convolutionCorrection) {
						if (actualGrainRadius < rLimit) {
							if (Settings::ConvolutionMode == E_LAPLACE) {
								const double yInterceptBottom = 0.83;
								const double cSlope = (1 - yInterceptBottom)
										/ cbrt(rLimit);
								radiuscorrection = (cSlope
										* cbrt(actualGrainRadius)
										+ yInterceptBottom);
							} else if (Settings::ConvolutionMode
									== E_LAPLACE_RITCHARDSON) {
								const double yInterceptBottom = 0.77;
								double const a = 4.7;
								radiuscorrection = 1
										- exp(-actualGrainRadius / a)
												* (1 - yInterceptBottom);
							} else if (Settings::ConvolutionMode
									== E_GAUSSIAN) {
								const double yInterceptBottom = 0.77;
								double const a = 4.5;
								radiuscorrection = 1.0 - exp(-actualGrainRadius / a) * (1.0 - yInterceptBottom);
							}
						}
					}
					weight = getWeight(i, j);
					m_outputDistance->setValueAt(i, j,
							val
									+ (m_outputDistance->getValueAt(i, j) - val)
											* weight * radiuscorrection);
				} else {
					m_outputDistance->setValueAt(i, j, -m_grainHandler->delta);
				}

				if (Settings::UseStoredElasticEnergy == true) {
					double f_StoredElasticEnergy =
							m_grainBoundary.get_f_StoredElasticEnergy(
									m_IDLocal.getValueAt(i, j).grainID);
					m_outputDistance->setValueAt(i, j,
						(m_outputDistance->getValueAt(i, j) - (f_StoredElasticEnergy * m_grainHandler->get_dt() / m_grainHandler->TimeSlope ))); //##MK
				}
				if (Settings::UseMagneticField) {
					double f_magneticEnergy =
							m_grainBoundary.get_f_magneticEnergy(
									m_IDLocal.getValueAt(i, j).grainID)
									* m_grainHandler->get_dt() / m_grainHandler->TimeSlope;  //##MK

					m_outputDistance->setValueAt(i, j,
							(m_outputDistance->getValueAt(i, j)
									- f_magneticEnergy));
				}

			}
		}
	}
	//plot_box(true, 1, "BeforeConvoGaussian", true);
	//plot_box(true, 2, "AfterConvoGaussian", true);
	reizeIDLocalToDistanceBuffer();
	m_IDLocal.clear();
}
#ifdef USE_FFTW
void LSbox::destroyFFTWs() {
	fftw_destroy_planp(m_forwardPlan);
	fftw_destroy_planp(m_backwardsPlan);
}
#endif

double LSbox::getGBEnergyTimesGBMobility(int i, int j) {
	//LSbox* neighbour = IDLocal.getValueAt(i, j).getElementAt(0);
	LSbox* neighbour = m_grainHandler->getGrainByID(
			m_IDLocal.getValueAt(i, j).grainID);
	characteristics& found = m_grainBoundary.getDirectNeighbourCaracteristic(
			neighbour);
	if (found.energyDensity * found.mobility > 0.01)
		return found.energyDensity * found.mobility;
	else
		return 0.01;
}

double LSbox::getGBEnergyTimesGBMobility(LSbox* neighbour) {
	characteristics& found = m_grainBoundary.getDirectNeighbourCaracteristic(
			neighbour);
	if (found.energyDensity * found.mobility > 0.01)
		return found.energyDensity * found.mobility;
	else
		return 0.01;
}

double LSbox::getGBEnergy(LSbox* neighbour) {
	characteristics& found = m_grainBoundary.getDirectNeighbourCaracteristic(
			neighbour);
	return found.energyDensity;
}

void LSbox::reizeIDLocalToDistanceBuffer() {
	int xmaxId = m_outputDistance->getMaxX();
	int xminId = m_outputDistance->getMinX();
	int ymaxId = m_outputDistance->getMaxY();
	int yminId = m_outputDistance->getMinY();

	m_IDLocal.resize(xminId, yminId, xmaxId, ymaxId);
}
#ifdef USE_FFTW
void LSbox::makeFFTPlans(double *in, double* out, fftw_complex *fftTemp,
		fftw_plan *fftplan1, fftw_plan *fftplan2) { /* creates plans for FFT and IFFT */
	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	*fftplan1 = fftw_plan_dft_r2c_2d(n, n, in, fftTemp, FFTW_ESTIMATE);
	*fftplan2 = fftw_plan_dft_c2r_2d(n, n, fftTemp, out, FFTW_ESTIMATE);
	/*
	 The flags argument is usually either FFTW_MEASURE or FFTW_ESTIMATE. FFTW_MEASURE
	 instructs FFTW to run and measure the execution time of several FFTs in order to find the
	 best way to compute the transform of size n. This process takes some time (usually a few
	 seconds), depending on your machine and on the size of the transform. FFTW_ESTIMATE,
	 on the contrary, does not run any computation and just builds a reasonable plan that is
	 probably sub-optimal. In short, if your program performs many transforms of the same size
	 and initialization time is not important, use FFTW_MEASURE; otherwise use the estimate. */
}

//how to compile both single and double precision FFTW libs?
//two-step process
void LSbox::makeFFTPlans(float *in, float* out, fftwf_complex *fftTemp,
		fftwf_plan *fftplan1, fftwf_plan *fftplan2) { /* creates plans for FFT and IFFT */
	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	*fftplan1 = fftwf_plan_dft_r2c_2d(n, n, in, fftTemp, FFTW_ESTIMATE);
	*fftplan2 = fftwf_plan_dft_c2r_2d(n, n, fftTemp, out, FFTW_ESTIMATE);

}

void LSbox::makeFFTPlans(ExpandingVector<char>& memory_dump) {
	fftwp_complex *fftTemp = (fftwp_complex*) &memory_dump[0];

	makeFFTPlans(m_inputDistance->getRawData(), m_outputDistance->getRawData(),
			fftTemp, &m_forwardPlan, &m_backwardsPlan);
}
#endif

void LSbox::preallocateMemory(ExpandingVector<char>& memory_dump) {
	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	int desired_size = 0;
#ifdef USE_FFTW
	desired_size = n * (floor(n / 2) + 1) * sizeof(fftwp_complex);
#elif defined USE_MKL
	desired_size = n * (floor(n / 2) + 1) * sizeof(MKL_Complex16);
#endif
	memory_dump.expand(desired_size);
}
#ifdef USE_FFTW
void LSbox::convolutionGeneratorFFTW(fftwp_complex *fftTemp, fftwp_plan fftplan1,
		fftwp_plan fftplan2) {
	/* Function returns in u the updated value of u as described below..
	 u -> (G_{dt})*u
	 Assumptions:
	 fftplan1 converts u to it's FT (in fftTemp), and
	 fftplan2 converts the FT (after pointwise multiplication with G)
	 (*outputDistance) to a real-valued level set function at u.
	 Memory is already allocated in fftTemp
	 (necessary to create the plans) */

	int n = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
	double dt = m_grainHandler->get_dt();
	int n2 = floor(n / 2) + 1;
	int nn = (*m_grainHandler).get_realDomainSize();
	double nsq = nn * nn;
	double k = 2.0 * PI / n;
	double G;
	double coski;
	int j2;
	int i2;

	executeFFTW(fftplan1);
	//	Forward DFT

	switch (Settings::ConvolutionMode) {
		case E_LAPLACE: {
			for (int i = 0; i < n2; i++) {
				coski = cos(k * i);
				for (int j = 0; j < n; j++) {
					G = 2.0 * (2.0 - coski - cos(k * j)) * nsq;
					G = 1.0 / (1.0 + (dt * G)) / (n * n);
					fftTemp[i + n2 * j][0] = fftTemp[i + n2 * j][0] * G;
					fftTemp[i + n2 * j][1] = fftTemp[i + n2 * j][1] * G;
				}
			}
			break;
		}
		case E_LAPLACE_RITCHARDSON: {
			//			Ritchardson Extrapolation
			for (int i = 0; i < n2; i++) {
				coski = cos(k * i);
				for (int j = 0; j < n; j++) {
					G = 2.0 * (2.0 - coski - cos(k * j)) * nsq;
					G = (4.0 / pow(1 + 1.5 * (dt) / 40 * G, 40) - 1.0 / pow(
									1 + 3.0 * (dt) / 40 * G, 40)) / 3.0 / (double) (n * n);
					fftTemp[i + n2 * j][0] = fftTemp[i + n2 * j][0] * G;
					fftTemp[i + n2 * j][1] = fftTemp[i + n2 * j][1] * G;
				}
			}
			break;
		}
		case E_GAUSSIAN: {
			double n_nsq = n * n;

			//	Convolution with Normaldistribution
			for (int i = 0; i < n2; i++) {
				i2 = mymin(i,n-i);
				for (int j = 0; j < n; j++) {
					j2 = mymin(j,n-j);
					G = exp(-(i2 * i2 + j2 * j2) * 4.0 * dt * nsq / n_nsq * PI * PI) / n_nsq;
					fftTemp[i + n2 * j][0] = fftTemp[i + n2 * j][0] * G;
					fftTemp[i + n2 * j][1] = fftTemp[i + n2 * j][1] * G;
				}
			}
			break;
		}
		default:
		break;
	}

	executeFFTW(fftplan2);
	//	Inverse DFT
}
#elif defined USE_MKL

void LSbox::convolutionGeneratorMKL(MKL_Complex16* fftTemp)
{
	//	MKL_LONG dimensions[2];

	if (m_grainHandler->get_loop()==0) {
		m_dimensions[0]=0; m_dimensions[1]=0;
		m_b_input_strides[0] = 0; m_b_input_strides[2] = 1;
		m_b_output_strides[0] = 0; m_b_output_strides[2] = 1;
		m_f_input_strides[0] = 0; m_f_input_strides[2] = 1;
		m_f_output_strides[0] = 0; m_f_output_strides[2] = 1;
		m_handle =0; m_b_handle = 0;
	}

	//	m_handle = 0;

	DFTI_CONFIG_VALUE precision;

#if PRECISION > 0
	precision = DFTI_SINGLE;
#else
	precision = DFTI_DOUBLE;
#endif
	bool update_backward_plan = false;
	if (m_dimensions[0] != m_outputDistance->getMaxX() - m_outputDistance->getMinX()) {
		update_backward_plan = true;
		m_dimensions[0] = m_outputDistance->getMaxX() - m_outputDistance->getMinX();
		m_dimensions[1] = m_outputDistance->getMaxY() - m_outputDistance->getMinY();
		DftiFreeDescriptor(&m_handle);
		DftiCreateDescriptor(&m_handle, precision, DFTI_REAL, 2, m_dimensions);
		DftiSetValue(m_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		DftiSetValue(m_handle, DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX);
		m_f_input_strides[1] = m_dimensions[0];
		DftiSetValue(m_handle, DFTI_INPUT_STRIDES, m_f_input_strides);
		m_f_output_strides[1] = m_dimensions[0]/2 + 1;
		DftiSetValue(m_handle, DFTI_OUTPUT_STRIDES, m_f_output_strides);
		DftiCommitDescriptor(m_handle);
	}

	DftiComputeForward(m_handle, m_inputDistance->getRawData(), fftTemp);

	int n = m_dimensions[0];
	double dt = m_grainHandler->get_dt();
	int n2 = floor(n / 2) + 1;
	int nn = (*m_grainHandler).get_realDomainSize();
	double nsq = nn * nn;
	double k = 2.0 * PI / n;
	double G;
	double coski;
	int j2;
	int i2;

	switch (Settings::ConvolutionMode) {
		case E_LAPLACE: {
			for (int i = 0; i < n2; i++) {
				coski = cos(k * i);
				for (int j = 0; j < n; j++) {
					G = 2.0 * (2.0 - coski - cos(k * j)) * nsq;
					G = 1.0 / (1.0 + (dt * G)) / (n * n);
					fftTemp[i + n2 * j].real = fftTemp[i + n2 * j].real * G;
					fftTemp[i + n2 * j].imag = fftTemp[i + n2 * j].imag * G;
				}
			}
			break;
		}
		case E_LAPLACE_RITCHARDSON: {
			//			Ritchardson Extrapolation
			for (int i = 0; i < n2; i++) {
				coski = cos(k * i);
				for (int j = 0; j < n; j++) {
					G = 2.0 * (2.0 - coski - cos(k * j)) * nsq;
					G = (4.0 / pow(1 + 1.5 * (dt) / 40 * G, 40) - 1.0 / pow(
									1 + 3.0 * (dt) / 40 * G, 40)) / 3.0 / (double) (n * n);
					fftTemp[i + n2 * j].real = fftTemp[i + n2 * j].real * G;
					fftTemp[i + n2 * j].imag = fftTemp[i + n2 * j].imag * G;
				}
			}
			break;
		}
		case E_GAUSSIAN: {
			//cout << n * n << "   " << nsq << endl;
			double n_nsq = n * n;
			//			Convolution with Normaldistribution
			for (int i = 0; i < n2; i++) {
				i2 = mymin(i,n-i);
				for (int j = 0; j < n; j++) {
					j2 = mymin(j,n-j);
					G = exp(-(i2 * i2 + j2 * j2) * 4.0 * dt * nsq / n_nsq * PI * PI) / n_nsq;
					fftTemp[i + n2 * j].real = fftTemp[i + n2 * j].real * G;
					fftTemp[i + n2 * j].imag = fftTemp[i + n2 * j].imag * G;
				}
			}
			break;
		}
		default:
		break;
	}
	if (update_backward_plan) {
		DftiFreeDescriptor(&m_b_handle);
		DftiCreateDescriptor(&m_b_handle, precision, DFTI_REAL, 2, m_dimensions);
		DftiSetValue(m_b_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		DftiSetValue(m_b_handle, DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX);
		m_b_input_strides[1] = m_dimensions[0]/2 + 1;
		DftiSetValue(m_b_handle, DFTI_INPUT_STRIDES, m_b_input_strides);
		m_b_output_strides[1] = m_dimensions[0];
		DftiSetValue(m_b_handle, DFTI_OUTPUT_STRIDES, m_b_output_strides);
		//DFTI_BACKWARD_SCALE(m_b_handle, 1/n_nsq);
		DftiCommitDescriptor(m_b_handle);
	}
	DftiComputeBackward(m_b_handle, fftTemp, m_outputDistance->getRawData());

}
#endif

#ifdef USE_FFTW
void LSbox::executeFFTW(fftw_plan fftplan) {
	fftw_execute(fftplan);
}

void LSbox::executeFFTW(fftwf_plan fftplan) {
	fftwf_execute(fftplan);
}
#endif
/**************************************/
/**************************************/

/**************************************/
/**************************************/

void LSbox::switchInNOut() {
	DimensionalBufferReal* temp;

	temp = m_inputDistance;
	m_inputDistance = m_outputDistance;
	m_outputDistance = temp;
}

// Comparison + Helperfunctions
/**************************************/
/**************************************/

void LSbox::executeSetComparison() {
	m_newXMin = m_outputDistance->getMaxX();
	m_newXMax = m_outputDistance->getMinX();
	m_newYMin = m_outputDistance->getMaxY();
	m_newYMax = m_outputDistance->getMinY();
	for (int i = m_outputDistance->getMinY(); i < m_outputDistance->getMaxY();
			i++) {
		for (int j = m_outputDistance->getMinX();
				j < m_outputDistance->getMaxX(); j++) {
			if (abs(m_inputDistance->getValueAt(i, j))
					< 0.7 * m_grainHandler->delta) {
				m_outputDistance->setValueAt(i, j,
						0.5
								* (m_inputDistance->getValueAt(i, j)
										- m_outputDistance->getValueAt(i, j)));
			} else
				m_outputDistance->setValueAt(i, j,
						m_inputDistance->getValueAt(i, j));

			if (m_outputDistance->getValueAt(i, j) >= 0) {
				if (i < m_newYMin)
					m_newYMin = i;
				if (i > m_newYMax)
					m_newYMax = i;
				if (j < m_newXMin)
					m_newXMin = j;
				if (j > m_newXMax)
					m_newXMax = j;
			}
		}
	}
	m_newXMin -= m_grainHandler->get_grid_blowup();
	m_newXMax += m_grainHandler->get_grid_blowup();
	m_newYMin -= m_grainHandler->get_grid_blowup();
	m_newYMax += m_grainHandler->get_grid_blowup();
	if (m_newXMin >= m_newXMax || m_newYMin >= m_newYMax) {
		m_newXMin = m_outputDistance->getMinX();
		m_newXMax = m_outputDistance->getMaxX();
		m_newYMin = m_outputDistance->getMinY();
		m_newYMax = m_outputDistance->getMaxY();
	}
}

bool LSbox::checkIntersection(LSbox* box2) {
	if (m_inputDistance->getMinX() > box2->m_inputDistance->getMaxX()
			|| m_inputDistance->getMaxX() < box2->m_inputDistance->getMinX()
			|| m_inputDistance->getMinY() > box2->m_inputDistance->getMaxY()
			|| m_inputDistance->getMaxY() < box2->m_inputDistance->getMinY())
		return false;
	return true;
}

void LSbox::executeComparison() {
	if (grainExists() != true)
		return;

	m_outputDistance->clearValues(-1.0);

	std::vector<LSbox*>::iterator it_nn;

	m_secondOrderNeighbours = m_comparisonList;
	for (it_nn = m_secondOrderNeighbours.begin();
			it_nn != m_secondOrderNeighbours.end(); it_nn++) {
		int x_min_new, x_max_new, y_min_new, y_max_new;

		if (m_inputDistance->getMinX() < (**it_nn).m_inputDistance->getMinX())
			x_min_new = (**it_nn).m_inputDistance->getMinX();
		else
			x_min_new = m_inputDistance->getMinX();

		if (m_inputDistance->getMaxX() > (**it_nn).m_inputDistance->getMaxX())
			x_max_new = (**it_nn).m_inputDistance->getMaxX();
		else
			x_max_new = m_inputDistance->getMaxX();

		if (m_inputDistance->getMinY() < (**it_nn).m_inputDistance->getMinY())
			y_min_new = (**it_nn).m_inputDistance->getMinY();
		else
			y_min_new = m_inputDistance->getMinY();

		if (m_inputDistance->getMaxY() > (**it_nn).m_inputDistance->getMaxY())
			y_max_new = (**it_nn).m_inputDistance->getMaxY();
		else
			y_max_new = m_inputDistance->getMaxY();

		for (int i = y_min_new; i < y_max_new; i++) {
			for (int j = x_min_new; j < x_max_new; j++) {
				//				if (abs(inputDistance->getValueAt(i, j)) < m_grainHandler->delta) {
				double dist = (**it_nn).getDistanceFromInputBuff(i, j);
				//				if (abs(dist) < m_grainHandler->delta) {
				if (dist > m_outputDistance->getValueAt(i, j)) {
					m_outputDistance->setValueAt(i, j, dist);
					//IDLocal.getValueAt(i, j).insertAtPosition(E_FIRST_POSITION,*it_nn);
					m_IDLocal.getValueAt(i, j).grainID = (**it_nn).getID();
				}
			}
		}
	}
	if (BoundaryIntersection()) {
		m_intersectsBoundaryGrain = true;
		boundaryCondition();
	} else
		m_intersectsBoundaryGrain = false;
	//plot_box(true,2,"Com",true);
	//plot_box(true,3,"IDLocalContour", true);
}

bool LSbox::BoundaryIntersection() {
	int xMinBoundary = m_grainHandler->get_grid_blowup()
			+ m_grainHandler->getBoundaryGrainTube();
	int yMinBoundary = xMinBoundary;

	int xMaxBoundary = m_grainHandler->get_ngridpoints()
			- m_grainHandler->get_grid_blowup()
			- m_grainHandler->getBoundaryGrainTube();
	int yMaxBoundary = xMaxBoundary;

	if (m_outputDistance->getMinX() > xMinBoundary
			&& m_outputDistance->getMaxX() < xMaxBoundary
			&& m_outputDistance->getMinY() > yMinBoundary
			&& m_outputDistance->getMaxY() < yMaxBoundary)
		return false;
	else
		return true;
}

double LSbox::getDistanceFromInputBuff(int i, int j) {
	return m_inputDistance->getValueAt(i, j);
}

void LSbox::boundaryCondition() {
	int grid_blowup = m_grainHandler->get_grid_blowup();
	double h = m_grainHandler->get_h();
	int m = m_grainHandler->get_ngridpoints();
	int distXMin, distXMax, distX;
	int distYMin, distYMax, distY;
	double dist = 0;

	for (int i = m_inputDistance->getMinY(); i < m_inputDistance->getMaxY();
			i++) {
		for (int j = m_inputDistance->getMinX(); j < m_inputDistance->getMaxX();
				j++) {
			distXMin = -(j - grid_blowup);
			distYMin = -(i - grid_blowup);
			distXMax = (j - (m - grid_blowup));
			distYMax = (i - (m - grid_blowup));

			if (abs(distXMin) < abs(distXMax))
				distX = distXMin;
			else
				distX = distXMax;
			if (abs(distYMin) < abs(distYMax))
				distY = distYMin;
			else
				distY = distYMax;

			if (distX > 0 && distY > 0)
				dist = sqrt((double) distX * distX + distY * distY);

			else if (distX < 0 && distY > 0)
				dist = distY;
			else if (distX > 0 && distY < 0)
				dist = distX;
			else if (distX < 0 && distY < 0)
				dist = max(distX, distY);
			else if (distX == 0) {
				if (distY == 0)
					dist = 0;
				else if (distY < 0)
					dist = 0;
				else if (distY > 0)
					dist = distY;
			} else if (distY == 0) {
				if (distX < 0)
					dist = 0;
				else if (distX > 0)
					dist = distX;
			}

			if (dist * h > m_outputDistance->getValueAt(i, j)) {
				m_outputDistance->setValueAt(i, j, dist * h);
				//IDLocal.getValueAt(i, j).insertAtPosition(E_FIRST_POSITION, boundary);
				m_IDLocal.getValueAt(i, j).grainID = 0;
			}
			//			}
		}
	}
}

void LSbox::computeSecondOrderNeighbours() {
	vector<LSbox*> neighbourCandidates;
	if (grainExists() != true)
		return;
	//	vector<characteristics>::iterator it, it_ngC;
	//	vector<LSbox*>::iterator it_com, it_nC;
	bool just_in;
	//	neighbors_2order.clear();
	m_comparisonList.clear();
	// neighbors_2order gets a copy of ToCompare in comparison_box, so that it can be used here for reference for other objects.

	for (auto it = m_secondOrderNeighbours.begin();
			it != m_secondOrderNeighbours.end(); it++) {
		if ((*it)->grainExists() == true)
			m_comparisonList.push_back((*it));
		for (auto it_ngC = (*it)->m_secondOrderNeighbours.begin();
				it_ngC != (*it)->m_secondOrderNeighbours.end(); it_ngC++) {
			if ((*it_ngC)->grainExists() == true)
				if (checkIntersection((*it_ngC))) {
					neighbourCandidates.push_back((*it_ngC));
				}
		}
	}
	for (auto it_nC = neighbourCandidates.begin();
			it_nC != neighbourCandidates.end(); it_nC++) {
		just_in = false;
		if ((*it_nC) == this)
			continue;
		if ((*it_nC) == m_grainHandler->boundary)
			continue;
		for (auto it_com = m_comparisonList.begin();
				it_com != m_comparisonList.end(); it_com++) {
			if ((*it_com) == (*it_nC)) {
				just_in = true;
				break;
			}
		}
		if ((!just_in))
			m_comparisonList.push_back((*it_nC));
	}
	neighbourCandidates.clear();
}

void LSbox::setIDLocal(int ID) {
	IDChunkMinimal SetThisID;
	SetThisID.grainID = ID;
	m_IDLocal.clearValues(SetThisID);
}

/**************************************/
// end of Comparison
/**************************************/

// Find Contour operates on inputDistance
/**************************************/
/**************************************/

void LSbox::marchingSquares(DimensionalBufferReal* which)
{
	m_exists = m_grainBoundary.extractBoundaryAndJunctions(*which, m_IDLocal);
	m_volume = computeVolume();
}


void LSbox::extractContour()
{
	m_isMotionRegular = true;
	double newVolume = -1;
	int old_contour_len = m_grainBoundary.getBoundarySegmentCount();
	m_exists = m_grainBoundary.extractBoundaryAndJunctions(*m_inputDistance, m_IDLocal);
	if (!grainExists()) {
		return;
	}
	newVolume = computeVolume();
	if (isMotionRegular(old_contour_len, m_grainBoundary.getBoundarySegmentCount(), getVolume(), newVolume) == false) {
		do {
			clearContourGrainArea();
			m_exists = m_grainBoundary.extractBoundaryAndJunctions(*m_inputDistance, m_IDLocal, m_grainHandler->loop);
			if (!grainExists()) {
				cout << "we cleared a small grain " << m_ID << endl;
				return;
			}
			newVolume = computeVolume();
		} while (isMotionRegular(old_contour_len, m_grainBoundary.getBoundarySegmentCount(), getVolume(), newVolume) == false);
	}

	m_centroid = m_grainBoundary.calculateCentroid();
	int m = m_grainHandler->get_ngridpoints();
	bool out = false;
	if (m_newXMin < 0) {
		m_newXMin = 0;
		out = true;
	}
	if (m_newYMin < 0) {
		m_newYMin = 0;
		out = true;
	}
	if (m_newXMax > m) {
		m_newXMax = m;
		out = true;
	}
	if (m_newYMax > m) {
		m_newYMax = m;
		out = true;
	}
	if (out == true) {
		int loop = m_grainHandler->loop;
		cout << endl << "Timestep: " << m_grainHandler->loop << endl << endl;
		cout << "WARNING - undefined Boxsize in Box: " << m_ID
				<< " in Timestep: " << loop << "!!" << endl;
		cout << "Number of gridpoints: " << m << endl;
		cout << m_newYMin << " || " << m_newXMin << " || " << m_newYMax
				<< " || " << m_newXMax << endl;
	}
	m_outputDistance->resize(m_newXMin, m_newYMin, m_newXMax, m_newYMax);
	m_outputDistance->resizeToSquare(m_grainHandler->get_ngridpoints());

	m_perimeter = m_grainBoundary.computePerimeter();
	return;
}


bool LSbox::isMotionRegular(int old_contour, int new_contour, double old_volume, double new_volume)
{
	//Formula approximating von-Neumann-Mullins rule.
	//return (new_volume - old_volume)/m_grainHandler->get_dt() * (3/PI) >= -20.0;
	//	Formula approximating von-Neumann-Mullins rule.
	//return (new_volume - old_volume)/m_grainHandler->get_dt() * (3/PI) >= -20.0;
	double mullinsCriterion = (new_volume - old_volume)
			/ m_grainHandler->get_dt() * (3 / PI);
	double contourRatio = ((double) new_contour) / old_contour;
	double volumeRatio = new_volume / old_volume;

	if (mullinsCriterion > -6.0)
		// we are quiet save to detect a regular motion:
		return true;

	if (volumeRatio < 0.10) {
		//grain shrinks much faster than expected
		return false;
	}
	if (contourRatio < 0.35) {
		// we are quiet save to detect a spike if:
		return false;
	}
	if (contourRatio > 0.7) {
		return true;
	}

	// this is the grey zone: Mullins is < -6; 0.35 < contourRatio < 0.7
	// we are not save in neglecting,
	// so we do a flag and let the grain pass:
	m_isMotionRegular = false;
//	cout << "flagged ID " << m_ID << endl;
	return true;

	//	return (((double)new_contour)/old_contour >= 0.3) /*&& ((new_volume - old_volume)/m_grainHandler->get_dt() * (3/PI) >= -16.0) */&& ((double)new_volume / old_volume >= 0.3);
}

void LSbox::updateFirstOrderNeigbors() {
	if (grainExists() != true)
		return;
	m_grainBoundary.buildDirectNeighbours(*m_inputDistance, m_IDLocal);
}
double LSbox::computeVolume() {
	return m_grainBoundary.computeVolume();
}

void LSbox::computeVolumeAndEnergy() {
	if (grainExists() != true || m_isMotionRegular == false)
		return;

	m_energy = 0;
	vector<characteristics>::iterator it;

	double newVolume = abs(computeVolume()); //getVolume()
	m_grainBoundary.buildDirectNeighbours(*m_inputDistance, m_IDLocal);
	m_energy = m_grainBoundary.computeEnergy();
	//!
	//! Evaluating the area variation in the current time step and
	//! saving this variation together with the current number
	//! of neighbours in a vector for further analyses.
	//! The area variation is normalized by a factor coming from the
	//! Neumann-Mullins equation.
	//!
	m_volume = newVolume;
}

/**************************************/
//  Redistancing
/**************************************/

void LSbox::executeExactRedist() {
	if (grainExists() != true)
		return;
	double h = m_grainHandler->get_h();
	if (m_grainBoundary.getBoundarySegmentCount() == 0)
		marchingSquares(m_outputDistance);
	// 	resize the outputDistance array. be careful because during this part of algorithm both arrays have not the same size!!
	double delta = m_grainHandler->get_grid_blowup();
	for (int j = m_outputDistance->getMinX(); j < m_outputDistance->getMaxX();
			j++) {
		for (int i = m_outputDistance->getMinY();
				i < m_outputDistance->getMaxY() - 1; i++) {
			if (abs(m_outputDistance->getValueAt(i, j)) > delta * h) {
				m_outputDistance->setValueAt(i, j,
						sgn(m_outputDistance->getValueAt(i, j)) * delta * h);
				continue;
			}
			SPoint x(i, j, 0, 0);
			double Dist = m_grainBoundary.DistanceToGrainBondary(x);
			Dist *= sgn(m_outputDistance->getValueAt(i, j)) * h;
			m_outputDistance->setValueAt(i, j, Dist);
		}
	}
}

void LSbox::executeRedistancing() {
	if (grainExists() != true)
		return;
	double h = m_grainHandler->get_h();
	double candidate, i_slope, distToZero;

	m_outputDistance->clearValues(-1.0);

	// 	resize the outputDistance array. be careful because during this part of algorithm both arrays have not the same size!!
	int intersec_xmin, intersec_xmax, intersec_ymin, intersec_ymax;

	if (m_inputDistance->getMinX() < m_outputDistance->getMinX())
		intersec_xmin = m_outputDistance->getMinX();
	else
		intersec_xmin = m_inputDistance->getMinX();

	if (m_inputDistance->getMinY() < m_outputDistance->getMinY())
		intersec_ymin = m_outputDistance->getMinY();
	else
		intersec_ymin = m_inputDistance->getMinY();

	if (m_inputDistance->getMaxX() < m_outputDistance->getMaxX())
		intersec_xmax = m_inputDistance->getMaxX();
	else
		intersec_xmax = m_outputDistance->getMaxX();

	if (m_inputDistance->getMaxY() < m_outputDistance->getMaxY())
		intersec_ymax = m_inputDistance->getMaxY();
	else
		intersec_ymax = m_outputDistance->getMaxY();

	// y-direction forward
	for (int j = intersec_xmin; j < m_outputDistance->getMaxX(); j++) {
		for (int i = intersec_ymin; i < m_outputDistance->getMaxY() - 1; i++) {
			if (j < intersec_xmax && i < intersec_ymax - 1) {
				if (m_inputDistance->getValueAt(i, j)
						* m_inputDistance->getValueAt(i + 1, j) <= 0.0) {
					// interpolate
					i_slope = (m_inputDistance->getValueAt(i + 1, j)
							- m_inputDistance->getValueAt(i, j)) / h;
					distToZero = -m_inputDistance->getValueAt(i, j) / i_slope;
					if (abs(m_outputDistance->getValueAt(i, j))
							> abs(distToZero))
						m_outputDistance->setValueAt(i, j,
								-distToZero * sgn(i_slope));
				}
				// calculate new distance candidate and assign if appropriate
				candidate = m_outputDistance->getValueAt(i, j)
						+ (sgn(m_inputDistance->getValueAt(i + 1, j)) * h);
				if (abs(candidate)
						< abs(m_outputDistance->getValueAt(i + 1, j)))
					m_outputDistance->setValueAt(i + 1, j, candidate);
			} else {
				candidate = m_outputDistance->getValueAt(i, j)
						+ (sgn(m_outputDistance->getValueAt(i + 1, j)) * h);
				if (abs(candidate)
						< abs(m_outputDistance->getValueAt(i + 1, j)))
					m_outputDistance->setValueAt(i + 1, j, candidate);
			}
		}
	}

	for (int j = intersec_xmin; j < m_outputDistance->getMaxX(); j++) {
		for (int i = intersec_ymax - 1; i > m_outputDistance->getMinY(); i--) {
			if (j < intersec_xmax && i > intersec_ymin) {
				// calculate new distance candidate and assign if appropriate
				candidate = m_outputDistance->getValueAt(i, j)
						+ (sgn(m_inputDistance->getValueAt(i - 1, j)) * h);
				if (abs(candidate)
						< abs(m_outputDistance->getValueAt(i - 1, j)))
					m_outputDistance->setValueAt(i - 1, j, candidate);
			} else {
				candidate = m_outputDistance->getValueAt(i, j)
						+ (sgn(m_outputDistance->getValueAt(i - 1, j)) * h);
				if (abs(candidate)
						< abs(m_outputDistance->getValueAt(i - 1, j)))
					m_outputDistance->setValueAt(i - 1, j, candidate);
			}
		}
	}
	for (int i = intersec_ymin; i < m_outputDistance->getMaxY(); i++) {
		for (int j = intersec_xmin; j < m_outputDistance->getMaxX() - 1; j++) {
			// x-direction forward
			if (j < intersec_xmax - 1 && i < intersec_ymax) {
				if (m_inputDistance->getValueAt(i, j)
						* m_inputDistance->getValueAt(i, j + 1) <= 0.0) {
					// interpolate
					i_slope = (m_inputDistance->getValueAt(i, j + 1)
							- m_inputDistance->getValueAt(i, j)) / h;
					distToZero = -m_inputDistance->getValueAt(i, j) / i_slope;
					if (abs(m_outputDistance->getValueAt(i, j))
							> abs(distToZero))
						m_outputDistance->setValueAt(i, j,
								-distToZero * sgn(i_slope));
				}
				candidate = m_outputDistance->getValueAt(i, j)
						+ (sgn(m_inputDistance->getValueAt(i, j + 1)) * h);
				if (abs(candidate)
						< abs(m_outputDistance->getValueAt(i, j + 1)))
					m_outputDistance->setValueAt(i, j + 1, candidate);
			} else {
				candidate = m_outputDistance->getValueAt(i, j)
						+ (sgn(m_outputDistance->getValueAt(i, j + 1)) * h);
				if (abs(candidate)
						< abs(m_outputDistance->getValueAt(i, j + 1)))
					m_outputDistance->setValueAt(i, j + 1, candidate);
			}
		}
	}

	for (int i = intersec_ymin; i < m_outputDistance->getMaxY(); i++) {
		for (int j = intersec_xmax - 1; j > m_outputDistance->getMinX(); j--) {
			// x-direction outputDistanceward
			//check for sign change
			if (j > intersec_xmin && i < intersec_ymax) {
				// calculate new distance candidate and assign if appropriate
				candidate = m_outputDistance->getValueAt(i, j)
						+ (sgn(m_inputDistance->getValueAt(i, j - 1)) * h);
				if (abs(candidate)
						< abs(m_outputDistance->getValueAt(i, j - 1)))
					m_outputDistance->setValueAt(i, j - 1, candidate);
			} else {
				candidate = m_outputDistance->getValueAt(i, j)
						+ sgn(m_outputDistance->getValueAt(i, j - 1)) * h;
				if (abs(candidate)
						< abs(m_outputDistance->getValueAt(i, j - 1)))
					m_outputDistance->setValueAt(i, j - 1, candidate);
			}
		}
	}
	//	if(m_grainHandler->get_loop() > 1 ) executeExactRedist();
	m_outputDistance->clampValues(-m_grainHandler->delta,
			m_grainHandler->delta);

	m_inputDistance->resize(m_outputDistance->getMinX(),
			m_outputDistance->getMinY(), m_outputDistance->getMaxX(),
			m_outputDistance->getMaxY());
	// 	 set the references for the convolution step
	//	if(id==1) plot_box(true, 2, "after_resize", true);
}

/**************************************/
// end of redist
/**************************************/

/**************************************/
// plot the box and all its properties
/**************************************/
void LSbox::resizeGrid(int newSize) {

	int ngridpointsNew = newSize + 2 * m_grainHandler->get_grid_blowup();
	double h = m_grainHandler->get_h();
	double hn = 1.0 / (double) newSize;

	m_newXMin = m_outputDistance->getMinX() * (h / hn) + 0.5;
	m_newXMax = m_outputDistance->getMaxX() * (h / hn) + 0.5;
	m_newYMin = m_outputDistance->getMinY() * (h / hn) + 0.5;
	m_newYMax = m_outputDistance->getMaxY() * (h / hn) + 0.5;

	double xl, xr, yo, yu;

	double pointx, pointy;

	// resize to complete superposition
	//	while (minXnew * hn > outputDistance->getMinX() * h && minXnew > 0) {
	//		minXnew--;
	//	}
	//	while (minYnew * hn > outputDistance->getMinY() * h && minYnew > 0) {
	//		minYnew--;
	//	}
	//	while (maxXnew * hn < outputDistance->getMaxX() * h && maxXnew < ngridpointsNew) {
	//		maxXnew++;
	//	}
	//	while (maxYnew * hn < outputDistance->getMaxY() * h && maxYnew < ngridpointsNew) {
	//		maxYnew++;
	//	}
	m_inputDistance->resize(m_newXMin, m_newYMin, m_newXMax, m_newYMax);
	m_inputDistance->resizeToSquare(ngridpointsNew);
	//   plot_box(true, 2, "before_resize", true);

	for (int i = m_inputDistance->getMinY(); i < m_inputDistance->getMaxY();
			i++) {
		for (int j = m_inputDistance->getMinX(); j < m_inputDistance->getMaxX();
				j++) {
			pointx = j * (hn / h);
			pointy = i * (hn / h);

			xl = int(pointx);
			xr = int(pointx + 1);
			yo = int(pointy + 1);
			yu = int(pointy);

			if (xr > m_outputDistance->getMaxX() - 2
					|| yo > m_outputDistance->getMaxY() - 2
					|| yu < m_outputDistance->getMinY()
					|| xl < m_outputDistance->getMinX()) {
				m_inputDistance->setValueAt(i, j, -m_grainHandler->delta);
				continue;
			}
			double ro, ru, newDistVal;
			ro = 1 / (xr - xl)
					* ((xr - pointx) * m_outputDistance->getValueAt(yo, xl)
							+ (pointx - xl)
									* m_outputDistance->getValueAt(yo, xr));
			ru = 1 / (xr - xl)
					* ((xr - pointx) * m_outputDistance->getValueAt(yu, xl)
							+ (pointx - xl)
									* m_outputDistance->getValueAt(yu, xr));
			newDistVal = 1 / (yo - yu)
					* ((yo - pointy) * ru + (pointy - yu) * ro);
			if (newDistVal != newDistVal) {
				char waitbuffer;
				cerr << " nan " << endl;
				cin >> waitbuffer;
			}
			m_inputDistance->setValueAt(i, j, newDistVal);

		}
	}
	m_outputDistance->resize(m_newXMin, m_newYMin, m_newXMax, m_newYMax);
	m_outputDistance->resizeToSquare(ngridpointsNew);

	//	if(id==1) plot_box(true, 2, "after_resize", true);
	//plot_box for all boxes and compare with prior !

}

void LSbox::recalculateIDLocal()
{
	reizeIDLocalToDistanceBuffer();
	executeComparison();
}


double LSbox::computeMisorientation(LSbox* grain_2) {
	double result = 0.;
	if (Settings::LatticeType == 225 || Settings::LatticeType == 229 ) {
		result = m_orientationQuat->misorientationCubicQxQ(grain_2->m_orientationQuat);
	}
	if ( Settings::LatticeType == 194 ) {
		result = m_grainHandler->m_misOriHdl->calculateMisorientation_hexagonal(m_orientationQuat, grain_2->m_orientationQuat);
	}
	if (result >= DEG2RAD(1.)) {
		return result;
	}
	return DEG2RAD(1.);
}


double LSbox::computeMisorientation(unsigned int grainID)
{
	return computeMisorientation(m_grainHandler->getGrainByID(grainID));
}


double LSbox::GbEnergyModel(double theta, LSbox* candidate)
{
	double gamma_hagb = 1.; //relative value
	double theta_ref = DEG2RAD(15.);
	if (Settings::IsIsotropicNetwork == false) {
		if (theta > theta_ref) {
			return gamma_hagb;
		}
		if (theta < DEG2RAD(1.)) {
			theta = DEG2RAD(1.);
		}
		return gamma_hagb * (theta / theta_ref) * (1. - log(theta / theta_ref)); //TODO::catch NaN
	}
	return gamma_hagb;
}


double LSbox::GbMobilityModel(double thetaMis, LSbox* candidate)
{
	if (Settings::IsIsotropicNetwork == false) {
		//filter based on lattice typse// 8.66025 = 
		// thus 8.66 is the permissible deviation from the perfect Twinboundary SIGMA 3
		if (Settings::IdentifyTwins == true) {
			//Theta (15)* 1/ sqrt(SIGMA) here SiGMA is the number of coincidence points (3 for the Twinboundary)
			//TODO::but not all boundaries are coherent twin boundaries like assumed here
			if (MisoriToTwinBoundary(candidate) < DEG2RAD(8.66025)) {
				return 0.01;
			}
			else {
				return 1. - (0.99 * exp(-5.0 * (pow(thetaMis / DEG2RAD(15.), 9.))));
			}
		}
		else {
			return 1. - (0.99 * exp(-5.0 * (pow(thetaMis / DEG2RAD(15.), 9.))));
		}
	}
	return 1.;
}


bool LSbox::isNeighbour(LSbox* candidate) {
	return m_grainBoundary.isBoxDirectNeighbour(candidate);
}

void LSbox::constructBoundarySectors() {
	m_grainBoundary.buildBoundarySectors(m_IDLocal);
}

double LSbox::getWeight(int i, int j) {
	return m_grainBoundary.getWeight(i, j, m_IDLocal);
}

void LSbox::clearContourGrainArea() {
	if (m_grainBoundary.getBoundarySegmentCount() == 0)
		return;

	double minx = m_grainBoundary.getRawBoundary()[0].x, maxx =
			m_grainBoundary.getRawBoundary()[0].x, miny =
			m_grainBoundary.getRawBoundary()[0].y, maxy =
			m_grainBoundary.getRawBoundary()[0].y;
	for (unsigned int i = 0; i < m_grainBoundary.getRawBoundary().size(); i++) {
		minx = min(minx, m_grainBoundary.getRawBoundary()[i].x);
		maxx = max(maxx, m_grainBoundary.getRawBoundary()[i].x);

		miny = min(miny, m_grainBoundary.getRawBoundary()[i].y);
		maxy = max(maxy, m_grainBoundary.getRawBoundary()[i].y);
	}
	for (int i = miny - 1; i < maxy + 1; i++)
		for (int j = minx - 1; j < maxx + 1; j++) {
			if (m_inputDistance->getValueAt(i, j) >= 0) {
				m_inputDistance->setValueAt(i, j, -m_grainHandler->get_h());
			}
		}
}


double LSbox::getWeigthFromHandler(int i, int j) {
	return m_grainHandler->weightsMatrix[i][j];
}


LSbox* LSbox::getNeighbourAt(int i, int j) {
	if (i < 0 || i >= m_IDLocal.getMaxY() || j < 0 || j >= m_IDLocal.getMaxX())
		return NULL;
	else
		return m_grainHandler->getGrainByID(m_IDLocal.getValueAt(i, j).grainID);
}

void LSbox::outputMemoryUsage(ofstream& output) {
	output << m_inputDistance->getMaxX() - m_inputDistance->getMinX() << " "
			<< m_inputDistance->getMaxY() - m_inputDistance->getMinY() << " ";
	output << m_outputDistance->getMaxX() - m_outputDistance->getMinX() << " "
			<< m_outputDistance->getMaxY() - m_outputDistance->getMinY() << endl;
	output << m_IDLocal.getMaxX() - m_IDLocal.getMinX() << " "
			<< m_IDLocal.getMaxY() - m_IDLocal.getMinY() << endl;
	output << m_IDLocal.getMaxX() - m_IDLocal.getMinX() << " "
			<< m_IDLocal.getMaxY() - m_IDLocal.getMinY() << endl;

	output << m_outputDistance->getTotalMemoryUsed() << " "
			<< m_IDLocal.getTotalMemoryUsed() << endl;

}

vector<int> LSbox::getDirectNeighbourIDs() {
	return m_grainBoundary.getDirectNeighbours();
}

vector<double> LSbox::getGBLengths() {
	return m_grainBoundary.getGbSegmentLength();
}

void LSbox::computeDirectNeighbours(
		const RTree<unsigned int, int, 2, float>& tree) {
	int min[2], max[2];
	min[0] = getMinX();
	min[1] = getMinY();
	max[0] = getMaxX();
	max[1] = getMaxY();
	vector<unsigned int> intersectingGrains;
	tree.Search(min, max, intersectingGrains);
	for (unsigned int k = 0; k < intersectingGrains.size(); k++) {
		if (m_ID != intersectingGrains[k]) {
			m_secondOrderNeighbours.push_back(
					m_grainHandler->getGrainByID(intersectingGrains[k]));
			m_grainBoundary.addDirectNeighbourManual(
					m_grainHandler->getGrainByID(intersectingGrains[k]));
		}
	}
}
void LSbox::calculateMagneticEnergy() {
	double ND[3] = {0., 0., 1.};
	double *euler = new double[3];
	euler = m_orientationQuat->quaternion2Euler();
	double p1 = euler[0];
	double t = euler[1];
	double p2 = euler[2];
	double cosine = 0;

	//Rotated ND in order to represent the unit vector with c orientation
	double cAxis[3] = {0., 0., 0.};

	double rotMatrix[9] = { cos(p1) * cos(p2) - sin(p1) * sin(p2) * cos(t), sin(
			p1) * cos(p2) + cos(p1) * sin(p2) * cos(t), sin(p2) * sin(t), -cos(
			p1) * sin(p2) - sin(p1) * cos(p2) * cos(t), -sin(p1) * sin(p2)
			+ cos(p1) * cos(p2) * cos(t), cos(p2) * sin(t), sin(p1) * sin(t),
			-cos(p1) * sin(t), cos(t) };

	double trans[3][3];

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			trans[j][i] = rotMatrix[i * 3 + j];

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			cAxis[i] += trans[i][j] * ND[j];

	cosine = (Settings::MagneticParams.MagneticVector_x * cAxis[0])
			+ (Settings::MagneticParams.MagneticVector_y * cAxis[1])
			+ (Settings::MagneticParams.MagneticVector_z * cAxis[2]);
	if (cosine > 1.) {
		cerr << "cosine > 1 in calculate magneticEnergy" << "\n";
	}
	m_magneticEnergy = 0.5 * Settings::MagneticParams.VacuumPermeability * Settings::MagneticParams.deltaMagSys
			* SQR(Settings::MagneticParams.MagneticForceField) * SQR(cosine);
	m_magneticEnergy = m_magneticEnergy / Settings::HAGB_Energy * Settings::Physical_Domain_Size;
	delete[] euler;
}


double LSbox::MisoriToTwinBoundary(LSbox* candidate) {
	Quaternion newQuat;
	newQuat = m_orientationQuat->misorientationQuaternionCubic(
			candidate->m_orientationQuat);
	double angle = newQuat.misorientationCubicQxQ(m_grainHandler->TwinBoundary);
	return angle;
}


size_t LSbox::get_contour_vertices( vector<double> & vrts )
{
	//vertices
	//facets
	//absCoordinates
	//x = CLAMP((iterator.x - m_grainHandler->get_grid_blowup()) * m_grainHandler->get_h());
	//y = CLAMP((iterator.y - m_grainHandler->get_grid_blowup()) * m_grainHandler->get_h());
	//for ( const auto& iterator : m_grainBoundary.getRawBoundary()) {
	for ( vector<SPoint>::iterator it = m_grainBoundary.getRawBoundary().begin();
		it != m_grainBoundary.getRawBoundary().end(); it++ ) {
		vrts.push_back( it->x );
		vrts.push_back( it->y );
	}
	return m_grainBoundary.getRawBoundary().size();
	//aforementioned loops follows the sequence
	/*
	unsigned int prv = 0;
	unsigned int nxt = 1;
	size_t last = m_grainBoundary.getRawBoundary().size();
	for ( vector<SPoint>::iterator it = m_grainBoundary.getRawBoundary().begin();
		it != m_grainBoundary.getRawBoundary().end(); it++ ) {
		fcts.push_back(prv);
		fcts.push_back(nxt);
		prv++;
		fcts.push_back( fct_id );
	}
	*/
}


void LSbox::get_contour_xdmf_topology( const unsigned int vrts_offset, vector<unsigned int> & inds )
{
	inds.push_back( 3 ); //XDMF polygon topology
	inds.push_back( m_grainBoundary.getRawBoundary().size() );
	unsigned int i = 0;
	for ( const auto& iterator : m_grainBoundary.getRawBoundary()) {
		inds.push_back( vrts_offset + i ); //+1 missing ?
		i++;
	}
}


void LSbox::get_contour_xdmf_info( vector<double> & ifo )
{
	for (const auto& iterator : m_grainBoundary.getRawBoundary()) {
		ifo.push_back( iterator.energy * iterator.mob );
	}
}


vector<double> LSbox::get_quaternion() {
	/*
	double* newori = this->m_orientationQuat->quaternion2Euler();
	double phi1 = newori[0];
	double Phi = newori[1];
	double phi2 = newori[2];
	delete [] newori;
	return vector<double> { phi1, Phi, phi2 };
	*/
	return vector<double> {
		m_orientationQuat->get_q0(),
		m_orientationQuat->get_q1(),
		m_orientationQuat->get_q2(),
		m_orientationQuat->get_q3() };
}


vector<double> LSbox::get_barycentre() {
	double x = (getMaxX() - getMinX()) / 2 + getMinX();
	x *= m_grainHandler->get_h();
	double y = (getMaxY() - getMinY()) / 2 + getMinY();
	y *= m_grainHandler->get_h();
	return vector<double>{ x, y };
}

#define I_DONT_KNOW_YET		(-100.0)