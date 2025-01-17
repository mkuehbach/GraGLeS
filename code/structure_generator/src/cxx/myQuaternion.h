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

//These functionlities were implemented by Christian Miessen.


#define SWAP(A, B) {struct tempStruct { char C[sizeof(A)];} swap_tmp;\
    swap_tmp = *( struct tempStruct*) &A;\
    *( struct tempStruct*) &A = *( struct tempStruct*) &B;\
    *( struct tempStruct*) &B = swap_tmp;}
#define PI 3.14159265358979323846


//IPFZ Coloring
#define FAIL_MYMATH_NUMERIC		(-1.0)
#define EPS_PROXIMITY_IPFZ		(0.01)
#define IPF_COLOR_STRETCH_R		(0.5)
#define IPF_COLOR_STRETCH_G		(0.5)
#define IPF_COLOR_STRETCH_B		(0.5)
#define EPS_ENVIRONMENT			(1e-7)
#define RGBRANGE				255
#define SYMMETRIES_IN_CUBIC		24
#define DOUBLE_ACCURACY 		(5e-16)


#ifndef myQuaternion_H_
#define myQuaternion_H_

class mathMethods;
class randomClass;

//myQuaternionen class
class myQuaternion {
private:
	double q0;
	double q1;
	double q2;
	double q3;

	double x;	//position in stereographic standard triangle
	double y;
public:
	myQuaternion(void);
	myQuaternion(double q0, double q1, double q2, double q3);
	myQuaternion(double alpha, double x, double y, double z, bool radiant);
	~myQuaternion(void);
	myQuaternion& operator=(const myQuaternion &rhs);
	//myQuaternion& myQuaternion::operator+ (Quarternion const& lhs, Quarternion const& rhs);
	//myQuaternion& myQuaternion::operator- (Quarternion const& lhs, Quarternion const& rhs);
	myQuaternion operator*(const myQuaternion &rhs);

	bool operator ==(myQuaternion const& rhs);
	bool operator !=(myQuaternion const& rhs);

	myQuaternion Inverse();
	double getNorm(void);
	void Normalize();
	void Invert();
	void Sort(void);
	myQuaternion misorientationQuaternionCubic(myQuaternion* p);
	double rotationAngleBetweenQuaternions(myQuaternion *p);
	void euler2Quaternion(double *euler);
	double* Quaternion2EulerConst(void) const;
	double* Quaternion2Euler(void);
	void randomOriShoemakeQuat(randomClass& r);
	double misorientationCubicQxQ(myQuaternion* p);
	myQuaternion* specificallyDisorientednewOriFromReference(
			double sigma_rayl_dis2bunge, randomClass& mymath);
	//myQuaternion* newOrientationFromReference( randomClass& mymath);
	myQuaternion* randomMisorientationShoemake( double theta, randomClass& thernd );
	myQuaternion* randomRotationSO3Shoemake( randomClass& thernd );

	//IPF coloring for cubic in stereographic projection
	void QuatOnVector3D( double *q, double* v, double* r );
	void multiplyQuaternions( double *q, double* p, double* r);
	void project2fundamentalregion_ipfz ( void );
	void quat2ipfz( unsigned char *rgb );




	inline double get_q0() {
		return q0;
	}
	;
	inline double get_q1() {
		return q1;
	}
	;
	inline double get_q2() {
		return q2;
	}
	;
	inline double get_q3() {
		return q3;
	}
	inline double get_x() {
		return x;
	}
	inline double get_y() {
		return y;
	}
	;
};

#endif /* myQuaternion_H_ */
