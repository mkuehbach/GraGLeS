/*
 Utility math library
 Copyright (C) 2015  Markus Kuehbach

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
#include "MyMath.h"
#include "Assertions.h"


//###LB: Unless indicated, all functions in radians.
//###MK: Quaternion algebra: q = -q define an equivalent rotation.


mathMethods::mathMethods()
{
}

mathMethods::~mathMethods()
{
}


double mathMethods::misorientationCubicQxQ(double q01, double q11, double q21,
		double q31, double q02, double q12, double q22, double q32) {
	int i;

	Real p[4] = { q01, q11, q21, q31 };
	Real q[4] = { q02, q12, q22, q32 };

	Real qm1[4]; //Inverse of quaternion q

	for (i = 0; i < 4; i++) {
		qm1[i] = q[i];
		if (i > 0)
			qm1[i] *= -1;
	}
	//qm1[0] *= -1;

	Real r[4]; //Resulting quaternion, rotation of the two previous quaternions pq-1

	multiplyQuaternions(qm1, p, r);

	//Now, we have to determine the smallest angle.

	Real r0[6][4]; //There are 12 possible angles

	//Note: The notation r0 is due to the definition of the quaternion which lie
	//in the fundamental zone, this vector possesses the smallest angle, in such a way
	//that r0 is actually the scalar part of this quaternion

	double a = r[0];
	double b = r[1];
	double c = r[2];
	double d = r[3];

	Real fac = 1.0 / sqrt(2.0); //0.70710678;

	//six fundamental quaternions
	r0[0][0] = (a + b) * fac;
	r0[0][1] = (a - b) * fac;
	r0[0][2] = (c + d) * fac;
	r0[0][3] = (c - d) * fac;
	r0[1][0] = (a + c) * fac;
	r0[1][1] = (a - c) * fac;
	r0[1][2] = (b + d) * fac;
	r0[1][3] = (b - d) * fac;
	r0[2][0] = (a + d) * fac;
	r0[2][1] = (a - d) * fac;
	r0[2][2] = (b + c) * fac;
	r0[2][3] = (b - c) * fac;
	r0[3][0] = (a + b + c + d) * 0.5;
	r0[3][1] = (a + b - c - d) * 0.5;
	r0[3][2] = (a - b + c - d) * 0.5;
	r0[3][3] = (a - b - c + d) * 0.5;
	r0[4][0] = (a + b + c - d) * 0.5;
	r0[4][1] = (a + b - c + d) * 0.5;
	r0[4][2] = (a - b + c + d) * 0.5;
	r0[4][3] = (a - b - c - d) * 0.5;
	r0[5][0] = a;
	r0[5][1] = b;
	r0[5][2] = c;
	r0[5][3] = d;

	Real omega = 0.0;

	for (i = 0; i < 6; i++)
		for (int j = 0; j < 4; j++)
			if (fabs(r0[i][j]) > omega)
				omega = fabs(r0[i][j]);

	QUICKASSERT( omega < 1.01 );

	if (omega > 1.0) //avoid singularity of acos function
		omega = (Real) (int) omega;

	omega = 2 * acos(omega);
	QUICKASSERT( omega <= 1.099 );
	return omega;
}


void mathMethods::multiplyQuaternions(double *q, double* p, double* r) {
	//mathematically multiplies quaternions q and p, active or passive rotation convention does not matter
	//verified via http://mathworld.wolfram.com/Quaternion.html as well as http://www.mathworks.de/de/help/aeroblks/quaternionmultiplication.html
	//equivalent to the multiplication given in Grimmer, H, 1974 Acta Cryst A30, 685-688
	r[0] = +q[0] * p[0] - q[1] * p[1] - q[2] * p[2] - q[3] * p[3];
	r[1] = +q[1] * p[0] + q[0] * p[1] - q[3] * p[2] + q[2] * p[3];
	r[2] = +q[2] * p[0] + q[3] * p[1] + q[0] * p[2] - q[1] * p[3];
	r[3] = +q[3] * p[0] - q[2] * p[1] + q[1] * p[2] + q[0] * p[3];
}

void mathMethods::multiplyQuaternions2(double *p, double *q, double* r) {
	//what has been implemented in misorientationCubic in COReV2 is not a the standard definition of a quaternion multiplication!!!
	//but with qm1 as -q0, q1, q2, q3 applied to p it is the same as MTex
	r[0] = +p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
	r[1] = +p[1] * q[0] + p[0] * q[1] + p[3] * q[2] - p[2] * q[3];
	r[2] = +p[2] * q[0] - p[3] * q[1] + p[0] * q[2] + p[1] * q[3];
	r[3] = +p[3] * q[0] + p[2] * q[1] - p[1] * q[2] + p[0] * q[3];
}
