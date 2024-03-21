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
#pragma once
#ifndef _MYMATH_H_
#define _MYMATH_H_

#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "Quaternion.h"


#define PI 3.14159265358979323846
#define QUAT2EUL_ETA 	1e-20

#define SWAP(A, B) {struct tempStruct { char C[sizeof(A)];} swap_tmp;\
    swap_tmp = *( struct tempStruct*) &A;\
    *( struct tempStruct*) &A = *( struct tempStruct*) &B;\
    *( struct tempStruct*) &B = swap_tmp;}


typedef double Real;

struct Point{
	double x;
	double y;
	double z;
};


class mathMethods
{
public:
	mathMethods( void );
	~mathMethods( void );

	void multiplyQuaternions( double *q, double* p, double* r );
	void multiplyQuaternions2( double *q, double* p, double* r ); //this function is what has been implemented in COReV2

	//calculate disorientation among two orientations in various parametrizations
	Quaternion* misorientationQuaternionCubic( double* p, double* q);
	double misorientationCubicQxQ( double q01, double q11, double q21, double q31, double q02, double q12, double q22, double q32 );
};

#endif
