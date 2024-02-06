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

#ifndef _random_h_
#define _random_h_

#include <math.h>
#include <time.h>
#include <omp.h>
#include <stdint.h>
#include <stdlib.h>
#include <float.h>

//Park Miller constants
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


#define _PI_ 3.1415926535897932384626433832795
#define EPSILON 0.01
#define EPS1 0.001
#define EPS2 1.0e-8
#define SQR(a) ((a)*(a))
#define CUBE(a) ((a)*(a)*(a))
#define MIN(X,Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X,Y) (((X) > (Y)) ? (X) : (Y))
#define MINIMUM_ANGULAR_SPREAD 8.726646259971647884618e-3


//Marsaglia's generator constants
#define R4LEN		128
#define R4DN		(3.442619855899)
#define R4M1		(2147483648.0)
#define R4TN		(3.442619855899)
#define R4VN		(9.91256303526217e-03)
#define R4R			(3.442620)


//Mersenne Twister Parameters
#define NMT 624
#define MMT 397
#define MATRIX_A 0x9908b0df
#define UPPER_MASK 0x80000000
#define LOWER_MASK 0x7fffffff

#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y) (y >> 11)
#define TEMPERING_SHIFT_S(y) (y << 7)
#define TEMPERING_SHIFT_T(y) (y << 15)
#define TEMPERING_SHIFT_L(y) (y >> 18)

#define RANDMAX (1.0-DBL_EPSILON)

class randomClass
{
public:
	randomClass(long seed = 0);
	~randomClass();

	void initPM(long sd) { seed = sd; }
	void initR4Uni( uint32_t sd ) { jsr_r4 = sd; }
	void initSHR3( uint32_t sd ) { jsr_shr3 = sd; }
	void initMT( uint32_t sd ) { seedMT = sd; sgenrand(); warmupMT( 500000 ); }

	void r4_nor_setup( void );


	double parkMiller( void );
	double r4_nor( double mu, double sigma );
	float r4_nor( void );
	float r4_uni( void );
	uint32_t shr3_seeded ( void );

	double MersenneTwister( void );
	void sgenrand( void );
	void warmupMT( unsigned int n ) {
		double toss = 0.0;
		for ( unsigned int r = 0; r < n; r++ ) {toss = this->MersenneTwister();}
	};

  private:
		long seed;			//Park-Miller minimal generator
		long iy;
		long* iv;

		uint32_t jsr_shr3;	//seed for Marsaglia's shr3 generator
		uint32_t jsr_r4;	//seed for Marsaglia's r4_uni generator

		uint32_t* kn;		//Ziggurat state set for Marsaglia's Ziggurat standard normal distributed generator
		float* fn;
		float* wn;

		uint32_t seedMT;	//state set for the MersenneTwister
		uint32_t* mt;
		uint32_t mti;
		uint32_t mag01[2] = {0x0, MATRIX_A};
};

typedef randomClass *randomClassP;

#endif
