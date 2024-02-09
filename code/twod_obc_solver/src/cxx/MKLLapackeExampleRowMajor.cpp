#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
//extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
//extern void print_int_vector( char* desc, MKL_INT n, MKL_INT* a );

/* Parameters */
#define N 3
#define NRHS 1
#define LDA N
#define LDB NRHS

/* Auxiliary routine: printing a matrix */
/*void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}
// Auxiliary routine: printing a vector of integers/
void print_int_vector( char* desc, MKL_INT n, MKL_INT* a ) {
        MKL_INT j;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
        printf( "\n" );
}
*/

void print_matrix_2( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
		MKL_INT j = 0; //always because b has only one column
		for( MKL_INT i = 0; i < m; i++ ) {
			printf( " %6.2f\n", a[i*lda+j] );
		}
}

int main(int argc, char *argv[]) {
		/* Locals */
		MKL_INT n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
		/* Local arrays */
		MKL_INT ipiv[N];
		double a[LDA*N] = {
			6.80, -6.05, -0.45,
			-2.11, -3.30,  2.58,
			5.66, 5.36, -2.70 
			//no solution example
			//2.f, -4.f, 6.f,
			//-1.f, 3.f, -2.f,
			//1.f, -2.f, 3.f
		};
		double b[LDB*N] = {
			4.02,
			6.19,
			-8.22
			//no solution example
			//5.f,
			//-1.f,
			//1.f
		};
		/* Executable statements */
		printf( "LAPACKE_dgesv (row-major, high-level) Example Program Results\n" );
		//Solve A*X = B */
		info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb );
		//Check for the exact singularity
		if( info > 0 ) { //no solution case
			//printf( "The diagonal element of the triangular factor of A,\n" );
			//printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
			printf( "the solution could not be computed.\n" );
			exit( 1 );
		}
		//see how to get results
		print_matrix_2( "Solution", n, nrhs, b, ldb );
		/* Print details of LU factorization */
		//print_matrix( "Details of LU factorization", n, n, a, lda );
		/* Print pivot indices */
		//print_int_vector( "Pivot indices", n, ipiv );
		exit( 0 );
   /* End of LAPACKE_dgesv Example */
//row major example
