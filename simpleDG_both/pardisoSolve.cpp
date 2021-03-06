#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mkl.h"

//#define  MKL_INT int
//#define INT int

//void PardisoSolve(long long mtype, long long n, long long *ia, long long *ja, void *A, long long nz,
//										 void *b, void *x, void *pt, long long phase, long long *iparm)
void PardisoSolve(MKL_INT mtype, MKL_INT n, MKL_INT *ia, MKL_INT *ja, void *A, int nz,
										 void *b, void *x, void *pt, MKL_INT phase, MKL_INT *iparm)
{  
	MKL_INT i, nrhs;
	MKL_INT maxfct, mnum, msglvl, error;       
	MKL_INT idum;    // Integer dummy 
	double ddum;  // Double dummy 
	nrhs = 1;    // number of rhs
	maxfct = 1;  // Maximum number of numerical factorizations
	mnum = 1;    // Which factorization to use.
	msglvl = 0;  // 1 = Print statistical information in file 
	error = 0;   // Initialize error flag 
  
/*  MKL_INT *ia, *ja;
  ia = (MKL_INT *)malloc(sizeof(MKL_INT)*nz);
  ja = (MKL_INT *)malloc(sizeof(MKL_INT)*nz);
  for( int i=0; i<nz; i++ ) {
    ia[i] = (MKL_INT)ia0[i];
    ja[i] = (MKL_INT)ja0[i];
  };*/



  printf( " PardisoSolve :" );
  switch( phase ) 
  {
    case( -11 ):
      printf( " Setup Phase \n" );
      /*
      printf("\nINTERNAL: The matrix is given by: ");fflush(stdout);
      for (i = 0; i < nz; i+=1) {
        printf("\n i=%d, j=%d, A(i,j)=%f", ia[i], ja[i], ((double*)A)[i] );fflush(stdout);
      };	printf ("\n");
      */
      
	    // Pardiso control parameters.  
    	for (i = 0; i < 64; i++) iparm[i] = 0;             

	    iparm[0] = 1; // No solver default  
	    iparm[1] = 3;//3; // Fill-in reordering from METIS  
	    iparm[2] = 1; // Numbers of processors, value of OMP_NUM_THREADS  
	    iparm[3] = 0; // No iterative-direct algorithm  
	    iparm[4] = 0; // No user fill-in reducing permutation  
	    iparm[5] = 0; // 0 = Write solution into x  
	    iparm[6] = 0; // Not in use  
	    iparm[7] = 200; // Max numbers of iterative refinement steps  
	    iparm[8] = 0; // Not in use  
	    iparm[9] = 10; // Perturb the pivot elements with 1E-13  
	    iparm[10] = 1; // 1 = Use nonsymmetric permutation and scaling MPS  
	    iparm[11] = 0; // Not in use  
	    iparm[12] = 0; // 1 = Maximum weighted matching algorithm is switched-on (default for non-symmetric)  
	    iparm[13] = 0; // Output: Number of perturbed pivots  
	    iparm[14] = 0; // Not in use  
	    iparm[15] = 0; // Not in use  
	    iparm[16] = 0; // Not in use  
	    iparm[17] = -1; // Output: Number of nonzeros in the factor LU  
	    iparm[18] = -1; // Output: Mflops for LU factorization  
	    iparm[19] = 0; // Output: Numbers of CG Iterations  
	    iparm[34] = 1; // 0= one-based indices, 1= zero-besed indices  
//	    iparm[60] = 1; // added by mohsen, due to problem with RAM
//	    iparm[59] = 1;
//	    iparm[51] = 16;

    break;
    case( 11 ):    // Symbolic Factorization
      printf( " Symbolic Factorization \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (MKL_INT *)ia, (MKL_INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, &ddum, &ddum, &error);
	    if (error != 0) {
		    printf("\nERROR during symbolic factorization: %d\n", error);
		    exit(1);
	    }
    break;
    case( 22 ):    // Numerical Factorization
      printf( " Numerical Factorization \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (MKL_INT *)ia, (MKL_INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, &ddum, &ddum, &error);
	    if (error != 0) {
		    printf("\nERROR during numerical factorization: %d\n", error);
		    exit(2);
	    }
    break;
    case( 33 ):   //Back substitution and iterative refinement
      printf( " Back Substitution and Iterative Refinement \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (MKL_INT *)ia, (MKL_INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, b, x, &error);
	    if (error != 0) {
		    printf("\nERROR during solution: %d\n", error);
		    exit(3);
	    }
    break;
    case( -1 ): // Release internal memory. 
      printf( " Release Internal Memory \n" );
	    pardiso(pt, &maxfct, &mnum, &mtype, &phase, 
        &n, A, (MKL_INT *)ia, (MKL_INT *)ja, &idum, &nrhs,
		    iparm, &msglvl, &ddum, &ddum, &error);
		break;
	};
};

/*************************************************************************************

#include <iostream>
#include <Eigen/SparseCore>
using Eigen::SparseMatrix;
using Eigen::MatrixXd;
typedef SparseMatrix<double,Eigen::RowMajor> SpMatrix;
using namespace Eigen;
using namespace std;


// Simple Test Routine
int pardisoTest( void ) 
{ 
  int n = 5;

  MatrixXd b(n,1), x(n,1);
  b(0) = 1;
  b(1) = 2;
  b(2) = 3;
  b(3) = 4;
  b(4) = 5;
  
	cout << "\nThe right hand side is:\n" << b << endl;

  SpMatrix SpM(n,n);
  SpM.coeffRef(0,1) = 1.0; 
  SpM.coeffRef(1,1) = -2.0;
  SpM.coeffRef(1,2) = 2.0;
  SpM.coeffRef(2,0) = 1.0;
  SpM.coeffRef(2,2) = -1.0;
  SpM.coeffRef(3,3) = -10.0;
  SpM.coeffRef(3,2) = -3.0;
  SpM.coeffRef(4,4) = 5.0;
  SpM.makeCompressed();
  
	cout << "\nThe sparse matrix is given by: \n";
  for (int k=0; k<SpM.outerSize(); ++k)
    for (SpMatrix::InnerIterator it(SpM,k); it; ++it ) {
      cout << "row: " << it.row() << " column: " << it.col() << " value: " << it.value() << endl;  
    };
  cout << "number of non-zeros " << SpM.nonZeros() << endl;
    
	int iparm[64];
	int *pt[64];
  for(int i=0; i<64; i++ ) pt[i]=0;
  PardisoSolve(11, SpM.rows(), SpM.outerIndexPtr(), SpM.innerIndexPtr(), SpM.valuePtr(), SpM.nonZeros(), 
               &b(0), &x(0), pt, -11, iparm);
  PardisoSolve(11, SpM.rows(), SpM.outerIndexPtr(), SpM.innerIndexPtr(), SpM.valuePtr(), SpM.nonZeros(), 
               &b(0), &x(0), pt, 11, iparm);
  PardisoSolve(11, SpM.rows(), SpM.outerIndexPtr(), SpM.innerIndexPtr(), SpM.valuePtr(), SpM.nonZeros(), 
               &b(0), &x(0), pt, 22, iparm);
  PardisoSolve(11, SpM.rows(), SpM.outerIndexPtr(), SpM.innerIndexPtr(), SpM.valuePtr(), SpM.nonZeros(), 
               &b(0), &x(0), pt, 33, iparm);  
  
	cout << "\nThe solution of the system is:\n" << x << endl;

	return 0;
}
/*************************************************************************************/


