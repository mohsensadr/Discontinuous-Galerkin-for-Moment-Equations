/* -------------------------------------------------------------------- */
/*      Example program to show the use of the "PARDISO" routine        */
/*      on for unsymmetric linear systems                               */
/* -------------------------------------------------------------------- */
/*      This program can be downloaded from the following site:         */
/*      http://www.pardiso-project.org                                  */
/*                                                                      */
/*  (C) Olaf Schenk, Institute of Computational Science                 */
/*      Universita della Svizzera italiana, Lugano, Switzerland.        */
/*      Email: olaf.schenk@usi.ch                                       */
/* -------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "papi.h"
#include <omp.h>


#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif

/* PARDISO prototype. */
#ifdef __cplusplus
extern "C" {
#endif

int solving_phi(int *ia, int *ja, double *a, int n, double *b, double *x, int ) ;

/* PARDISO prototype. */
void pardisoinit (void   *, MKL_INT    *,   MKL_INT *, MKL_INT *, double *, MKL_INT *);
void pardiso     (void   *, MKL_INT    *,   MKL_INT *, MKL_INT *,    MKL_INT *, MKL_INT *, 
                  double *, MKL_INT    *,    MKL_INT *, MKL_INT *,   MKL_INT *, MKL_INT *,
                     MKL_INT *, double *, double *, MKL_INT *, double *);
void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec     (int *, int *, double *, int *);
void pardiso_printstats (int *, int *, double *, int *, int *, int *, double *, int *);
}


void PardisoSolve_original(MKL_INT mtype, MKL_INT n, MKL_INT *ia, MKL_INT *ja, double *a, int nnz,
										 double *b, double *x, void *pt, MKL_INT phase, MKL_INT *iparm)
{
  float  rtime, ptime, mflops;
  long long flpops, flpops_old;
  double start, end;
  int retval;


//    MKL_INT     n =  215712;
//    int      nnz = 23163773;
//    MKL_INT  *ia = (MKL_INT *) malloc((n+1)*sizeof(MKL_INT));
//    MKL_INT  *ja = (MKL_INT *)malloc((nnz)*sizeof(MKL_INT));
//    double *a = (  double *)malloc((nnz)*sizeof( double));
//    double *b = (  double *)malloc((n)*sizeof( double));
//    double *x = ( double *)malloc((n)*sizeof( double));


MKL_INT i;


    //int      mtype = 11;        /* Real unsymmetric matrix */

    /* RHS and solution vectors. */

    MKL_INT      nrhs = 1;          /* Number of right hand sides. */

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
        /* Set right hand side to i. */


    
   // void    *pt[64];

    /* Pardiso control parameters. */
    //int      iparm[64];
    double   dparm[64];
    MKL_INT      solver;
    MKL_INT      maxfct, mnum, error, msglvl;

    /* Number of processors. */
    MKL_INT      num_procs;

    /* Auxiliary variables. */
    char    *var;
    MKL_INT       k;

    double   ddum;              /* Double dummy */
    MKL_INT      idum;              /* Integer dummy. */

/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters and initialize the solvers      */
/*     internal adress pointers. This is only necessary for the FIRST   */
/*     call of the PARDISO solver.                                      */
/* ---------------------------------------------------------------------*/
      
    error = 0;
    solver = 0; /* use sparse direct solver */


    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
        // return 1;
    }
    else
        printf("[PARDISO]: License check was successful ... \n");
 

    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    
    iparm[1] = 2; // Fill-in reordering from METIS  
    iparm[2]  = num_procs;
    iparm[7] = 200; // Max numbers of iterative refinement steps  
    iparm[9] = 9;//10 // Perturb the pivot elements with 1E-13  
    iparm[10] = 1; /* no scaling  */
    iparm[12] = 0; /* no matching */
    
 
	    iparm[3] = 0; // No iterative-direct algorithm  
	    iparm[4] = 0; // No user fill-in reducing permutation  
	    iparm[5] = 0; // 0 = Write solution into x  
	    iparm[6] = 0; // Not in use  

	    iparm[8] = 0; // Not in use  

 
	    iparm[11] = 0; // Not in use  
  
	    iparm[13] = 0; // Output: Number of perturbed pivots  
	    iparm[14] = 0; // Not in use  
	    iparm[15] = 0; // Not in use  
	    iparm[16] = 0; // Not in use  
	    iparm[17] = -1; // Output: Number of nonzeros in the factor LU  
	    iparm[18] = -1; // Output: Mflops for LU factorization  
	    iparm[19] = 0; // Output: Numbers of CG Iterations  



    maxfct = 1;         /* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 0;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */


/* -------------------------------------------------------------------- */    
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */ 
    for (i = 0; i < n+1; i++) {
        ia[i] += 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] += 1;
    }

/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

/*    
    pardiso_chkmatrix  (&mtype, &n, a, ia, ja, &error);
    if (error != 0) {
        printf("\nERROR in consistency of matrix: %d", error);
        exit(1);
    }
*/

/* -------------------------------------------------------------------- */
/* ..  pardiso_chkvec(...)                                              */
/*     Checks the given vectors for infinite and NaN values             */
/*     Input parameters (see PARDISO user manual for a description):    */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */
/*
    pardiso_chkvec (&n, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR  in right hand side: %d", error);
        exit(1);
    }
*/
/* -------------------------------------------------------------------- */
/* .. pardiso_printstats(...)                                           */
/*    prints information on the matrix to STDOUT.                       */
/*    Use this functionality only for debugging purposes                */
/* -------------------------------------------------------------------- */
/*
    pardiso_printstats (&mtype, &n, a, ia, ja, &nrhs, b, &error);
    if (error != 0) {
        printf("\nERROR right hand side: %d", error);
        exit(1);
    }
 */

/* -------------------------------------------------------------------- */    
/* ..  Reordering and Symbolic Factorization.  This step also allocates */
/*     all memory that is necessary for the factorization.              */
/* -------------------------------------------------------------------- */ 
    phase = 11; 
//fprintf(gp, "\n#%19s%20s%20s%20s","step","time", "mflops", "iparm[18]");
//if((retval=PAPI_flops( &rtime, &ptime, &flpops, &mflops))<PAPI_OK)
//	printf("\n OOPS!");
//start = omp_get_wtime();
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//end = omp_get_wtime();
//PAPI_flops( &rtime, &ptime, &flpops, &mflops );
//fprintf(gp, "\n%20c%20lf%20f%20d",'1', end-start, mflops,  iparm[18]);
//printf("\nstep= 1 time= %lf mflops= %f iparm[18]= %d", end-start, mflops,  iparm[18]);
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
   
/* -------------------------------------------------------------------- */    
/* ..  Numerical factorization.                                         */
/* -------------------------------------------------------------------- */    
    phase = 22;
iparm[18] = -1;
//PAPI_flips( &rtime, &ptime, &flpops, &mflops ); 
//start = omp_get_wtime();
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
//end = omp_get_wtime();
//PAPI_flops( &rtime, &ptime, &flpops, &mflops );
//fprintf(gp, "\n%20c%20lf%20f%20d",'2', end-start, mflops,  iparm[18]);
//printf("\nstep= 2 time= %lf mflops= %f iparm[18]= %d", end-start, mflops,  iparm[18]);
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    printf("\nFactorization completed ...\n ");

/* -------------------------------------------------------------------- */    
/* ..  Back substitution and iterative refinement.                      */
/* -------------------------------------------------------------------- */    
    phase = 33;

iparm[18] = -1;
//PAPI_flips( &rtime, &ptime, &flpops, &mflops ); 
//start = omp_get_wtime();
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);
//end = omp_get_wtime();
//PAPI_flops( &rtime, &ptime, &flpops, &mflops );
//fprintf(gp, "\n%20c%20lf%20f%20d",'3', end-start, mflops,  iparm[18]);
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }

    printf("\nSolve completed ... ");







/* -------------------------------------------------------------------- */    
/* ..  Convert matrix back to 0-based C-notation.                       */
/* -------------------------------------------------------------------- */ 
    for (i = 0; i < n+1; i++) {
        ia[i] -= 1;
    }
    for (i = 0; i < nnz; i++) {
        ja[i] -= 1;
    }

/* -------------------------------------------------------------------- */    
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */

 
//    phase = -1;                 /* Release internal memory. */
//iparm[18] = -1;
//PAPI_flips( &rtime, &ptime, &flpops, &mflops ); 
//start = omp_get_wtime();
//    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
//             &n, &ddum, ia, ja, &idum, &nrhs,
//             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
//end = omp_get_wtime();
//PAPI_flops( &rtime, &ptime, &flpops, &mflops );
//fprintf(gp, "\n%20c%20lf%20f%20d",'4', end-start, mflops,  iparm[18]);
//printf("\nstep= 4 time= %lf mflops= %f iparm[18]= %d", end-start, mflops,  iparm[18]);
//}
//fclose(gp);
//fp = fopen("results.txt", "w");
//    for(i=0; i<n; i++)
//    {
//      		fprintf(fp, "%lf\n", x[i]);
//    }
//    fclose(fp);
//    printf("results wrote in results.txt\n");
//    return 0;
} 
