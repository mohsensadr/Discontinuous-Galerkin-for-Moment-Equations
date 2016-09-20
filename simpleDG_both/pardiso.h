// -------------------------------------------------------------------- 
// Define PARDISO prototype and matrix type. 

typedef void * _MKL_DSS_HANDLE_t;
typedef int    _INTEGER_t;

void PARDISO(
             _MKL_DSS_HANDLE_t, _INTEGER_t *, _INTEGER_t *,
             _INTEGER_t *, _INTEGER_t *, _INTEGER_t *,
             void *,  _INTEGER_t *,  _INTEGER_t *,
             _INTEGER_t *, _INTEGER_t *, _INTEGER_t *,
             _INTEGER_t *, void *, void *, _INTEGER_t *);


typedef struct {	
  double re;	
  double i;
} COMPLEX;

// Done: Define PARDISO prototype. 
// -------------------------------------------------------------------- 


// simplified interface to pardiso
void modPardisoSolve(int mtype, int n, int *ia, int *ja, void *A, int nz,
										 void *b, void *x, void *pt, int phase, int *iparm);
void PardisoSolve_original(MKL_INT mtype, MKL_INT n, MKL_INT *ia, MKL_INT *ja, void *a, int nnz,
										 void *b, void *x, void *pt, MKL_INT phase, MKL_INT *iparm);



