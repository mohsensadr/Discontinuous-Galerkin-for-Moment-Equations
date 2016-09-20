#include<stdio.h>
#include <iostream>
#include <math.h>
//#include "mkl.h"

#include "EigenSetup.h"
#include "Mesh.h"
#include "System.h"
#include "Numerics.h"

#include "Tools.h"

#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT int
#endif
#if defined(MKL_ILP64)
#define which 1
#else
#define which 2
#endif

/**************************** PARDISO SOLVE INTERFACE ******************
simple sparse Matrix-Data:
(n,ia,ja,A,nz) = (dimension,outer-indices,innner-indices,values,#non-zeros)
b: right hand side will be over-written with solution
x: solution (output)
/**********************************************************************/
//void PardisoSolve(MKL_INT mtype, MKL_INT n, MKL_INT *ia, MKL_INT *ja, void *A, int nz,
//		 void *b, void *x, void *pt, MKL_INT phase, MKL_INT *iparm);
void PardisoSolve(MKL_INT mtype, MKL_INT n, MKL_INT *ia, MKL_INT *ja, void *A, int nz,
										 void *b, void *x, void *pt, MKL_INT phase, MKL_INT *iparm);
		
void PardisoSolve_original(MKL_INT mtype, MKL_INT n, MKL_INT *ia, MKL_INT *ja, double *a, int nnz,
										 double *b, double *x, void *pt, MKL_INT phase, MKL_INT *iparm);

void Numerics::InfoText( char *text )
{
  sprintf( text, "Method = DG, order = %d, nBasis = %d, nTot = %d, nnz = %ld, nnzRatio = %.2f, res = %.2e", 
	   order, nBasis, nTot, nnz, 1.0*nnz/nTot, res );
	   
};

		 
Numerics::Numerics( MeshCore *mesh, System *sys, int Order )
{
  startTimer(">> Numerics Setup\n");
  Sys = sys;
  Mesh = mesh;  
 
  order = Order;
 
  if( Mesh->Type == 3 ) nBasis = (order+1)*(order+2)/2;
  if( Mesh->Type == 4 ) nBasis = (order+1)*(order+1);
  nTot = Mesh->Elements.size()*Sys->nEqn*nBasis;
  
  nnz = 0;
  res = 0.0;
  stopTimer();
  
  startTimer(" allocate stencil size\n");
  VectorXi StencilSize(nTot);// the knowledge of the sparsity pattern accelerates the matrix build-up
  StencilSize.setZero(nTot);
  StencilSize = VectorXi::Constant(nTot,4*Sys->nEqn*nBasis);// precise band width needs to be improved!!
  stopTimer();
  
  startTimer(" allocate sparse matrices \n");
  Atotal.resize(nTot,nTot);
  Atotal.reserve(StencilSize);
  sol.resize(nTot,1);
  sol.setZero(nTot,1);
  rhs.resize(nTot,1);
  rhs.setZero(nTot,1);
  stopTimer();  
  
  cout << "Order of Elements " << order << "\n";
  cout << "Number of basis functions " << nBasis << "\n";
  cout << "Total number of unknowns " << nTot << "\n"; 
}

void Numerics::Assemblation()
{
  startTimer(">> Assemblation - element-wise\n"); 
  
  for (Range::iterator elem = Mesh->Elements.begin(); elem != Mesh->Elements.end(); elem++) { 
    std::vector<EntityHandle> vertices;
    int type;
    Mesh->get_connectivity(&(*elem), 1, vertices, true );
    type = vertices.size();
    vertices.resize(type+1);
    vertices[type] = vertices[0];
    
    MatrixXd pos(3,nBasis);
    Mesh->get_coords(&(vertices[0]),type,&pos(0));
    LagrangePts( pos, order, type );
    
    MatrixXd AssemMatX = AssembleDiffX(pos,order,type);
    MatrixXd AssemMatY = AssembleDiffY(pos,order,type);
    MatrixXd AssemMass = AssembleMass(pos,order,type);
    
    int offset = Mesh->get_ID(*elem)*Sys->nEqn*nBasis;
    
    for( int i=0; i<nBasis; i++ ) {  // To Do: switch i <-> j
      for( int j=0; j<nBasis; j++ ) {
	
      	for( int k=0; k<Sys->Ax.outerSize(); ++k)
      	  for (SpMatrix::InnerIterator it(Sys->Ax,k); it; ++it ) 
      	    Atotal.coeffRef(offset+it.row()*nBasis+i,offset+it.col()*nBasis+j) += it.value()*AssemMatX(i,j);
    	  for( int k=0; k<Sys->Ay.outerSize(); ++k)
    	    for (SpMatrix::InnerIterator it(Sys->Ay,k); it; ++it ) 
    	      Atotal.coeffRef(offset+it.row()*nBasis+i,offset+it.col()*nBasis+j) += it.value()*AssemMatY(i,j);
  	    for( int k=0; k<Sys->P.outerSize(); ++k)
  	      for (SpMatrix::InnerIterator it(Sys->P,k); it; ++it ) 
        		Atotal.coeffRef(offset+it.row()*nBasis+i,offset+it.col()*nBasis+j) += it.value()*AssemMass(i,j);
      	      
  	    MatrixXd pos0 = pos.col(j);
  	    MatrixXd f = Sys->Force(pos0);
    	  for( int k=0; k<Sys->nEqn; k++)
    	    rhs(offset+k*nBasis+i) += AssemMass(i,j)*f(k);
      };  
    };    
    
    VectorXi ID(type);
    Mesh->tag_get_data(Mesh->ID,&vertices[0],type,&ID(0));
    
    int nBasisE = order+1;
    MatrixXi idx(type,nBasisE),idxN(type,nBasisE);
    if( (order == 2) && (type == 3) ) {
      idx << 0,1,3, 1,2,4, 2,0,5;
      idxN << 0,2,5, 1,0,3, 2,1,4;    
    };
    if( (order == 1) && (type == 3) ) {
      idx << 0,1, 1,2, 2,0;
      idxN << 0,2, 1,0, 2,1;    
    };
    if( (order == 2) && (type == 4) ) {
      idx <<  0,1,4, 1,2,5, 2,3,6, 3,0,7;
      idxN << 0,3,7, 1,0,4, 2,1,5, 3,2,6;    
    };
    if( (order == 1) && (type == 4) ) {
      idx <<  0,1, 1,2, 2,3, 3,0;
      idxN << 0,3, 1,0, 2,1, 3,2;    
    };
    
    for( int nn = 0; nn < type; nn++ ) {
      MatrixXd posE(3,nBasisE);
      for( int ix=0; ix < nBasisE; ix++ ) posE.col(ix) = pos.col(idx(nn,ix));
      
      MatrixXd normal = cross(posE.col(0)-posE.col(1));
      MatrixXd Am = Sys->Aminus(normal);
      MatrixXd AssemMassE = AssembleMassE(posE,order);
      
      Range neighbor;
      Mesh->get_adjacencies(&(vertices[nn]), 2, 2, true, neighbor );
      neighbor -= Range(*elem,*elem);
      
      int offsetN =  Mesh->get_ID(neighbor.front())*Sys->nEqn*nBasis;
            
      if(neighbor.size() == 1 ) {
        const EntityHandle *connectN;
        int num_connectN;	
        Mesh->get_connectivity(neighbor.front(), connectN, num_connectN, true);
        VectorXi IDN(type);                                           // assuming the neighbar has the same type!
        Mesh->tag_get_data(Mesh->ID,connectN,num_connectN,&IDN(0));
        int pp = 0; while( pp < type && ID[nn] != IDN[pp]) { pp++; }
        
        for( int i=0; i<nBasisE; i++ )
          for( int j=0; j<nBasisE; j++ )
            for( int k=0; k<Sys->nEqn; k++ )
              for( int l=0; l<Sys->nEqn; l++ ) {
                Atotal.coeffRef(offset+k*nBasis+idx(nn,i),offset+l*nBasis+idx(nn,j)) += 0.5*Am(k,l)*AssemMassE(i,j);
                Atotal.coeffRef(offset+k*nBasis+idx(nn,i),offsetN+l*nBasis+idxN(pp,j)) -= 0.5*Am(k,l)*AssemMassE(i,j);
              };
      };
      
      if(neighbor.size() == 0 ) {
        Range edge;
        Mesh->get_adjacencies(&(vertices[nn]), 2, 1, true, edge );
        int boundaryflag = Mesh->tag_get_integer(Mesh->BoundaryID,edge.front());
      
        MatrixXd normals(3,2);
        Mesh->tag_get_data(Mesh->BoundaryNormal,&(vertices[nn]),2,&normals(0)); 
        
        for( int i=0; i<nBasisE; i++ ) {
          for( int j=0; j<nBasisE; j++ ) { 
            if( j == 2 ) normal = (0.5*normals.col(0)+0.5*normals.col(1)).normalized();
            else normal = normals.col(j);
            MatrixXd pos0 = posE.col(j);
            Sys->setBCData(pos0,normal,boundaryflag);
            
            MatrixXd AmBC = 0.5*Am*Sys->XBC;
            for( int k=0; k<Sys->nEqn; k++ )
              for( int l=0; l<Sys->nEqn; l++ )
                Atotal.coeffRef(offset+k*nBasis+idx(nn,i),offset+l*nBasis+idx(nn,j)) += AmBC(k,l)*AssemMassE(i,j);	     
              
            MatrixXd AmBCrhs = 0.5*Am*Sys->BCrhs;	    
            for( int k=0; k<Sys->nEqn; k++)
              rhs(offset+k*nBasis+idx(nn,i)) += -AssemMassE(i,j)*AmBCrhs(k);	    
          };  
        };	
      };      
    };
    
  };
  stopTimer();
  
  startTimer(" compress sparse matrix\n"); 
  Atotal.prune(1,1e-13);  
  Atotal.makeCompressed();  
  stopTimer();
  
  nnz = Atotal.nonZeros();
  printf( "Non zeros: %ld\n", nnz );
  printf( "Non zeros per row: %f\n", 1.0*nnz/nTot );
  
};

void printMatrix( const MatrixXd& Tmp, int i0, int j0, int i1, int j1 )
{
  printf( "\n {\n" );
  for( int i=i0; i<i1; i++ ) {
    if( i>i0 ) printf( ",\n");
    printf( " { " );
    for( int j=j0;j<j1; j++ ) {
      if( j>j0 ) printf(", ");
      printf( "%20.12f ", Tmp.coeffRef(i,j) );
    };
    printf( "}" );
  };
  printf( " }\n\n" );
  cout << Tmp.block(i0,j0,i1-i0,j1-j0);
  exit(0);
};  

void Numerics::LinearSolve()
{
  MKL_INT iparm[64];
  void *pt[64];
  for(int i=0; i<64; i++ ) pt[i]=0;
  
/*
  int  n, *ia, *ja, nz;
  double *A, *b, *x;
  n = Atotal.rows();
  nz = Atotal.nonZeros();

//  ia = (long long *) malloc( (n+1)*sizeof(long long) );
//  ja = (long long *) malloc(   nz *sizeof(long long) );
//  A  = (double *) malloc(   nz *sizeof(double) );
//  b  = (double *) malloc(   n *sizeof(double) );

  ia = (int *) malloc( (n+1)*sizeof(int) );
  ja = (int *) malloc(   nz *sizeof(int) );
  A  = (double *) malloc(   nz *sizeof(double) );
  b  = (double *) malloc(   n *sizeof(double) );
  for(int i=0; i<n+1; i++ )
	ia[i] =  Atotal.outerIndexPtr()[i];
  for(int i=0; i<nz; i++ )
	ja[i] =  Atotal.innerIndexPtr()[i];
  for(int i=0; i<nz; i++ )
	A[i] =  Atotal.valuePtr()[i];
  for(int i=0; i<n; i++ )
	b[i] =  rhs(i); 

  for(int i=0; i<n+1; i++){
     if(ia[i] < 0|| ia[i]!=ia[i]){
	cout<<"ia["<<i<<"] = "<< ia[i]<< endl;
	break;
	}
  }
  for(int i=0; i<nz; i++){
     if(ja[i] < 0 || ja[i]!=ja[i]){
	cout<<"ja["<<i<<"] = "<< ja[i]<< endl;
	break;
	}
  }
// there is nan in A and b

  for(int i=0; i<nz; i++){
     if( A[i] != A[i] ){
	cout<<"A["<<i<<"] = "<< A[i]<< endl;
	break;
	}
  }
  for(int i=0; i<n; i++){
     if( b[i] != b[i] ){
	cout<<"b["<<i<<"] = "<< b[i]<< endl;
	break;
	}
  }
  for(int i=0; i<n; i++){
     if( rhs(i) != rhs(i) ){
	cout<<"rhs["<<i<<"] = "<<rhs(i)<< endl;
	break;
	}
  }
*/

  //printMatrix( Atotal, 10*Sys->nEqn*nBasis, 10*Sys->nEqn*nBasis, 11*Sys->nEqn*nBasis, 11*Sys->nEqn*nBasis ); 
  
  if(which ==1){

  startTimer(">> LU decomposition\n"); 

  PardisoSolve(11, Atotal.rows(), (MKL_INT *)Atotal.outerIndexPtr(), (MKL_INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
	       Atotal.nonZeros(), &rhs(0), &sol(0), pt, -11, iparm);
  stopTimer(); 
  startTimer("");
  PardisoSolve(11, Atotal.rows(), (MKL_INT *)Atotal.outerIndexPtr(), (MKL_INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
	       Atotal.nonZeros(), &rhs(0), &sol(0), pt, 11, iparm);
  stopTimer(); 
  startTimer("");
  PardisoSolve(11, Atotal.rows(), (MKL_INT *)Atotal.outerIndexPtr(), (MKL_INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
	       Atotal.nonZeros(), &rhs(0), &sol(0), pt, 22, iparm);  
  stopTimer(); 

  startTimer(" Solve\n");  
  PardisoSolve(11, Atotal.rows(), (MKL_INT *)Atotal.outerIndexPtr(), (MKL_INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
	       Atotal.nonZeros(), &rhs(0), &sol(0), pt, 33, iparm);    
  stopTimer(); 

  // free memory
  PardisoSolve(11, Atotal.rows(), (MKL_INT *)Atotal.outerIndexPtr(), (MKL_INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
  	       Atotal.nonZeros(), &rhs(0), &sol(0), pt, -1, iparm);    
  }
  else
       PardisoSolve_original(11, Atotal.rows(), (MKL_INT *)Atotal.outerIndexPtr(), (MKL_INT *)Atotal.innerIndexPtr(), Atotal.valuePtr(), 
	       Atotal.nonZeros(), &rhs(0), &sol(0), pt, 11, iparm);
  startTimer(" Compute residuum\n");  
  rhs = Atotal*sol-rhs;
  res = rhs.norm();
  stopTimer();
  
  cout << "Residuum: " << res << endl;
};



void Numerics::PutSolution()
{
  Mesh->tag_get_handle("Temperature", 1, MB_TYPE_DOUBLE, Temperature, MB_TAG_DENSE | MB_TAG_CREAT );
  Mesh->tag_get_handle("Velocity", 2, MB_TYPE_DOUBLE, Velocity, MB_TAG_DENSE | MB_TAG_CREAT );
  for (Range::iterator vert = Mesh->Vertices.begin(); vert != Mesh->Vertices.end(); vert++) {
    int idx = Mesh->get_ID(*vert);
    Range Elems;
    double theta = 0.0;
    MatrixXd s;
    s = MatrixXd::Zero(2,1);
    Mesh->get_adjacencies(&(*vert), 1, 2, true, Elems );
    for (Range::iterator elem = Elems.begin(); elem != Elems.end(); elem++) {
      const EntityHandle *vertices;
      int type;	
      Mesh->get_connectivity(*elem, vertices, type, true);
      VectorXi ID(type);
      Mesh->tag_get_data(Mesh->ID,vertices,type,&ID(0));
      int pp = 0; while( pp < 3 && ID[pp] != idx) { pp++; }
      int offset =  Mesh->get_ID(*elem)*Sys->nEqn*nBasis;
      
      theta += sol(offset+3*nBasis+pp);
      s(0) += sol(offset+1*nBasis+pp);
      s(1) += sol(offset+2*nBasis+pp);
    };
    theta = theta / Elems.size();
    s = s / Elems.size();
      
    Mesh->tag_set_data(Temperature,&(*vert),1,&theta);      
    Mesh->tag_set_data(Velocity,&(*vert),1,&s(0));      
  };
};

char *Numerics::Convert( double a, char *s )
{
  int i = fabs(a+1e-7);
  int j = fabs(10*(a-i+1e-7));
  int k = fabs(100*(a-i-0.1*j+1e-7));
  int l = fabs(1000*(a-i-0.1*j-0.01*k+1e-7));
  int m = fabs(10000*(a-i-0.1*j-0.01*k-0.001*l+1e-7));

  if( (m != 0)&&(j == 0)&&(k == 0) ) sprintf( s, "%dp%d%d%d%d", i, j, k, l, m );
  else if( l != 0 ) sprintf( s, "%dp%d%d%d", i, j, k, l );
  else if( k != 0 ) sprintf( s, "%dp%d%d", i, j, k );
  else sprintf( s, "%dp%d", i, j );
  return( s );
};

void Numerics::OutputSimple()
{
  char s1[20], s2[20], s3[20], s4[20], name[300], solv[300];
  FILE *fp;

  //sprintf( solv, "HybridFDFEM_%s", Type );   
  //sprintf( name, "%s_nst%d_Kn%s_Om%s_t%s.dat", solv, 
	//   FullGridSize( B, bMax ), Convert(Kn,s1), Convert(Om,s3), s2 );
 
  sprintf( name, "../output/output_nst%d_h%s.dat", Mesh->Vertices.size(), Convert(Mesh->hMax,s1) );

  fp = fopen( name, "w" );
  
  Tag Temperature; 
  Tag Velocity; 
  Mesh->tag_get_handle("Temperature", Temperature );
  Mesh->tag_get_handle("Velocity", Velocity );

  for (Range::iterator vert = Mesh->Vertices.begin(); vert != Mesh->Vertices.end(); vert++) {
    MatrixXd pos(3,1), temp(1,1), velocity(2,1);
    Mesh->get_coords(&(*vert),1,&pos(0));
    Mesh->tag_get_data(Temperature,&(*vert),1,&temp(0));      
    Mesh->tag_get_data(Velocity,&(*vert),1,&velocity(0));
    
    fprintf( fp, "%20.8f", pos(0) );
    fprintf( fp, "%20.8f", pos(1) );
    fprintf( fp, "%20.8f", temp(0) );
    fprintf( fp, "%20.8f", velocity(0) );
    fprintf( fp, "%20.8f", velocity(1) );
    fprintf( fp, "\n" );          
  };
  
  fclose( fp );
  printf( ">> File %s writing done\n", name );

};




