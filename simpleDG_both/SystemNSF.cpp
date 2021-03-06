#include <stdio.h>
#include <iostream>
//#include "mkl.h"

#include <EigenSetup.h>
#include <Eigen/Eigenvalues>
#include <math.h>

#include "System.h"

#include "Tools.h"


// SYSTEM CLASS
//----------------------------------------------

void System::InfoText( char *text )
{
  sprintf( text, "System = %s, nEqn = %d, tau = %.2f, chi = %.2f, theta0 = %.2f, theta1 = %.2f, v0 = %.2f, A0 = %.2f, A1 = %.2f, A2 = %.2f, eps = %.2e",
	   SystemText.c_str(), nEqn, tau, chi, theta0, theta1, v0, A0, A1, A2, eps );
  
};


System::System( double Tau )
{
  startTimer(">> System Setup \n"); 
  VectorXi ColSize;
  
  nEqn = 9;
  nBC = 3;
  SystemText = "NSF(heat+flow/BCodd)"; 
  
  tau = Tau;
  chi = 1.0;
  theta0 = 1.0;
  theta1 = 2.0;
  v0 = 1.0;
  
  A0 = 2.0*tau;
  A1 = 0.0;
  A2 = 0.0;
  
  eps = 1e-5;
    
  int Annz = 11;
  VectorXi Ai(Annz), Aj(Annz);
  MatrixXd Aval(Annz,1);
  Ai << 0,1,1,1,2,3,3,4,5,6,7;
  Aj << 1,4,3,0,5,7,1,1,2,1,3;
  Aval << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0/3, 1.0/2, -1.0/3, 1.0;
  
  int BCnnz = 7;
  VectorXi BCi(BCnnz), BCj(BCnnz);
  MatrixXd BCval(BCnnz,1);
  BCi << 0,0,0,1,1,2,2;
  BCj << 0,1,4,2,5,3,7;
  BCval << 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0;

  OddVar.resize(nBC);
  OddVar << 1,5,7; // caution: C-indices

  //***************************************//
  
  Ax.resize(nEqn,nEqn);
  ColSize.resize(nEqn);
  ColSize.setZero();
  for( int ix=0; ix < Annz; ix++ ) ColSize(Ai(ix))++; 
  Ax.reserve(ColSize);
  for( int ix=0; ix < Annz; ix++ ) Ax.coeffRef(Ai(ix),Aj(ix)) = Aval(ix,0);
  Ax.prune(1,1e-13);
  
  Ay.resize(nEqn,nEqn);
  Matrix<double,2,1> normal(0.0,1.0);
  Ay = invProjector(normal) * Ax * Projector(normal);
  Ay.prune(1,1e-13);

  P.resize(nEqn,nEqn);
  ColSize << 0,0,0,0,1,1,1,1,1;
  P.reserve(ColSize);
  P.coeffRef(4,4) = 1.0/tau;
  P.coeffRef(5,5) = 1.0/tau;
  P.coeffRef(6,6) = 1.0/tau;
  P.coeffRef(7,7) = 1.0/tau;
  P.coeffRef(8,8) = 1.0/tau;
  P.prune(1,1e-13);  
  
  BC.resize(nBC,nEqn);
  ColSize.resize(nBC);
  ColSize.setZero();
  for( int ix=0; ix < BCnnz; ix++ ) ColSize(BCi(ix))++; 
  BC.reserve(ColSize);
  for( int ix=0; ix < BCnnz; ix++ ) BC.coeffRef(BCi(ix),BCj(ix)) = BCval(ix,0);
  BC.prune(1,1e-13);

  stopTimer();
  
  startTimer(" characteristic setup \n"); 
  EigenSolver<MatrixXd> ES(Ax);
  MatrixXd vecs = ES.pseudoEigenvectors();
  VectorXd vals = ES.pseudoEigenvalueMatrix().diagonal();
  double maxEV = vals.cwiseAbs().maxCoeff();
  Aminus1D = vecs*(vals.cwiseAbs()-vals).asDiagonal()*vecs.inverse();

  X0.resize(nEqn,nBC);
  for( int i=0; i<nBC; i++ ) X0.col(i) = VectorXd::Unit(nEqn,OddVar(i)); // odd BCs
  
  //Ordering Order(vals);
  //for( int i=0; i<nBC; i++ ) X0.col(i) = vecs.col(Order.index(i)); // Riemann BCs
  
  XBC.resize(nEqn,nEqn);
  BCrhs.resize(nEqn,1); 
  stopTimer();

  cout << "Number of equations " << nEqn << "\n";
  cout << "Number of boundary conditions " << nBC << "\n";
  cout << "maximal charac. velocity " << maxEV << "\n";
  
};

MatrixXd System::Force( MatrixXd& pos )
{
  MatrixXd res(nEqn,1); 
  res = MatrixXd::Zero(nEqn,1);
  double r = pos.norm();
  
  res(3) = A0 + A2*r*r + A1*pos(0)/r;
  return( res );  
};

void System::setBCData( MatrixXd& pos, MatrixXd& normal, int boundaryID )
{
  MatrixXd g(nBC,1);
  double norm = pos.norm(), chiW, eps_p, eps_n, vtW, vnW, pW, thetaW;
  
  switch(boundaryID) {
  case 10:// inner circle
    chiW = chi;
    //eps_n = eps;
    eps_p = 1.0/eps;
    vnW = 0.0;
    vtW = 0.0;
    pW = 0.0;
    thetaW = theta0;
    break;
  case 20:// outer circle
    chiW = chi;
    //eps_n = 1.0;
    eps_p = eps;
    vnW = v0*normal(0);
    vtW = -v0*normal(1);
    pW = -0.23*normal(0);
    thetaW = theta1;
    break;
  }; 
 
  g << vnW*eps_p - pW*chiW, -vtW*chiW, -thetaW*chiW;
  
  XBC = chiW*BC;
  for( int ix=0; ix<OddVar.size(); ix++ ) XBC.coeffRef(ix,OddVar(ix)) = -1.0;
  XBC.coeffRef(0,1) = -eps_p;
  //XBC.coeffRef(0,0) = eps_n;
  //XBC.coeffRef(0,4) = eps_n;

  BCrhs = invProjector(normal)*X0*(XBC*X0).inverse()*g;
  XBC = invProjector(normal)*X0*(XBC*X0).inverse()*XBC*Projector(normal);
};

static inline void SpBlock( int idx, const MatrixXd& P, SpMatrix& T )
{
  for( int i=0; i < P.rows(); i++ )
    for( int j=0; j < P.cols(); j++ )
      T.coeffRef(idx+i,idx+j) = P(i,j); 
};
 
SpMatrix System::Projector( const MatrixXd& normal )
{
  double nx = normal(0), ny = normal(1);
  double nxnx = nx*nx, nyny = ny*ny;
  int idx;
  MatrixXd P0(1,1), P1(2,2), P2(3,3), P3(4,4);
  P0 << 1.0;
  P1 << nx, ny, -ny, nx;
  P2 << nxnx, 2*nx*ny, nyny, -nx*ny, nxnx-nyny, nx*ny, nyny, -2*nx*ny, nxnx;
  P3 << nx*nxnx, 3*ny*nxnx, 3*nx*nyny, ny*nyny, -ny*nxnx, nx*nxnx - 2*nx*nyny, 2*ny*nxnx - ny*nyny, nx*nyny, 
        nx*nyny, -2*ny*nxnx + ny*nyny, nx*nxnx - 2*nx*nyny, ny*nxnx, -ny*nyny, 3*nx*nyny, -3*ny*nxnx, nx*nxnx;
  
  SpMatrix T(nEqn,nEqn);
  VectorXi ColSize(nEqn);
  ColSize << 1,2,2,1,3,3,3,2,2;

  SpBlock( 0, P0, T );
  SpBlock( 1, P1, T );
  SpBlock( 3, P0, T );
  SpBlock( 4, P2, T );
  SpBlock( 7, P1, T );
  T.makeCompressed();
  return( T );
};

MatrixXd System::mirrow( const MatrixXd& vector )
{
  MatrixXd mirrow(3,1);
  mirrow << vector(0,0),-vector(1,0),0;
  return( mirrow );
};

SpMatrix System::invProjector( const MatrixXd& normal )
{
  return( Projector( mirrow(normal) ) );
};

//------------------------------------------------
// Upwinding Flux Matrix
//------------------------------------------------
MatrixXd System::Aminus( const MatrixXd& normal )
{
  return(invProjector(normal)*Aminus1D*Projector(normal)  );
};


