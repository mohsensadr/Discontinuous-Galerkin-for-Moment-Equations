
class Numerics  
{
  public:
    Numerics( MeshCore *Mesh, System *Sys, int Order );
    
    void Assemblation();
    void LinearSolve();
    void PutSolution();
    void OutputSimple();
    void ComputeError( char *Text );
    void InfoText( char *Text );
    
    MeshCore *Mesh;
    System *Sys;
    Tag Temperature, Velocity; 
    int nBasis, order, nTot, nnz;
    double res;
    
  private:
    char *Convert( double a, char *s );
    MatrixXd AssembleDiffX( const MatrixXd& pos, int order, int type );
    MatrixXd AssembleDiffY( const MatrixXd& pos, int order, int type );
    MatrixXd AssembleMass( const MatrixXd& pos, int order, int type );
    MatrixXd AssembleMassE( const MatrixXd& pos, int order );
    void LagrangePts( MatrixXd& Pts, int order, int type );
    MatrixXd cross( const MatrixXd& vector );
    
    SpMatrixL Atotal;
    MatrixXd sol, rhs;
};
