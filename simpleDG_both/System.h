  

class System 
{
  public:
    System( double Tau );
    int nEqn, nBC;
    SpMatrix Ax, Ay, P, BC;
    VectorXi OddVar;
    MatrixXd XBC, BCrhs;
    
    MatrixXd Force( MatrixXd& pos );
    void setBCData( MatrixXd& pos, MatrixXd& normal, int boundaryID );
    
    MatrixXd Aminus( const MatrixXd& normal );
    
    double A0, A1, A2;
    double tau;
    double chi, theta0, theta1, uW, v0, eps;
    std::string SystemText;
    
    void InfoText( char *Text );
    
  private:
    MatrixXd mirrow( const MatrixXd& vector );
    
    MatrixXd X0, Aminus1D;
    SpMatrix Projector( const MatrixXd& normal );
    SpMatrix invProjector( const MatrixXd& normal );
    
};

