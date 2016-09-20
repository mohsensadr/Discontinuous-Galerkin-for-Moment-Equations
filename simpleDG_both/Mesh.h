#include "moab/Core.hpp"
#include "moab/Range.hpp"
using namespace moab;

class MeshCore : public Core  
{
  public:
    MeshCore();
    MeshCore( const char *filename );
    void setup( const char *filename );
    int get_ID( const EntityHandle& entity_handle );
    int tag_get_integer( const Tag tag_handle, const EntityHandle& entity_handle );
    void InfoText( char *Text );
    
    double hMax;
    int Type;
    Range Vertices, Edges, Elements, BCVertices;
    Tag ID, BoundaryTag, BoundaryNormal, BoundaryID;
    char MeshFile[100];
    
  private:
    MatrixXd NormalApprox( const MatrixXd& P0, const MatrixXd& Pts );
};
