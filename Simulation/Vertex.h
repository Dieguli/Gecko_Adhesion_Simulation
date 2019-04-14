//---------------------------------------------------------------------------

#ifndef VertexH
#define VertexH
//---------------------------------------------------------------------------

#define MAX_NEIGHBORS 100

#include "Vector3d.h"
#include "Matrix3d.h"

class Vertex;
class TetraMesh;
class Tetra;

//-----------------------------------------------
class VertexNeighbor
//-----------------------------------------------
{
public:
  Vertex*  _vertex;
  Matrix3d _stiffness;
  Matrix3d _equation_A;	// linear equation Ax = b
  int      _refCount;
};

//-----------------------------------------------
class Vertex
//-----------------------------------------------
{
public:
  Vertex(float x, float y, float z);
  ~Vertex();

  // public methods
  void draw(float size);

  void transferInfo(Vertex *v); // for vertex split
  void transferInfo(Vertex *v1, Vertex *v2, float f); // for edge split

  void addNeighborReference(Vertex *v);
  void subNeighborReference(Vertex *v);
  void removeNeighbors();

  void addStiffness(Vertex *v, Matrix3d &m);
  void subStiffness();

  void setupDynamicEquation(float delta_t, float damping);
  void setupStaticEquation();
  void multiply(bool with_p);  // conjugate gradients
  void gaussSeidel();          // gauss seidel

  bool preventOverstreching();
  void elasticForce(Vector3d &f);

  int  numNeighbors() { return _numNeighbors; }

  // representation
  bool      _animated;
  bool      _frontera;
  Vector3d  _originalCoord;
  Vector3d  _worldCoord;
  Vector3d  _velocity;
  float     _mass;

  Vector3d  _normal;            // for smooth normals

  float     _forceDistr;        // grab force distribution
  Vector3d  _fExt;              // external force acting on this vertex
  Vector3d  _fColl;             // stored separately for derivative computation
  Vector3d  _f0;                // elastic force at rest state (plastic mats)

  Vector3d  _solution;          // used in conjugate gradients
  Vector3d  _u, _p, _r;         // "
  Vector3d  _equation_b;        // linear equation Ax = b

  int       _mark;
  int       _fractureMark;
  bool      _fractureTip;

  Vector3d  _collDisp;          // collision information
  bool      _colliding;
  TetraMesh*_collMesh;
  Tetra*    _collTetra;

private:
  int             _numNeighbors;
  VertexNeighbor* _neighbor[MAX_NEIGHBORS];
};

#endif
