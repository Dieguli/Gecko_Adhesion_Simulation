//---------------------------------------------------------------------------

#ifndef TetraH
#define TetraH

#include "Matrix3d.h"
#include "Vertex.h"
#include "Bounds3d.h"

//---------------------------------------------------------------------------

class Triangle;
class TetraMesh;

//--------------------------------------------------------------------
class Tetra {
//--------------------------------------------------------------------
public:
  Tetra(TetraMesh *mesh, Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4);
  ~Tetra();

  Vertex* sideVertex(int side, int nr);
  Vertex* oppositeVertex(int side);
  int     sideNr(Vertex *v1, Vertex *v2, Vertex *v3);
  bool    contains(Vertex *v);
  void    replaceVertex(Vertex *v, Vertex *newV);
  void    computeStress(float &eigenVal, Vector3d &eigenVec);

  float   height(int side, bool originalCoords);
  void    barycentricCoords(bool originalCoords, Vector3d &pos,
            float &b0, float &b1, float &b2, float &b3);

  void    getCenter(Vector3d &center);

  void    draw();
  void    updateBounds();
  void    getBounds(bool originalCoords, Bounds3d &bounds);
  bool    pointInside(Vector3d &p, bool originalCoords);

  void    updateOrientation();
  void    resetOrientation() { _orientation.setUnity(); }
  void    getOrientation(Matrix3d &o) { o = _orientation; }
  void    setOrientation(Matrix3d &o) { _orientation = o; }

  void    addElasticCoeffs();

  void    addNeighborReferences();
  void    subNeighborReferences();

  // representation
  Vertex*   _vertex[4];
  Tetra*    _neighbor[4];
  Triangle* _triangle[4];

  Bounds3d  _bounds;
  TetraMesh *_itsMesh;
  int      _mark;
  int      _fractureSide;      // -1 or +1
  bool     _deleted;
  float    _mass;

private:
  float worldVolume();
  float originalVolume();
  void  getStrainStressRelation(float &a, float &b, float &c);
  void  computeN(Vector3d N[], float &det);
  void  computeStiffness();
  void  computeStrain(Matrix3d &strain);
  void  computeStress(Matrix3d &stress);
  void  elasticForce(int vertNr, Vector3d &f);
  void  computeTransformationMatrix(Matrix3d &A);
  void  computeBasis(bool originalCoords, Matrix3d &b);

  // FEM data -------------------

  Matrix3d  _stiffness[4][4];

  Matrix3d  _orientation;
};


#endif
