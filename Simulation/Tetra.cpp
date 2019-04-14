//---------------------------------------------------------------------------

#include "Draw.h"

#include <math.h>
#include "Tetra.h"
#include "Vertex.h"
#include "PhysicsParams.h"
#include <stdio.h>

#define COLOR_RED 40
#define COLOR_GREEN 160
#define COLOR_BLUE 50

//-----------------------------------------------------------------------
Tetra::Tetra(TetraMesh *mesh, Vertex *v1, Vertex *v2, Vertex *v3, Vertex *v4)
//-----------------------------------------------------------------------
{
  _vertex[0] = v1;
  _vertex[1] = v2;
  _vertex[2] = v3;
  _vertex[3] = v4;

  _itsMesh = mesh;
  _deleted = false;

  _mark = -1;
  _fractureSide = -1;
  _orientation.setUnity();

  computeStiffness();
  addNeighborReferences();
  updateBounds();
}


//--------------------------------------------------------------------
Tetra::~Tetra()
//--------------------------------------------------------------------
{
}


//--------------------------------------------------------------------
Vertex* Tetra::sideVertex(int side, int nr)
//--------------------------------------------------------------------
{
  switch (side) {
    case 0:
      switch(nr) {
        case 0 : return _vertex[0];
        case 1 : return _vertex[3];
        case 2 : return _vertex[1];
      }; break;
    case 1:
      switch(nr) {
        case 0 : return _vertex[1];
        case 1 : return _vertex[3];
        case 2 : return _vertex[2];
      }; break;
    case 2:
      switch(nr) {
        case 0 : return _vertex[0];
        case 1 : return _vertex[2];
        case 2 : return _vertex[3];
      }; break;
    case 3:
      switch(nr) {
        case 0 : return _vertex[0];
        case 1 : return _vertex[1];
        case 2 : return _vertex[2];
    }; break;
  }
  return 0;
}


//--------------------------------------------------------------------
Vertex* Tetra::oppositeVertex(int side)
//--------------------------------------------------------------------
{
  switch (side) {
    case 0: return _vertex[2];
    case 1: return _vertex[0];
    case 2: return _vertex[1];
    case 3: return _vertex[3];
  }
  return 0;
}


//--------------------------------------------------------------------
int Tetra::sideNr(Vertex *v1, Vertex *v2, Vertex *v3)
//--------------------------------------------------------------------
{
  int i, sideCount[4];

  for (i = 0; i < 4; i++)
    sideCount[i] = 0;

  if (v1 == _vertex[0] || v2 == _vertex[0] || v3 == _vertex[0]) {
    sideCount[0]++; sideCount[2]++; sideCount[3]++;
  }
  if (v1 == _vertex[1] || v2 == _vertex[1] || v3 == _vertex[1]) {
    sideCount[0]++; sideCount[1]++; sideCount[3]++;
  }
  if (v1 == _vertex[2] || v2 == _vertex[2] || v3 == _vertex[2]) {
    sideCount[1]++; sideCount[2]++; sideCount[3]++;
  }
  if (v1 == _vertex[3] || v2 == _vertex[3] || v3 == _vertex[3]) {
    sideCount[0]++; sideCount[1]++; sideCount[2]++;
  }
  for (i = 0; i < 4; i++)
    if (sideCount[i] == 3) return i;

  return -1;
}


//--------------------------------------------------------------------
bool Tetra::contains(Vertex *v)
//--------------------------------------------------------------------
{
  for (int i = 0; i < 4; i++)
    if (_vertex[i] == v) return true;
  return false;
}

//--------------------------------------------------------------------
void Tetra::replaceVertex(Vertex *v, Vertex *newV)
//--------------------------------------------------------------------
{
  subNeighborReferences();

  for (int i = 0; i < 4; i++) {
    if (_vertex[i] == v)
      _vertex[i] = newV;
  }
  if (v->_originalCoord.dist(newV->_originalCoord) > 1e-8) {
    computeStiffness();
  }
  addNeighborReferences();
}

//--------------------------------------------------------------------
void Tetra::updateBounds()
//--------------------------------------------------------------------
{
  _bounds.clear();
  for (int i = 0; i < 4; i++)
    _bounds.include(_vertex[i]->_worldCoord);
}


//--------------------------------------------------------------------
void Tetra::getBounds(bool originalCoords, Bounds3d &bounds)
//--------------------------------------------------------------------
{
  bounds.clear();

  if (originalCoords) {
    for (int i = 0; i < 4; i++)
      bounds.include(_vertex[i]->_originalCoord);
  }
  else {
    for (int i = 0; i < 4; i++)
      bounds.include(_vertex[i]->_worldCoord);
  }
}



//--------------------------------------------------------------------
float Tetra::worldVolume()
//--------------------------------------------------------------------
{
  float x0,x1,x2,x3, y0,y1,y2,y3, z0,z1,z2,z3;

  x0 = _vertex[0]->_worldCoord.x;
  y0 = _vertex[0]->_worldCoord.y;
  z0 = _vertex[0]->_worldCoord.z;
  x1 = _vertex[1]->_worldCoord.x - x0;
  y1 = _vertex[1]->_worldCoord.y - y0;
  z1 = _vertex[1]->_worldCoord.z - z0;
  x2 = _vertex[2]->_worldCoord.x - x0;
  y2 = _vertex[2]->_worldCoord.y - y0;
  z2 = _vertex[2]->_worldCoord.z - z0;
  x3 = _vertex[3]->_worldCoord.x - x0;
  y3 = _vertex[3]->_worldCoord.y - y0;
  z3 = _vertex[3]->_worldCoord.z - z0;

  float det =   x1*y2*z3 + x2*y3*z1 + x3*y1*z2
              - x1*y3*z2 - x2*y1*z3 - x3*y2*z1;
  return fabs(det / 6.0);
}

//--------------------------------------------------------------------
float Tetra::originalVolume()
//--------------------------------------------------------------------
{
  float x0,x1,x2,x3, y0,y1,y2,y3, z0,z1,z2,z3;

  x0 = _vertex[0]->_originalCoord.x;
  y0 = _vertex[0]->_originalCoord.y;
  z0 = _vertex[0]->_originalCoord.z;
  x1 = _vertex[1]->_originalCoord.x - x0;
  y1 = _vertex[1]->_originalCoord.y - y0;
  z1 = _vertex[1]->_originalCoord.z - z0;
  x2 = _vertex[2]->_originalCoord.x - x0;
  y2 = _vertex[2]->_originalCoord.y - y0;
  z2 = _vertex[2]->_originalCoord.z - z0;
  x3 = _vertex[3]->_originalCoord.x - x0;
  y3 = _vertex[3]->_originalCoord.y - y0;
  z3 = _vertex[3]->_originalCoord.z - z0;

  float det =   x1*y2*z3 + x2*y3*z1 + x3*y1*z2
              - x1*y3*z2 - x2*y1*z3 - x3*y2*z1;
  return fabs(det / 6.0);
}


//--------------------------------------------------------------------
void Tetra::getStrainStressRelation(float &a, float &b, float &c)
//--------------------------------------------------------------------
{
  // returns coeffs of 6 x 6 matrix E,  stress = E * strain
  //     a b b 0 0 0
  //     b a b 0 0 0
  // E = b b a 0 0 0
  //     0 0 0 c 0 0
  //     0 0 0 0 c 0
  //     0 0 0 0 0 c

  float d = physicsParams.youngModulus / (1.0 + physicsParams.poissonRatio) /
               (1.0 - 2.0*physicsParams.poissonRatio);

  a = (1.0 - physicsParams.poissonRatio) * d;
  b = physicsParams.poissonRatio * d;
  c = physicsParams.youngModulus / 2.0 / (1.0 + physicsParams.poissonRatio);
}


//--------------------------------------------------------------------
void Tetra::computeN(Vector3d N[], float &det)
//--------------------------------------------------------------------
{
  float x0,x1,x2,x3, y0,y1,y2,y3, z0,z1,z2,z3;

  x0 = _vertex[0]->_originalCoord.x;
  y0 = _vertex[0]->_originalCoord.y;
  z0 = _vertex[0]->_originalCoord.z;
  x1 = _vertex[1]->_originalCoord.x - x0;
  y1 = _vertex[1]->_originalCoord.y - y0;
  z1 = _vertex[1]->_originalCoord.z - z0;
  x2 = _vertex[2]->_originalCoord.x - x0;
  y2 = _vertex[2]->_originalCoord.y - y0;
  z2 = _vertex[2]->_originalCoord.z - z0;
  x3 = _vertex[3]->_originalCoord.x - x0;
  y3 = _vertex[3]->_originalCoord.y - y0;
  z3 = _vertex[3]->_originalCoord.z - z0;

  det =   x1*y2*z3 + x2*y3*z1 + x3*y1*z2
        - x1*y3*z2 - x2*y1*z3 - x3*y2*z1;
  float det1 = 1.0 / det;
  N[1].x = (z2*y3 - y2*z3) * det1;
  N[2].x = (y1*z3 - z1*y3) * det1;
  N[3].x = (z1*y2 - y1*z2) * det1;
  N[0].x = -N[1].x - N[2].x - N[3].x;

  N[1].y = (x2*z3 - z2*x3) * det1;
  N[2].y = (z1*x3 - x1*z3) * det1;
  N[3].y = (x1*z2 - z1*x2) * det1;
  N[0].y = -N[1].y - N[2].y - N[3].y;

  N[1].z = (y2*x3 - x2*y3) * det1;
  N[2].z = (x1*y3 - y1*x3) * det1;
  N[3].z = (y1*x2 - x1*y2) * det1;
  N[0].z = -N[1].z - N[2].z - N[3].z;
}



//--------------------------------------------------------------------
void Tetra::computeStiffness()
//--------------------------------------------------------------------
{
  // K = vol * B^T*E*B

  float a,b,c;
  getStrainStressRelation(a,b,c);

  Vector3d N[4];
  float det;
  computeN(N, det);
  float vol = fabs(det / 6.0);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      Matrix3d &Kij = _stiffness[i][j];
      Vector3d &Ni = N[i];
      Vector3d &Nj = N[j];
      Kij.a00 = a*Ni.x*Nj.x + c*(Ni.y*Nj.y + Ni.z*Nj.z);
      Kij.a01 = b*Ni.x*Nj.y + c*(Ni.y*Nj.x);
      Kij.a02 = b*Ni.x*Nj.z + c*(Ni.z*Nj.x);

      Kij.a10 = b*Ni.y*Nj.x + c*(Ni.x*Nj.y);
      Kij.a11 = a*Ni.y*Nj.y + c*(Ni.x*Nj.x + Ni.z*Nj.z);
      Kij.a12 = b*Ni.y*Nj.z + c*(Ni.z*Nj.y);

      Kij.a20 = b*Ni.z*Nj.x + c*(Ni.x*Nj.z);
      Kij.a21 = b*Ni.z*Nj.y + c*(Ni.y*Nj.z);
      Kij.a22 = a*Ni.z*Nj.z + c*(Ni.y*Nj.y + Ni.x*Nj.x);
      Kij.scale(vol);
    }
  }
}


//--------------------------------------------------------------------
void Tetra::computeStrain(Matrix3d &strain)
//--------------------------------------------------------------------
{
  // e = B * u

  Vector3d N[4];
  float det;
  computeN(N, det);

  strain.setZero();
  for (int i = 0; i < 4; i++) {
    Vertex *v = _vertex[i];
    Vector3d u;
    u.multTranspose(_orientation, v->_worldCoord);
    u.sub(v->_originalCoord);
    Vector3d &Ni = N[i];
    strain.a00 += Ni.x * u.x;
    strain.a11 += Ni.y * u.y;
    strain.a22 += Ni.z * u.z;
    strain.a01 += Ni.y * u.x + Ni.x * u.y;
    strain.a12 += Ni.z * u.y + Ni.y * u.z;
    strain.a20 += Ni.z * u.x + Ni.x * u.z;
  }
  strain.a02 = strain.a20;
  strain.a21 = strain.a12;
  strain.a10 = strain.a01;
}


//--------------------------------------------------------------------
void Tetra::computeStress(Matrix3d &stress)
//--------------------------------------------------------------------
{
  // s = E * e

  float a,b,c;
  getStrainStressRelation(a,b,c);

  Matrix3d strain; computeStrain(strain);
  strain.scale(-1.0);  // debug: why necessary?

  stress.a00 = a * strain.a00 + b * strain.a11 + b * strain.a22;
  stress.a11 = b * strain.a00 + a * strain.a11 + b * strain.a22;
  stress.a22 = b * strain.a00 + b * strain.a11 + a * strain.a22;
  stress.a01 = c * strain.a01;
  stress.a02 = c * strain.a02;
  stress.a12 = c * strain.a12;

  stress.a10 = stress.a01;
  stress.a20 = stress.a02;
  stress.a21 = stress.a12;
}


// ------------------------------------------------------------------------------
void Tetra::computeStress(float &eigenVal, Vector3d &eigenVec)
// ------------------------------------------------------------------------------
{
  Matrix3d stress;
  computeStress(stress);  
  stress.largestEigenvalue(eigenVal, eigenVec);
  eigenVec.mult(_orientation, eigenVec);
}


//--------------------------------------------------------------------
void Tetra::elasticForce(int vertNr, Vector3d &f)
//--------------------------------------------------------------------
{
  f.setZero();
  for (int j = 0; j < 4; j++) {
    Vertex *vj = _vertex[j];
    Vector3d p;
    p.multTranspose(_orientation, vj->_worldCoord);
    p.sub(vj->_originalCoord);
    p.mult(_stiffness[vertNr][j], p);
    f.add(p);
  }
  f.mult(_orientation, f);
}


//--------------------------------------------------------------------
void Tetra::addElasticCoeffs()
//--------------------------------------------------------------------
{
  for (int i = 0; i < 4; i++) {
    Vertex *vi = _vertex[i];
    Vector3d f, f0; f0.setZero();

    for (int j = 0; j < 4; j++) {
      Vertex *vj = _vertex[j];
      f.mult(_stiffness[i][j], vj->_originalCoord);
      f0.add(f);

      if (j >= i) {
        Matrix3d Ki;
        Ki.mult(_orientation, _stiffness[i][j]);
        Ki.multTranspose(_orientation);
        _vertex[i]->addStiffness(_vertex[j], Ki);
        if (j > i) {
          Ki.transpose();
          _vertex[j]->addStiffness(_vertex[i], Ki);
        }
      }
    }
    f0.mult(_orientation, f0);
    vi->_f0.sub(f0);
  }
}


//--------------------------------------------------------------------
void Tetra::addNeighborReferences()
//--------------------------------------------------------------------
{
  _mass = originalVolume() * physicsParams.density;
  for (int i = 0; i < 4; i++) {
    _vertex[i]->_mass += 0.25*_mass;
    for (int j = 0; j < 4; j++) {
      _vertex[i]->addNeighborReference(_vertex[j]);
    }
  }
}


//--------------------------------------------------------------------
void Tetra::subNeighborReferences()
//--------------------------------------------------------------------
{
  for (int i = 0; i < 4; i++) {
    _vertex[i]->_mass -= 0.25*_mass;
    for (int j = 0; j < 4; j++) {
      _vertex[i]->subNeighborReference(_vertex[j]);
    }
  }
}


// ------------------------------------------------------------------------
void Tetra::barycentricCoords(bool originalCoords, Vector3d &pos,
       float &b0, float &b1, float &b2, float &b3)
// ------------------------------------------------------------------------
{
  // to transform any point as the tetra compute
  // p' = _vertex[0]->_worldCoord + A * (p - _vertex[0]->_originalCoord)

  float x1,x2,x3, y1,y2,y3, z1,z2,z3;

  Vector3d v,v0;
  if (originalCoords) {
    v0 = _vertex[0]->_originalCoord;
    v = _vertex[1]->_originalCoord;
    x1 = v.x - v0.x; y1 = v.y - v0.y; z1 = v.z - v0.z;
    v = _vertex[2]->_originalCoord;
    x2 = v.x - v0.x; y2 = v.y - v0.y; z2 = v.z - v0.z;
    v = _vertex[3]->_originalCoord;
    x3 = v.x - v0.x; y3 = v.y - v0.y; z3 = v.z - v0.z;
  }
  else {
    v0 = _vertex[0]->_worldCoord;
    v = _vertex[1]->_worldCoord;
    x1 = v.x - v0.x; y1 = v.y - v0.y; z1 = v.z - v0.z;
    v = _vertex[2]->_worldCoord;
    x2 = v.x - v0.x; y2 = v.y - v0.y; z2 = v.z - v0.z;
    v = _vertex[3]->_worldCoord;
    x3 = v.x - v0.x; y3 = v.y - v0.y; z3 = v.z - v0.z;
  }

  float  det =   x1*y2*z3 + x2*y3*z1 + x3*y1*z2
               - x1*y3*z2 - x2*y1*z3 - x3*y2*z1;
  float det1 = 1.0 / det;

  float q10 = (y2*z3 - y3*z2) * det1;
  float q11 = (x3*z2 - x2*z3) * det1;
  float q12 = (x2*y3 - x3*y2) * det1;
  float q20 = (z1*y3 - y1*z3) * det1;
  float q21 = (x1*z3 - z1*x3) * det1;
  float q22 = (y1*x3 - x1*y3) * det1;
  float q30 = (y1*z2 - z1*y2) * det1;
  float q31 = (z1*x2 - x1*z2) * det1;
  float q32 = (x1*y2 - y1*x2) * det1;

  Vector3d p = pos; p.sub(v0);
  b1 = q10 * p.x + q11 * p.y + q12 * p.z;
  b2 = q20 * p.x + q21 * p.y + q22 * p.z;
  b3 = q30 * p.x + q31 * p.y + q32 * p.z;
  b0 = 1.0 - b1 - b2 - b3;
}


// ------------------------------------------------------------------------
void Tetra::computeTransformationMatrix(Matrix3d &A)
// ------------------------------------------------------------------------
{
  // to transform any point as the tetra compute
  // p' = _vertex[0]->_worldCoord + A * (p - _vertex[0]->_originalCoord)

  float x1,x2,x3, y1,y2,y3, z1,z2,z3;

  Vector3d v,v0;
  v0 = _vertex[0]->_originalCoord;
  v = _vertex[1]->_originalCoord;
  x1 = v.x - v0.x;
  y1 = v.y - v0.y;
  z1 = v.z - v0.z;
  v = _vertex[2]->_originalCoord;
  x2 = v.x - v0.x;
  y2 = v.y - v0.y;
  z2 = v.z - v0.z;
  v = _vertex[3]->_originalCoord;
  x3 = v.x - v0.x;
  y3 = v.y - v0.y;
  z3 = v.z - v0.z;

  float  det =   x1*y2*z3 + x2*y3*z1 + x3*y1*z2
               - x1*y3*z2 - x2*y1*z3 - x3*y2*z1;
  float det1 = 1.0 / det;

  float q10 = (y2*z3 - y3*z2) * det1;
  float q11 = (x3*z2 - x2*z3) * det1;
  float q12 = (x2*y3 - x3*y2) * det1;
  float q20 = (z1*y3 - y1*z3) * det1;
  float q21 = (x1*z3 - z1*x3) * det1;
  float q22 = (y1*x3 - x1*y3) * det1;
  float q30 = (y1*z2 - z1*y2) * det1;
  float q31 = (z1*x2 - x1*z2) * det1;
  float q32 = (x1*y2 - y1*x2) * det1;

  v0 = _vertex[0]->_worldCoord;
  v = _vertex[1]->_worldCoord;
  x1 = v.x - v0.x;
  y1 = v.y - v0.y;
  z1 = v.z - v0.z;
  v = _vertex[2]->_worldCoord;
  x2 = v.x - v0.x;
  y2 = v.y - v0.y;
  z2 = v.z - v0.z;
  v = _vertex[3]->_worldCoord;
  x3 = v.x - v0.x;
  y3 = v.y - v0.y;
  z3 = v.z - v0.z;

  A.a00 = x1 * q10 + x2 * q20 + x3 * q30;
  A.a01 = x1 * q11 + x2 * q21 + x3 * q31;
  A.a02 = x1 * q12 + x2 * q22 + x3 * q32;

  A.a10 = y1 * q10 + y2 * q20 + y3 * q30;
  A.a11 = y1 * q11 + y2 * q21 + y3 * q31;
  A.a12 = y1 * q12 + y2 * q22 + y3 * q32;

  A.a20 = z1 * q10 + z2 * q20 + z3 * q30;
  A.a21 = z1 * q11 + z2 * q21 + z3 * q31;
  A.a22 = z1 * q12 + z2 * q22 + z3 * q32;
}


// ------------------------------------------------------------------------
void Tetra::computeBasis(bool originalCoords, Matrix3d &b)
// ------------------------------------------------------------------------
{
  int i;
  Vector3d x[4];
  if (originalCoords)
    for (i = 0; i < 4; i++) x[i] = _vertex[i]->_originalCoord;
  else
    for (i = 0; i < 4; i++) x[i] = _vertex[i]->_worldCoord;

  Vector3d v, b1,b2,b3;

  b1.set(0,0,0);
  v.sub(x[1],x[0]); v.normalize(); b1.add(v);
  v.sub(x[2],x[0]); v.normalize(); b1.add(v);
  v.sub(x[3],x[0]); v.normalize(); b1.add(v);
  b1.normalize();

  b2.set(0,0,0);
  v.sub(x[0],x[1]); v.normalize(); b2.add(v);
  v.sub(x[2],x[1]); v.normalize(); b2.add(v);
  v.sub(x[3],x[1]); v.normalize(); b2.add(v);
  b2.normalize();

  b3.cross(b1, b2); b3.normalize();
  b2.cross(b1, b3); b2.scale(-1.0);

  b.setColumns(b1,b2,b3);
}


// ------------------------------------------------------------------------
void Tetra::updateOrientation()
// ------------------------------------------------------------------------
{
  Matrix3d A; computeTransformationMatrix(A);
  A.orthonormalize();
//  A.gramSchmidt();
  _orientation.set(A);

// another method to compute the rotation
/*
  Matrix3d b0,b1;
  computeBasis(true,  b0); b0.transpose();
  computeBasis(false, b1);
  _orientation.mult(b1, b0);
*/  
}


//--------------------------------------------------------------------
void Tetra::getCenter(Vector3d &center)
//--------------------------------------------------------------------
{
  center.setZero();
  for (int i = 0; i < 4; i++) center.add(_vertex[i]->_worldCoord);
  center.scale(0.25);
}


//--------------------------------------------------------------------
float Tetra::height(int side, bool originalCoords)
//--------------------------------------------------------------------
{
  Vector3d x0, x1, x2, a, n;
  if (originalCoords) {
    x0 = sideVertex(side, 0)->_originalCoord;
    x1 = sideVertex(side, 1)->_originalCoord;
    x2 = sideVertex(side, 2)->_originalCoord;
    a = oppositeVertex(side)->_originalCoord;
  }
  else {
    x0 = sideVertex(side, 0)->_worldCoord;
    x1 = sideVertex(side, 1)->_worldCoord;
    x2 = sideVertex(side, 2)->_worldCoord;
    a = oppositeVertex(side)->_worldCoord;
  }
  x1.sub(x0); x2.sub(x0);
  n.cross(x1, x2); n.normalize();
  float d = n.inner(x0);
  float an = a.inner(n);
  float nn = n.inner(n);
  if (nn == 0.0) return 0.0;
  return fabs((d - an) / nn);
}


/* ----------------------------------------------------------- */
bool Tetra::pointInside(Vector3d &p, bool originalCoords)
/* ----------------------------------------------------------- */
{
  float x1,x2,x3,x4, y1,y2,y3,y4, z1,z2,z3,z4, px, py, pz;
  float d0,d1,d2,d3,d4;

  int i;
  Vector3d x[4];
  if (originalCoords)
    for (i = 0; i < 4; i++) x[i] = _vertex[i]->_originalCoord;
  else
    for (i = 0; i < 4; i++) x[i] = _vertex[i]->_worldCoord;

  x4 = x[3].x; y4 = x[3].y; z4 = x[3].z;
  x1 = x[0].x - x4; y1 = x[0].y - y4; z1 = x[0].z - z4;
  x2 = x[1].x - x4; y2 = x[1].y - y4; z2 = x[1].z - z4;
  x3 = x[2].x - x4; y3 = x[2].y - y4; z3 = x[2].z - z4;

  px = p.x - x4; py = p.y - y4; pz = p.z - z4;

  d0 = x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1;
  d1 = px*y2*z3 - px*y3*z2 - x2*py*z3 + x2*y3*pz + x3*py*z2 - x3*y2*pz;
  d2 = x1*py*z3 - x1*y3*pz - px*y1*z3 + px*y3*z1 + x3*y1*pz - x3*py*z1;
  d3 = x1*y2*pz - x1*py*z2 - x2*y1*pz + x2*py*z1 + px*y1*z2 - px*y2*z1;

  x2 -= x1; y2 -= y1; z2 -= z1;
  x3 -= x1; y3 -= y1; z3 -= z1;
  px -= x1; py -= y1; pz -= z1;

  d4 = -x2*y3*pz + x2*py*z3 + x3*y2*pz - x3*py*z2 - px*y2*z3 + px*y3*z2;
  if (d0 < 0.0)
    return (d1 < 0.0) && (d2 < 0.0) && (d3 < 0.0) && (d4 < 0.0);
  else
    return (d1 > 0.0) && (d2 > 0.0) && (d3 > 0.0) && (d4 > 0.0);
}


//--------------------------------------------------------------------
void Tetra::draw()
//--------------------------------------------------------------------
{
  int i,j, side;
  Vector3d center; getCenter(center);
  float r = 0.5*_bounds.diagonal();

  if (draw3d.drawStress) {
    float eigenVal;
    Vector3d eigenVec;
    computeStress(eigenVal, eigenVec);
    eigenVec.scale(r* eigenVal/physicsParams.youngModulus);
    draw3d.arrow(center, eigenVec, 0.1*r, DRAW_RED);
  }

  if (draw3d.drawForces && !draw3d.drawForceSum) {
    for (i = 0; i < 4; i++) {
      Vector3d f; elasticForce(i, f); f.scale(-0.01);
      draw3d.arrow(_vertex[i]->_worldCoord, f, 0.1*r, DRAW_BLUE);
    }
  }
  
  if (draw3d.drawOrientations) {
    for (int i = 0; i < 3; i++) {
      Vector3d n;
      _orientation.getColumn(i, n);
      n.scale(1.0*r);
      draw3d.arrow(center, n, 0.05*r, DRAW_YELLOW);
    }
  }

  if (draw3d.drawTetras) {
    float color[4];
    float s = 1.0/256.0;
    color[0] = s*COLOR_RED;
    color[1] = s*COLOR_GREEN;
    color[2] = s*COLOR_BLUE;
    color[3] = 1.0;

    for (side = 0; side < 4; side++) {
      Vertex* v[3];
      for (i = 0; i < 3; i++)
        v[i] = sideVertex(side, i);

      Vector3d v1,v2,normal;
      v1.sub(v[1]->_worldCoord, v[0]->_worldCoord);
      v2.sub(v[2]->_worldCoord, v[0]->_worldCoord);
      normal.cross(v1,v2); normal.normalize();

      Vector3d pos[3];
      for (i = 0; i < 3; i++) {
        pos[i].sub(v[i]->_worldCoord, center);
        pos[i].scale(0.7); pos[i].add(center);
      }
      draw3d.triangle(pos[0], pos[1], pos[2], normal, color);
    }
  }
}




