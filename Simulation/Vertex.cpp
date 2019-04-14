//---------------------------------------------------------------------------

#include <math.h>
#include <assert.h>

#include "Vertex.h"
#include "Vector3d.h"
#include "Matrix3d.h"
#include "Draw.h"
#include "PhysicsParams.h"

//---------------------------------------------------------------------------

#define GS_OVER_RELAX 1.2 // Gauss Seidel overrelaxation

#define LOWER_DIST_FACTOR2 (0.5 * 0.5)  // preventing material
#define UPPER_DIST_FACTOR2 (1.5 * 1.5)  // from overstretching
#define CORRECTION_OVER_RELAX 1.4

#define F_COLLISION_RANGE 0.05          // [m] for computation of
                                        // derivative of collision force

//----------------------------------------------------
Vertex::Vertex(float x, float y, float z)
//----------------------------------------------------
{
  _originalCoord.set(x,y,z);
  _worldCoord.set(x,y,z);

  _animated = true;
  _velocity.setZero();
  _mass = 0.0;
  _forceDistr = 0.0;
  _fExt.setZero();
  _fColl.setZero();
  _f0.setZero();

  _numNeighbors = 0;

  _solution.setZero();
  _r.setZero(); _p.setZero(); _u.setZero();

  _mark = -1;
  _fractureMark = 0;
  _fractureTip = false;

  _colliding = false;
  _collDisp.setZero();
}

//----------------------------------------------------
Vertex::~Vertex()
//----------------------------------------------------
{
  removeNeighbors();
}


//----------------------------------------------------
void Vertex::transferInfo(Vertex *v)
//----------------------------------------------------
// for vertex split
{
  // cloned values
  _animated = v->_animated;
  _originalCoord = v->_originalCoord;
  _worldCoord = v->_worldCoord;
  _velocity = v->_velocity;

  _normal = v->_normal;
  _forceDistr = v->_forceDistr;
  _fExt = v->_fExt;
  _fColl = v->_fColl;
  _solution = v->_solution;
  _equation_b = v->_equation_b;
  _fractureMark = v->_fractureMark;

  // init values
  _mass = 0.0;
  _mark = -1;
  _numNeighbors = 0;
  _fractureTip = false;
}


//--------------------------------------------------------
void Vertex::transferInfo(Vertex *v1, Vertex *v2, float f)
//--------------------------------------------------------
// for edge split
{
  // cloned values
  _animated = v1->_animated || v2->_animated;
  _originalCoord.interpolate(v1->_originalCoord, v2->_originalCoord, f);
  _worldCoord.interpolate(v1->_worldCoord, v2->_worldCoord, f);
  _velocity.interpolate(v1->_velocity, v2->_velocity, f);

  _normal.interpolate(v1->_normal, v1->_normal, f);
  _forceDistr = (1.0 - f)*v1->_forceDistr + f*v2->_forceDistr;
  _fExt.interpolate(v1->_fExt, v2->_fExt, f);
  _fColl.interpolate(v1->_fColl, v2->_fColl, f);
  _solution.interpolate(v1->_solution, v2->_solution, f);
  _equation_b.interpolate(v1->_equation_b, v2->_equation_b, f);
  _fractureMark = v1->_fractureMark;

  // init values
  _mass = 0.0;
  _mark = -1;
  _numNeighbors = 0;
  _fractureTip = false;
}


//--------------------------------------------------------------------
void Vertex::elasticForce(Vector3d &f)
//--------------------------------------------------------------------
{
  f.setZero();
  for (int i = 0; i < _numNeighbors; i++) {
    VertexNeighbor *n = _neighbor[i];
    Vector3d p;
    p.mult(n->_stiffness, n->_vertex->_worldCoord);
    f.add(p);
  }
  f.add(_f0);
}


//----------------------------------------------------
void Vertex::addNeighborReference(Vertex *v)
//----------------------------------------------------
{
  int i = 0;
  while ((i < _numNeighbors) && (_neighbor[i]->_vertex != v)) i++;
  if (i < _numNeighbors) {
    _neighbor[i]->_refCount++;
  }
  else if (_numNeighbors < MAX_NEIGHBORS) {
    i = _numNeighbors;
    _numNeighbors++;
    VertexNeighbor *n = new VertexNeighbor();
    n->_vertex = v;
    n->_stiffness.setZero();
    n->_refCount = 1;
    _neighbor[i] = n;
  }
  else {
    assert(false);  // too many neighbors
  }
}

//-------------------------------------------------------
void Vertex::subNeighborReference(Vertex *v)
//-------------------------------------------------------
{
  int i = 0;
  while ((i < _numNeighbors) && (_neighbor[i]->_vertex != v)) i++;
  if (i < _numNeighbors) {
    VertexNeighbor *n = _neighbor[i];
    n->_refCount--;
    if (n->_refCount <= 0) {
      delete _neighbor[i];
      _neighbor[i] = _neighbor[_numNeighbors-1];
      _numNeighbors--;
    }
  }
  else
    assert(false);  // neighbor not found
}


//----------------------------------------------------
void Vertex::removeNeighbors()
//----------------------------------------------------
{
  for (int i = 0; i < _numNeighbors; i++)
    delete _neighbor[i];
  _numNeighbors = 0;
}


//----------------------------------------------------
void Vertex::addStiffness(Vertex *v, Matrix3d &m)
//----------------------------------------------------
{
  int i = 0;
  while ((i < _numNeighbors) && (_neighbor[i]->_vertex != v)) i++;
  assert(i < _numNeighbors);
  _neighbor[i]->_stiffness.add(m);
}


//----------------------------------------------------
void Vertex::subStiffness()
//----------------------------------------------------
{
  for (int i = 0; i < _numNeighbors; i++) {
    VertexNeighbor *n = _neighbor[i];
    n->_stiffness.setZero();
  }
}


//--------------------------------------------------------------------------
void Vertex::setupDynamicEquation(float delta_t, float damping)
//--------------------------------------------------------------------------
{
  if (_numNeighbors <= 0) return;
  Vector3d vec;

  _equation_b.setZero();                        // right hand side
  for (int j = 0; j < _numNeighbors; j++) {
    VertexNeighbor *n = _neighbor[j];
    Vertex *vj    = n->_vertex;
    Matrix3d &Aij = n->_equation_A;

    Aij.set(n->_stiffness);
    vec.mult(Aij, vj->_worldCoord);
    _equation_b.sub(vec);
    Aij.scale(delta_t * delta_t);
    if (vj == this) {
      Aij.addDiagonal(_mass + delta_t * damping *_mass);
    }
  }
  _equation_b.sub(_f0);
  _equation_b.add(_fExt);
  _equation_b.scale(delta_t);                   // dt * f + m * v

  vec = _velocity; vec.scale(_mass);
  _equation_b.add(vec);
}


//--------------------------------------------------------------------------
void Vertex::setupStaticEquation()
//--------------------------------------------------------------------------
{
  if (_numNeighbors <= 0) return;

  _equation_b.setZero();

  for (int j = 0; j < _numNeighbors; j++) {
    VertexNeighbor *n = _neighbor[j];
    n->_equation_A.set(n->_stiffness);
  }
  _equation_b.sub(_f0);
  _equation_b.add(_fExt);               // external forces
}


/* ------------------------------------------------------------------------ */
void Vertex::multiply(bool with_p)
/* ------------------------------------------------------------------------ */
{
  // used in conjugate gradients
  // with_p    : u = A p
  // otherwise : p = r = b - A * solution

  if (_numNeighbors <= 0) return;
  Vector3d sum, vec;

  sum.set(0,0,0);
  for (int j = 0; j < _numNeighbors; j++) {
    VertexNeighbor *n = _neighbor[j];
    Vertex *vj    = n->_vertex;
    Matrix3d &Aij = n->_equation_A;

    if (with_p)
      vec.mult(Aij, vj->_p);
    else
      vec.mult(Aij, vj->_solution);
    sum.add(vec);
  }

  if (with_p)
    _u = sum;
  else {
    _r.sub(_equation_b, sum);
    _p = _r;
  }
}

/* ------------------------------------------------------------------------ */
void Vertex::gaussSeidel()
/* ------------------------------------------------------------------------ */
{
  if (_numNeighbors <= 0) return;

  Vector3d vec;
  Vector3d sum = _equation_b;

  Matrix3d K;

  for (int j = 0; j < _numNeighbors; j++) {
    VertexNeighbor *n = _neighbor[j];
    Vertex *vj    = n->_vertex;

    if (vj == this)
      K.set(n->_equation_A);
    else {
      vec.mult(n->_equation_A, vj->_solution);
      sum.sub(vec);
    }
  }

  float x = sum.x;
  float y = sum.y;
  float z = sum.z;
  float d,d1,d2,d3;

  d  = K.a00*K.a11*K.a22 + K.a01*K.a12*K.a20 + K.a02*K.a10*K.a21
     - K.a02*K.a11*K.a20 - K.a01*K.a10*K.a22 - K.a00*K.a12*K.a21;

  d1 =     x*K.a11*K.a22 + K.a01*K.a12*z     + K.a02*y    *K.a21
     - K.a02*K.a11*z     - K.a01*y    *K.a22 - x    *K.a12*K.a21;

  d2 = K.a00*y  *K.a22 +   x*K.a12*K.a20 + K.a02*K.a10*z
     - K.a02*y  *K.a20 -   x*K.a10*K.a22 - K.a00*K.a12*z;

  d3 = K.a00*K.a11*z     + K.a01*y    *K.a20 +     x*K.a10*K.a21
     -     x*K.a11*K.a20 - K.a01*K.a10*z     - K.a00*y    *K.a21;

  float dx = d1/d - _solution.x;
  float dy = d2/d - _solution.y;
  float dz = d3/d - _solution.z;

  _solution.x += GS_OVER_RELAX * dx;
  _solution.y += GS_OVER_RELAX * dy;
  _solution.z += GS_OVER_RELAX * dz;
}



// ---- prevent material from overstretching -------------------------------

void computeCorrection(Vector3d &p1, Vector3d &p2, float d2, Vector3d &correction)
{
  correction.sub(p2,p1);
  float len = sqrt(correction.len2());
  float d = sqrt(d2);
  if (len == 0.0) { correction.set(d,0,0); return; }
  float k = 1.0 - d/len;
  correction.scale(k);
}

// ------------------------------------------------------------------------
bool Vertex::preventOverstreching()
// ------------------------------------------------------------------------
{
  if (_numNeighbors <= 0) return false;

  int num = 0;
  Vector3d sum; sum.set(0,0,0);

  for (int i = 0; i < _numNeighbors; i++) {
    Vertex *v = _neighbor[i]->_vertex;
    if (v == this) continue;

    float olen2 = _originalCoord.dist2(v->_originalCoord);
    float lower2 = LOWER_DIST_FACTOR2 * olen2;
    float upper2 = UPPER_DIST_FACTOR2 * olen2;
    float len2  = _worldCoord.dist2(v->_worldCoord);
    Vector3d correction;
    if (len2 < lower2) {
      computeCorrection(_worldCoord, v->_worldCoord, lower2, correction);
      sum.add(correction); num++;
    }
    else if (len2 > upper2) {
      computeCorrection(_worldCoord, v->_worldCoord, upper2, correction);
      sum.add(correction); num++;
    }
  }
  if (num > 0) {
    sum.scale(1.0/num * CORRECTION_OVER_RELAX);
    _worldCoord.add(sum);
  }
  return (num > 0);
}


/* ------------------------------------------------------------------------ */
void Vertex::draw(float size)
/* ------------------------------------------------------------------------ */
{
  if (draw3d.drawVertices) {
    if (!_animated)
      draw3d.cube(_worldCoord, size, DRAW_BLACK);
    else if (_fractureTip)
      draw3d.cube(_worldCoord, size, DRAW_BLUE);
    else
      draw3d.cube(_worldCoord, size, DRAW_RED);
  }
  if (draw3d.drawForces && draw3d.drawForceSum) {
    Vector3d f; elasticForce(f); f.scale(-0.01);
    draw3d.arrow(_worldCoord, f, 0.5*size, DRAW_BLUE);
  }
}


