//---------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include "Vector3d.h"
#include "Matrix3d.h"

//---------------------------------------------------------------------------

float _sqrt(float f)
{
  if (f < 0.0 || f > 1e20) return 0.0;
  return sqrt(f);
}

float Vector3d::inner(Vector3d &v)
  { return x*v.x + y*v.y + z*v.z; }

float Vector3d::dist2(Vector3d &v)
{
  float dx,dy,dz;
  dx = x-v.x; dy = y-v.y; dz = z-v.z;
  return dx*dx + dy*dy + dz*dz;
}

float Vector3d::dist(Vector3d &v)
{
  float dx,dy,dz;
  dx = x-v.x; dy = y-v.y; dz = z-v.z;
  return _sqrt(dx*dx + dy*dy + dz*dz);
}

float Vector3d::len2()
{
  return x*x + y*y + z*z;
}

float Vector3d::len()
{
  return _sqrt(x*x + y*y + z*z);
}

void Vector3d::mult(Matrix3d &m, Vector3d &v)
{
  float bx,by,bz;
  bx = m.a00 * v.x + m.a01 * v.y + m.a02 * v.z;
  by = m.a10 * v.x + m.a11 * v.y + m.a12 * v.z;
  bz = m.a20 * v.x + m.a21 * v.y + m.a22 * v.z;
  x = bx; y = by; z = bz;
}

void Vector3d::multTranspose(Matrix3d &m, Vector3d &v)
{
  float bx,by,bz;
  bx = m.a00 * v.x + m.a10 * v.y + m.a20 * v.z;
  by = m.a01 * v.x + m.a11 * v.y + m.a21 * v.z;
  bz = m.a02 * v.x + m.a12 * v.y + m.a22 * v.z;
  x = bx; y = by; z = bz;
}


void Vector3d::add(Vector3d &v)
{
  x += v.x; y += v.y; z += v.z;
}

void Vector3d::add(Vector3d &v1, Vector3d &v2)
{
  x = v1.x + v2.x;
  y = v1.y + v2.y;
  z = v1.z + v2.z;
}

void Vector3d::sub(Vector3d &v1, Vector3d &v2)
{
  x = v1.x - v2.x;
  y = v1.y - v2.y;
  z = v1.z - v2.z;
}

void Vector3d::sub(Vector3d &v)
{
  x -= v.x; y -= v.y; z -= v.z;
}

void Vector3d::cross(Vector3d &v1, Vector3d &v2)
{
  float bx,by,bz;
  bx = v1.y * v2.z - v2.y * v1.z;
  by = v1.z * v2.x - v2.z * v1.x;
  bz = v1.x * v2.y - v2.x * v1.y;
  x = bx; y = by; z = bz;
}

void Vector3d::scale(float f)
{
  x *= f; y *= f; z *= f;
}

void Vector3d::normalize()
{
  float n = _sqrt(len2());
  if (n == 0.0) {
    x = 0.0; y = -1.0; z = 0.0;
  }
  else scale(1.0/n);
}

void Vector3d::min(Vector3d &v)
{
  if (v.x < x) x = v.x;
  if (v.y < y) y = v.y;
  if (v.z < z) z = v.z;
}

void Vector3d::max(Vector3d &v)
{
  if (v.x > x) x = v.x;
  if (v.y > y) y = v.y;
  if (v.z > z) z = v.z;
}

void Vector3d::reflect(Vector3d &n, float e)
/* e = 0.0 inelastic, e = 1.0 perfectly elastic */
{
  Vector3d d = n; d.normalize();
  float s = inner(d);
  if (s >= 0.0) return;
  d.scale(s * (1.0 + e));
  sub(d);
}

void Vector3d::damp(Vector3d &n, float f)
/* f = 0.0: damp to zero, f = 1.0: don' damp */
{
  Vector3d d = n; d.normalize();
  float s = inner(d);
  d.scale(s);
  Vector3d p; p.sub(*this, d);
  p.scale(f);
  add(d, p);
}


void rotate(float &x, float &y, float angle)
{
  float sinf = sin(angle);
  float cosf = cos(angle);
  float x1 = cosf * x - sinf * y;
  float y1 = sinf * x + cosf * y;
  x = x1; y = y1;
}

void Vector3d::rotateX(float angle)
{
  rotate(y, z, angle);
}

void Vector3d::rotateY(float angle)
{
  rotate(z, x, angle);
}

void Vector3d::rotateZ(float angle)
{
  rotate(x, y, angle);
}

void Vector3d::interpolate(Vector3d &v1, Vector3d &v2, float f)
{
  if (f < 0.0) f = 0.0;
  if (f > 1.0) f = 1.0;
  float g = 1.0 - f;
  x = g * v1.x + f * v2.x;
  y = g * v1.y + f * v2.y;
  z = g * v1.z + f * v2.z;
}


void Vector3d::show() {
  printf("%6.3f %6.3f %6.3f\n\n", x,y,z);
}


//--------------------------------------------------------------------
void normalize(Vector3d &m0, Vector3d &m1, Vector3d &m2, Vector3d &m3)
//--------------------------------------------------------------------
{
  // translate m0 to origin
  // rotate    m1 to x - axis
  // rotate    m2 to x-y plane

  float x,y,z, r, cosf, sinf;
  float m1x, m1y, m1z, m2x, m2y, m2z, m3x, m3y, m3z;

  /* translate m0 to origin */

  x = m0.x; y = m0.y; z = m0.z;
  m1x = m1.x - x; m1y = m1.y - y; m1z = m1.z - z;
  m2x = m2.x - x; m2y = m2.y - y; m2z = m2.z - z;
  m3x = m3.x - x; m3y = m3.y - y; m3z = m3.z - z;

  /* xy rotation */

  r = _sqrt( m1x*m1x + m1y*m1y );
  if (r < EPSILON) { cosf = 1.0; sinf = 0.0; }
  else { cosf = m1x/r; sinf = -m1y/r; }
  m1x = r; m1y = 0.0;
  x    = cosf * m2x - sinf * m2y;
  y    = sinf * m2x + cosf * m2y;
  m2x = x; m2y = y;
  x    = cosf * m3x - sinf * m3y;
  y    = sinf * m3x + cosf * m3y;
  m3x = x; m3y = y;

  /* xz rotation */

  r = _sqrt( m1x*m1x + m1z*m1z );
  if (r < EPSILON) { cosf = 1.0; sinf = 0.0; }
  else { cosf = m1x/r; sinf = -m1z/r; }
  m1x = r; m1z = 0.0;
  x    = cosf * m2x - sinf * m2z;
  z    = sinf * m2x + cosf * m2z;
  m2x = x; m2z = z;
  x    = cosf * m3x - sinf * m3z;
  z    = sinf * m3x + cosf * m3z;
  m3x = x; m3z = z;

  /* yz rotation */

  r =  _sqrt( m2y*m2y + m2z*m2z );
  if (r < EPSILON) { cosf = 1.0; sinf = 0.0; }
  else { cosf = m2y/r; sinf = -m2z/r; }
  m2y = r; m2z = 0.0;
  y    = cosf * m3y - sinf * m3z;
  z    = sinf * m3y + cosf * m3z;
  m3y = y; m3z = z;

  m0.x = 0.0; m0.y = 0.0; m0.z = 0.0;
  m1.x = m1x; m1.y = m1y; m1.z = m1z;
  m2.x = m2x; m2.y = m2y; m2.z = m2z;
  m3.x = m3x; m3.y = m3y; m3.z = m3z;
}


//--------------------------------------------------------------------
bool triangleLineCut(Vector3d &a, Vector3d &b, Vector3d &c,
        Vector3d &p1, Vector3d &p2, float &f)
//--------------------------------------------------------------------
{
  // does line p1-p2 cut triangle a,b,c?
  // p1 + f *(p2-p1) is the intersection point

  bool cuts = true;

  Vector3d n, d1,d2,d3;
  f = 0.0;
  d1.sub(b,a); d2.sub(c,a);
  n.cross(d1,d2);

  Matrix3d m; m.setColumns(d1, d2, n);
  m.invert();

  d1.sub(a, p1); d2.sub(p2,p1);
  f = n.inner(d2);
  if (f == 0.0) return false;
  f = n.inner(d1)/f;
  if (f > 1.0 || f < 0.0) cuts = false;

  Vector3d p;
  d2.scale(f);
  p.add(p1, d2);
  p.sub(a);
  p.mult(m,p);
  if ((p.x <= 0.0) || (p.y <= 0.0) || (p.x + p.y >= 1.0)) cuts = false;
/*
  Vector3d n1,n2,n3, scp;
  d2.scale(f);
  scp.add(p1, d2);
  d1.sub(a, scp); d2.sub(b, scp); d3.sub(c, scp);
  n1.cross(d1,d2); n2.cross(d2,d3); n3.cross(d3,d1);
  if ((n1.inner(n2) < 0.0) || (n2.inner(n3) < 0.0)) cuts = false;
*/
  return cuts;
}

