//---------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include "Vector3d.h"
#include "Matrix3d.h"

//---------------------------------------------------
void Matrix3d::set(Matrix3d &m) {
//---------------------------------------------------
  a00 = m.a00; a01 = m.a01; a02 = m.a02;
  a10 = m.a10; a11 = m.a11; a12 = m.a12;
  a20 = m.a20; a21 = m.a21; a22 = m.a22;
}

//---------------------------------------------------
void Matrix3d::setZero() {
//---------------------------------------------------
  a00 = 0.0; a01 = 0.0; a02 = 0.0;
  a10 = 0.0; a11 = 0.0; a12 = 0.0;
  a20 = 0.0; a21 = 0.0; a22 = 0.0;
}

//---------------------------------------------------
void Matrix3d::setUnity() {
//---------------------------------------------------
  a00 = 1.0; a01 = 0.0; a02 = 0.0;
  a10 = 0.0; a11 = 1.0; a12 = 0.0;
  a20 = 0.0; a21 = 0.0; a22 = 1.0;
}

//---------------------------------------------------
void Matrix3d::setTilde(float x, float y, float z)
//---------------------------------------------------
{
  a00 =  0.0; a01 =   -z; a02 =    y;
  a10 =    z; a11 =  0.0; a12 =   -x;
  a20 =   -y; a21 =    x; a22 =  0.0;
}

//---------------------------------------------------
void Matrix3d::set(int i, int j, float d)
//---------------------------------------------------
{
  if (i == 0) {
    if (j == 0)      a00 = d;
    else if (j == 1) a01 = d;
    else if (j == 2) a02 = d;
  }
  else if (i == 1) {
    if (j == 0)      a10 = d;
    else if (j == 1) a11 = d;
    else if (j == 2) a12 = d;
  }
  else if (i == 2) {
    if (j == 0)      a20 = d;
    else if (j == 1) a21 = d;
    else if (j == 2) a22 = d;
  }
}

//---------------------------------------------------
void Matrix3d::add(int i, int j, float d)
//---------------------------------------------------
{
  if (i == 0) {
    if (j == 0)      a00 += d;
    else if (j == 1) a01 += d;
    else if (j == 2) a02 += d;
  }
  else if (i == 1) {
    if (j == 0)      a10 += d;
    else if (j == 1) a11 += d;
    else if (j == 2) a12 += d;
  }
  else if (i == 2) {
    if (j == 0)      a20 += d;
    else if (j == 1) a21 += d;
    else if (j == 2) a22 += d;
  }
}


//---------------------------------------------------
void Matrix3d::getColumn(int i, Vector3d &v)
//---------------------------------------------------
{
  if (i == 0) {
    v.x = a00; v.y = a10; v.z = a20;
  }
  if (i == 1) {
    v.x = a01; v.y = a11; v.z = a21;
  }
  if (i == 2) {
    v.x = a02; v.y = a12; v.z = a22;
  }
}

//---------------------------------------------------
void Matrix3d::setColumns(Vector3d &v1, Vector3d &v2, Vector3d &v3)
//---------------------------------------------------
{
  a00 = v1.x; a10 = v1.y; a20 = v1.z;
  a01 = v2.x; a11 = v2.y; a21 = v2.z;
  a02 = v3.x; a12 = v3.y; a22 = v3.z;
}

//---------------------------------------------------
void Matrix3d::mult(Matrix3d &m) {
//---------------------------------------------------
  float b00,b01,b02,b10,b11,b12,b20,b21,b22;
  b00 = a00 * m.a00 + a01 * m.a10 + a02 * m.a20;
  b01 = a00 * m.a01 + a01 * m.a11 + a02 * m.a21;
  b02 = a00 * m.a02 + a01 * m.a12 + a02 * m.a22;

  b10 = a10 * m.a00 + a11 * m.a10 + a12 * m.a20;
  b11 = a10 * m.a01 + a11 * m.a11 + a12 * m.a21;
  b12 = a10 * m.a02 + a11 * m.a12 + a12 * m.a22;

  b20 = a20 * m.a00 + a21 * m.a10 + a22 * m.a20;
  b21 = a20 * m.a01 + a21 * m.a11 + a22 * m.a21;
  b22 = a20 * m.a02 + a21 * m.a12 + a22 * m.a22;
  a00 = b00; a01 = b01; a02 = b02;
  a10 = b10; a11 = b11; a12 = b12;
  a20 = b20; a21 = b21; a22 = b22;
}

//---------------------------------------------------
void Matrix3d::multTranspose(Matrix3d &m) {
//---------------------------------------------------
  float b00,b01,b02,b10,b11,b12,b20,b21,b22;
  b00 = a00 * m.a00 + a01 * m.a01 + a02 * m.a02;
  b01 = a00 * m.a10 + a01 * m.a11 + a02 * m.a12;
  b02 = a00 * m.a20 + a01 * m.a21 + a02 * m.a22;

  b10 = a10 * m.a00 + a11 * m.a01 + a12 * m.a02;
  b11 = a10 * m.a10 + a11 * m.a11 + a12 * m.a12;
  b12 = a10 * m.a20 + a11 * m.a21 + a12 * m.a22;

  b20 = a20 * m.a00 + a21 * m.a01 + a22 * m.a02;
  b21 = a20 * m.a10 + a21 * m.a11 + a22 * m.a12;
  b22 = a20 * m.a20 + a21 * m.a21 + a22 * m.a22;
  a00 = b00; a01 = b01; a02 = b02;
  a10 = b10; a11 = b11; a12 = b12;
  a20 = b20; a21 = b21; a22 = b22;
}

//---------------------------------------------------
void Matrix3d::mult(Matrix3d &m1, Matrix3d &m2) {
//---------------------------------------------------
  float b00,b01,b02,b10,b11,b12,b20,b21,b22;
  b00 = m1.a00 * m2.a00 + m1.a01 * m2.a10 + m1.a02 * m2.a20;
  b01 = m1.a00 * m2.a01 + m1.a01 * m2.a11 + m1.a02 * m2.a21;
  b02 = m1.a00 * m2.a02 + m1.a01 * m2.a12 + m1.a02 * m2.a22;

  b10 = m1.a10 * m2.a00 + m1.a11 * m2.a10 + m1.a12 * m2.a20;
  b11 = m1.a10 * m2.a01 + m1.a11 * m2.a11 + m1.a12 * m2.a21;
  b12 = m1.a10 * m2.a02 + m1.a11 * m2.a12 + m1.a12 * m2.a22;

  b20 = m1.a20 * m2.a00 + m1.a21 * m2.a10 + m1.a22 * m2.a20;
  b21 = m1.a20 * m2.a01 + m1.a21 * m2.a11 + m1.a22 * m2.a21;
  b22 = m1.a20 * m2.a02 + m1.a21 * m2.a12 + m1.a22 * m2.a22;
  a00 = b00; a01 = b01; a02 = b02;
  a10 = b10; a11 = b11; a12 = b12;
  a20 = b20; a21 = b21; a22 = b22;
}


//---------------------------------------------------
float Matrix3d::det() {
//---------------------------------------------------
  return  a00*a11*a22 + a01*a12*a20 + a02*a10*a21
         -a02*a11*a20 - a01*a10*a22 - a00*a12*a21;
}

//---------------------------------------------------
void Matrix3d::invert() {
//---------------------------------------------------
  float b00,b01,b02,b10,b11,b12,b20,b21,b22;
  float d = det();
//    assert(fabs(d) > 0.0);
  if (d == 0.0) {	// hack
    setUnity();
    return;
  }

  d = 1.0/d;
  b00 = a11*a22-a12*a21;
  b01 = a02*a21-a01*a22;
  b02 = a01*a12-a02*a11;
  b10 = a12*a20-a10*a22;
  b11 = a00*a22-a02*a20;
  b12 = a02*a10-a00*a12;
  b20 = a10*a21-a11*a20;
  b21 = a01*a20-a00*a21;
  b22 = a00*a11-a01*a10;
  a00 = b00*d; a01 = b01*d; a02 = b02*d;
  a10 = b10*d; a11 = b11*d; a12 = b12*d;
  a20 = b20*d; a21 = b21*d; a22 = b22*d;
}

//---------------------------------------------------
void Matrix3d::transpose() {
//---------------------------------------------------
  float a;
  a = a01; a01 = a10; a10 = a;
  a = a02; a02 = a20; a20 = a;
  a = a12; a12 = a21; a21 = a;
}

//---------------------------------------------------
void rotate(float &x, float &y, float angle)
//---------------------------------------------------
{
  float sinf = sin(angle);
  float cosf = cos(angle);
  float x1 = cosf * x - sinf * y;
  float y1 = sinf * x + cosf * y;
  x = x1; y = y1;
}

//---------------------------------------------------
void Matrix3d::rotateX(float angle)
//---------------------------------------------------
{
  rotate(a10, a20, angle);
  rotate(a11, a21, angle);
  rotate(a12, a22, angle);
}

//---------------------------------------------------
void Matrix3d::rotateY(float angle)
//---------------------------------------------------
{
  rotate(a20, a00, angle);
  rotate(a21, a01, angle);
  rotate(a22, a02, angle);
}

//---------------------------------------------------
void Matrix3d::rotateZ(float angle)
//---------------------------------------------------
{
  rotate(a00, a10, angle);
  rotate(a01, a11, angle);
  rotate(a02, a12, angle);
}


//---------------------------------------------------
void Matrix3d::add(Matrix3d &m) {
//---------------------------------------------------
  a00 += m.a00; a01 += m.a01; a02 += m.a02;
  a10 += m.a10; a11 += m.a11; a12 += m.a12;
  a20 += m.a20; a21 += m.a21; a22 += m.a22;
}

//---------------------------------------------------
void Matrix3d::addDiagonal(float f) {
//---------------------------------------------------
  a00 += f; a11 += f; a22 += f;
}

//---------------------------------------------------
void Matrix3d::subDiagonal(float f) {
//---------------------------------------------------
  a00 -= f; a11 -= f; a22 -= f;
}

//---------------------------------------------------
void Matrix3d::sub(Matrix3d &m) {
//---------------------------------------------------
  a00 -= m.a00; a01 -= m.a01; a02 -= m.a02;
  a10 -= m.a10; a11 -= m.a11; a12 -= m.a12;
  a20 -= m.a20; a21 -= m.a21; a22 -= m.a22;
}

//---------------------------------------------------
void Matrix3d::scale(float f) {
//---------------------------------------------------
  a00 *=f; a01 *= f; a02 *=f;
  a10 *=f; a11 *= f; a12 *=f;
  a20 *=f; a21 *= f; a22 *=f;
}

//---------------------------------------------------
float Matrix3d::frobeniusNorm()
//---------------------------------------------------
{
  return sqrt(a00*a00 + a01*a01 + a02*a02 +
              a10*a10 + a11*a11 + a12*a12 +
              a20*a20 + a21*a21 + a22*a22);
}

//---------------------------------------------------
void Matrix3d::interpolate(Matrix3d &m1, Matrix3d &m2, float f)
//---------------------------------------------------
{
  if (f < 0.0) f = 0.0;
  if (f > 1.0) f = 1.0;
  float g = 1.0 - f;
  a00 = g * m1.a00 + f * m2.a00;
  a01 = g * m1.a01 + f * m2.a01;
  a02 = g * m1.a02 + f * m2.a02;
  a10 = g * m1.a10 + f * m2.a10;
  a11 = g * m1.a11 + f * m2.a11;
  a12 = g * m1.a12 + f * m2.a12;
  a20 = g * m1.a20 + f * m2.a20;
  a21 = g * m1.a21 + f * m2.a21;
  a22 = g * m1.a22 + f * m2.a22;
}

//---------------------------------------------------
void Matrix3d::show() {
//---------------------------------------------------
  printf("%6.3f %6.3f %6.3f\n", a00, a01, a02);
  printf("%6.3f %6.3f %6.3f\n", a10, a11, a12);
  printf("%6.3f %6.3f %6.3f\n\n", a20, a21, a22);
}

//---------------------------------------------------
void Matrix3d::orthonormalize()
//---------------------------------------------------
{
  Vector3d x,y,z;
  x.set(a00, a10, a20);
  y.set(a01, a11, a21);
  x.normalize();
  z.cross(x,y); z.normalize();
  y.cross(z,x); y.normalize();

  a00 = x.x; a01 = y.x; a02 = z.x;
  a10 = x.y; a11 = y.y; a12 = z.y;
  a20 = x.z; a21 = y.z; a22 = z.z;
}

//---------------------------------------------------
void Matrix3d::gramSchmidt()
//---------------------------------------------------
{
  Vector3d v1, v2, v3;
  Vector3d w1, w2, w3, s2, s3;
  getColumn(0, v1); getColumn(1, v2); getColumn(2, v3);
  w1 = v1; w1.normalize();
  s2 = w1; s2.scale(w1.inner(v2)); w2.sub(v2, s2); w2.normalize();
  s3 = w2; s3.scale(w2.inner(v3)); s3.add(s2); w3.sub(v3,s3); w3.normalize();
  setColumns(w1,w2,w3);
}


//---------------------------------------------------
void SolveCubic(float b, float c, float d, float *x1, float *x2, float *x3)
//---------------------------------------------------
{
/* solves x^3 + b*x^2 + c*x + d = 0, returns real part of solutions */
  float e,f,g,u,s, sq, cosf;
  float r,rx,phi, x,y;

  u = -b/3;
  e = 3*u*u + 2*u*b + c;
  f = u*u*u + u*u*b + u*c + d;
  s = -e/3;
  g = -e*e*e/27;

  sq = f*f - 4*g;
  if (sq >= 0) {
    sq = sqrt(sq);
    r = (-f + sq)/2; phi = 0.0;
    if (r < 0) { r = -r; phi = M_PI; }
    if (r > EPSILON) r = exp(1.0/3*log(r));
  }
  else {
    sq = sqrt(-sq);
    x = -f/2; y = sq/2;
    r = sqrt(x*x + y*y);
    if (r < EPSILON) { r = 0.0; phi = 0.0; }
    else {
      cosf = x/r;
      if (cosf > 1.0) cosf = 1.0;
      if (cosf < -1.0) cosf = -1.0;
      phi = acos(cosf)/3; r = exp(1.0/3*log(r));
    }
  }
  if (r > EPSILON) rx = r + s/r; else rx = 0.0;
  *x1 = rx*cos(phi) + u;
  *x2 = rx*cos(phi + 2.0*M_PI/3.0) + u;
  *x3 = rx*cos(phi + 4.0*M_PI/3.0) + u;
/* the y's are needed if the solutions are complex
   not needed here because the Eigenvalues of a symmetric matix
   are always real
  if (r > EPSILON) ry = r - s/r; else ry = 0.0;
  *y1 = ry*sin(phi);
  *y2 = ry*sin(phi + 2.0*M_PI/3.0);
  *y3 = ry*sin(phi + 4.0*M_PI/3.0);
*/
}

//---------------------------------------------------
void Matrix3d::largestEigenvalue(float &eigenValue, Vector3d &eigenVec)
//---------------------------------------------------
{
/* The proc returns lambda, the largest eigenvalue and
   a corresponing eigenvector */
  float inv1,inv2,inv3;
  float a[3][3];
  float x[3];
  float l1,l2,l3,l;
  int i,j, i0,i1, j0,j1;
  int mi0 = 0; int mi1 = 0; int mj0 = 0; int mj1 = 0; int mj2 = 0;
  float det,d0,d1, max, s;

  s = a00+a01+a02 + a10+a11+a12 + a20+a21+a22;
  if (fabs(s) < 1.0) s = 1.0;
  else s = 1.0/s;
  if (s < 0.0) s = -s;

  a[0][0] = s*a00; a[0][1] = s*a01; a[0][2] = s*a02;
  a[1][0] = s*a10; a[1][1] = s*a11; a[1][2] = s*a12;
  a[2][0] = s*a20; a[2][1] = s*a21; a[2][2] = s*a22;

  inv1 = a[0][0] + a[1][1] + a[2][2];
  inv2 = a[0][0]*a[1][1]-a[0][1]*a[1][0] +
         a[0][0]*a[2][2]-a[0][2]*a[2][0] +
         a[1][1]*a[2][2]-a[1][2]*a[2][1];
  inv3 = a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1]
        -a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2] - a[0][2]*a[1][1]*a[2][0];

  SolveCubic(-inv1,inv2,-inv3, &l1,&l2,&l3);

//  if (fabs(l1) > fabs(l2)) l = l1; else l = l2;
//  if (fabs(l3) > fabs(l)) l = l3;
  if (l1 > l2) l = l1; else l = l2;  /* tension only */
  if (l3 > l) l = l3;


  a[0][0] -= l; a[1][1] -= l; a[2][2] -= l;
  eigenValue = l/s;

  max = 0.0;
  i0 = 1; i1 = 2;
  for (i = 0 ; i < 3; i++) {
    if (i == 1) i0--; if (i == 2) i1--;
    j0 = 1; j1 = 2;
    for (j = 0; j < 3; j++) {
      if (j == 1) j0--; if (j == 2) j1--;
      det = fabs(a[i0][j0]*a[i1][j1] - a[i0][j1]*a[i1][j0]);
      if (det > max) {
        max = det;
        mi0 = i0; mi1 = i1;
        mj0 = j0; mj1 = j1; mj2 = 3-j0-j1;
      }
    }
  }
  if (max > EPSILON) {	/* single eigenvalue */
    x[mj2] = -1.0;
    det = a[mi0][mj0]*a[mi1][mj1] - a[mi0][mj1]*a[mi1][mj0];
    d0  = a[mi0][mj2]*a[mi1][mj1] - a[mi0][mj1]*a[mi1][mj2];
    d1  = a[mi0][mj0]*a[mi1][mj2] - a[mi0][mj2]*a[mi1][mj0];
    x[mj0] = d0/det;
    x[mj1] = d1/det;
  }
  else {
    max = 0.0;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        if (fabs(a[i][j]) > max) {
          mi0 = i; mj0 = j; max = fabs(a[i][j]);
        }
      }
    }
    if (max > EPSILON) { /* double eigenvalue */
      mj1 = mj0+1; if (mj1 > 2) mj1 = 0; mj2 = 3-mj0-mj1;
      x[mj1] = -1.0;
      x[mj0] = a[mi0][mj1]/a[mi0][mj0];
      x[mj2] = 0.0;
    }
    else {  /* triple eigenvalue */
      x[0] = 1.0; x[1] = 0.0; x[2] = 0.0;
    }
  }
  eigenVec.set(x[0], x[1], x[2]);
  eigenVec.normalize();
}

