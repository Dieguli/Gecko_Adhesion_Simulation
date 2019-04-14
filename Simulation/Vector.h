#ifndef VectorH
#define VectorH

#define VECTOR_EPS 1e-20
#include <math.h>
#include <stdio.h>

/*----------------- Matrix3d --------------------------------------- */

class Matrix3d {
public:

  void set_zero() {
    a00 = 0.0; a01 = 0.0; a02 = 0.0;
    a10 = 0.0; a11 = 0.0; a12 = 0.0;
    a20 = 0.0; a21 = 0.0; a22 = 0.0;
  }

  void set_unity() {
    a00 = 1.0; a01 = 0.0; a02 = 0.0;
    a10 = 0.0; a11 = 1.0; a12 = 0.0;
    a20 = 0.0; a21 = 0.0; a22 = 1.0;
  }

  void set_tilde(float x, float y, float z)
  {
    a00 =  0.0; a01 =   -z; a02 =    y;
    a10 =    z; a11 =  0.0; a12 =   -x;
    a20 =   -y; a21 =    x; a22 =  0.0;
  }

  Matrix3d() {
    set_zero();
  }

  ~Matrix3d() {}

  void mult(Matrix3d &m1, Matrix3d &m2) {
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

  float det() {
    return  a00*a11*a22 + a01*a12*a20 + a02*a10*a21
           -a02*a11*a20 - a01*a10*a22 - a00*a12*a21;
  }

  void invert() {
    float b00,b01,b02,b10,b11,b12,b20,b21,b22;
    float d = det();
//    assert(fabs(d) > 0.0);
    if (d == 0.0) {	// hack
      set_unity();
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

  void transpose() {
    float a;
    a = a01; a01 = a10; a10 = a;
    a = a02; a02 = a20; a20 = a;
    a = a12; a12 = a21; a21 = a;
  }

  void rotate(float &x, float &y, float angle)
  {
    float sinf = sin(angle);
    float cosf = cos(angle);
    float x1 = cosf * x - sinf * y;
    float y1 = sinf * x + cosf * y;
    x = x1; y = y1;
  }

  void rotate_x(float angle)
  {
    rotate(a10, a20, angle);
    rotate(a11, a21, angle);
    rotate(a12, a22, angle);
  }

  void rotate_y(float angle)
  {
    rotate(a20, a00, angle);
    rotate(a21, a01, angle);
    rotate(a22, a02, angle);
  }

  void rotate_z(float angle)
  {
    rotate(a00, a10, angle);
    rotate(a01, a11, angle);
    rotate(a02, a12, angle);
  }


  void add(Matrix3d &m) {
    a00 += m.a00; a01 += m.a01; a02 += m.a02;
    a10 += m.a10; a11 += m.a11; a12 += m.a12;
    a20 += m.a20; a21 += m.a21; a22 += m.a22;
  }

  void scale(float f) {
    a00 *=f; a01 *= f; a02 *=f;
    a10 *=f; a11 *= f; a12 *=f;
    a20 *=f; a21 *= f; a22 *=f;
  }

  void show() {
    printf("%6.3f %6.3f %6.3f\n", a00, a01, a02);
    printf("%6.3f %6.3f %6.3f\n", a10, a11, a12);
    printf("%6.3f %6.3f %6.3f\n\n", a20, a21, a22);
  }

  float a00,a01,a02;
  float a10,a11,a12;
  float a20,a21,a22;
};


/*-----------------FEMVec ------------------------------------------ */

class Vector3d {
public:
  Vector3d() {}
  ~Vector3d() {}
  Vector3d& operator=(const Vector3d &v) {
    x = v.x; y = v.y; z = v.z;
    return *this;
  }

  void set(float new_x, float new_y, float new_z)
  { x = new_x, y = new_y; z = new_z; }

  float inner(Vector3d &v)
  { return x*v.x + y*v.y + z*v.z; }

  float dist2(Vector3d &v)
  {
    float dx,dy,dz;
    dx = x-v.x; dy = y-v.y; dz = z-v.z;
    return dx*dx + dy*dy + dz*dz;
  }

  float dist(Vector3d &v)
  {
    float dx,dy,dz;
    dx = x-v.x; dy = y-v.y; dz = z-v.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
  }

  float len2()
  {
    return x*x + y*y + z*z;
  }

  float len()
  {
    return sqrt(x*x + y*y + z*z);
  }

  void mult(Matrix3d &m, Vector3d &v)
  {
    float bx,by,bz;
    bx = m.a00 * v.x + m.a01 * v.y + m.a02 * v.z;
    by = m.a10 * v.x + m.a11 * v.y + m.a12 * v.z;
    bz = m.a20 * v.x + m.a21 * v.y + m.a22 * v.z;
    x = bx; y = by; z = bz;
  }

  void add(Vector3d &v)
  {
    x += v.x; y += v.y; z += v.z;
  }

  void add(Vector3d &v1, Vector3d &v2)
  {
    x = v1.x + v2.x;
    y = v1.y + v2.y;
    z = v1.z + v2.z;
  }

  void sub(Vector3d &v1, Vector3d &v2)
  {
    x = v1.x - v2.x;
    y = v1.y - v2.y;
    z = v1.z - v2.z;
  }

  void sub(Vector3d &v)
  {
    x -= v.x; y -= v.y; z -= v.z;
  }

  void cross(Vector3d &v1, Vector3d &v2)
  {
    float bx,by,bz;
    bx = v1.y * v2.z - v2.y * v1.z;
    by = v1.z * v2.x - v2.z * v1.x;
    bz = v1.x * v2.y - v2.x * v1.y;
    x = bx; y = by; z = bz;
  }

  void scale(float f)
  {
    x *= f; y *= f; z *= f;
  }

  void normalize()
  {
    float n = sqrt(len2());
    if (n == 0.0) {
      x = 0.0; y = -1.0; z = 0.0;
      printf("Vector3d::normalize() vec is zero!\n");
    }
    else scale(1.0/n);
  }

  void min(Vector3d &v)
  {
    if (v.x < x) x = v.x;
    if (v.y < y) y = v.y;
    if (v.z < z) z = v.z;
  }

  void max(Vector3d &v)
  {
    if (v.x > x) x = v.x;
    if (v.y > y) y = v.y;
    if (v.z > z) z = v.z;
  }

  void reflect(Vector3d &n, float e)
  /* e = 0.0 inelastic, e = 1.0 perfectly elastic */
  {
    Vector3d d = n; d.normalize();
    float s = inner(d);
    if (s >= 0.0) return;
    d.scale(s * (1.0 + e));
    sub(d);
  }

  void rotate(float &x, float &y, float angle)
  {
    float sinf = sin(angle);
    float cosf = cos(angle);
    float x1 = cosf * x - sinf * y;
    float y1 = sinf * x + cosf * y;
    x = x1; y = y1;
  }

  void rotate_x(float angle)
  {
    rotate(y, z, angle);
  }

  void rotate_y(float angle)
  {
    rotate(z, x, angle);
  }

  void rotate_z(float angle)
  {
    rotate(x, y, angle);
  }

  void show() {
    printf("%6.3f %6.3f %6.3f\n\n", x,y,z);
  }

  float x,y,z;
};


#endif
