//---------------------------------------------------------------------------

#ifndef Matrix3dH
#define Matrix3dH

class Vector3d;

//---------------------------------------------------------------------------

class Matrix3d {
public:

  Matrix3d() { setZero(); }
  ~Matrix3d() {}

  void set(Matrix3d &m);
  void setZero();
  void setUnity();
  void setTilde(float x, float y, float z);
  void set(int i, int j, float d);
  void add(int i, int j, float d);

  float trace() { return a00 + a11 + a22; }
  void getColumn(int i, Vector3d &v);
  void setColumns(Vector3d &v1, Vector3d &v2, Vector3d &v3);

  void mult(Matrix3d &m);
  void multTranspose(Matrix3d &m);
  void mult(Matrix3d &m1, Matrix3d &m2);

  float frobeniusNorm();

  float det();
  void  invert();
  void transpose();

  void rotateX(float angle);
  void rotateY(float angle);
  void rotateZ(float angle);

  void add(Matrix3d &m);
  void addDiagonal(float f);
  void subDiagonal(float f);
  void sub(Matrix3d &m);
  void scale(float f);

  void orthonormalize();
  void gramSchmidt();
  void largestEigenvalue(float &eigenValue, Vector3d &eigenVec);
  void interpolate(Matrix3d &m1, Matrix3d &m2, float f);

  void show();

  // representation

  float a00,a01,a02;
  float a10,a11,a12;
  float a20,a21,a22;
};



#endif
