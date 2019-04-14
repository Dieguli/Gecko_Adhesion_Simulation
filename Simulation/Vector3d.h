//---------------------------------------------------------------------------

#ifndef Vector3dH
#define Vector3dH

#define EPSILON 1e-9

class Matrix3d;

//---------------------------------------------------------------------------


//-------------------------------------------------
class Vector3d {
//-------------------------------------------------
public:
  Vector3d() {}
  ~Vector3d() {}

  Vector3d& operator=(const Vector3d &v) {
    x = v.x; y = v.y; z = v.z;
    return *this;
  }

  void set(float new_x, float new_y, float new_z)
  { x = new_x, y = new_y; z = new_z; }

  void setZero() { x = 0.0; y = 0.0; z = 0.0; }
  bool isZero() { return x == 0.0 && y == 0.0 && z == 0.0; }

  float dist2(Vector3d &v);
  float dist(Vector3d &v);
  float len2();
  float len();

  void mult(Matrix3d &m, Vector3d &v);
  void multTranspose(Matrix3d &m, Vector3d &v);
  void add(Vector3d &v);
  void add(Vector3d &v1, Vector3d &v2);
  void sub(Vector3d &v1, Vector3d &v2);
  void sub(Vector3d &v);

  void cross(Vector3d &v1, Vector3d &v2);
  float inner(Vector3d &v);

  void scale(float f);
  void normalize();
  void min(Vector3d &v);
  void max(Vector3d &v);
  void reflect(Vector3d &n, float e);
  void damp(Vector3d &n, float f);
  void rotateX(float angle);
  void rotateY(float angle);
  void rotateZ(float angle);
  void interpolate(Vector3d &v1, Vector3d &v2, float f);
  void show();

  // representation
  float x,y,z;
};

void normalize(Vector3d &m0, Vector3d &m1, Vector3d &m2, Vector3d &m3);
bool triangleLineCut(Vector3d &a, Vector3d &b, Vector3d &c,
        Vector3d &p1, Vector3d &p2, float &f);

#endif
