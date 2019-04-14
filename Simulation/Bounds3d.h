//---------------------------------------------------------------------------

#ifndef Bounds3dH
#define Bounds3dH
//---------------------------------------------------------------------------

#include "Vector3d.h"

//-------------------------------------------------
class Bounds3d {
//-------------------------------------------------
public:
  Bounds3d() { _empty = true; }
  ~Bounds3d() {}

  void clear() { _empty = true; }
  bool isEmpty() { return _empty; }

  void include(Vector3d &v);
  void intersect(Bounds3d &b1, Bounds3d &b2);
  void combine(Bounds3d &b1, Bounds3d &b2);

  bool intersect(Bounds3d &b);
  bool contain(Vector3d &v);
  float diagonal() { return _min.dist(_max); };
  void scale(float s, bool centered);

  Vector3d _min, _max;

private:
  bool _empty;
};

#endif
