#include "Bounds3d.h"

//---------------------------------------------------------------------------

#include "Vector3d.h"

//-----------------------------------------------------------
void Bounds3d::include(Vector3d &v)
//-----------------------------------------------------------
{
  if (_empty) {
    _max = v;
    _min = v;
    _empty = false;
  }
  else {
    _max.max(v);
    _min.min(v);
  }
}


//-----------------------------------------------------------
void Bounds3d::intersect(Bounds3d &b1, Bounds3d &b2)
//-----------------------------------------------------------
{
  if (b1._min.x > b2._min.x) _min.x = b1._min.x; else _min.x = b2._min.x;
  if (b1._max.x < b2._max.x) _max.x = b1._max.x; else _max.x = b2._max.x;
  if (_max.x < _min.x) { _empty = true; return; }

  if (b1._min.y > b2._min.y) _min.y = b1._min.y; else _min.y = b2._min.y;
  if (b1._max.y < b2._max.y) _max.y = b1._max.y; else _max.y = b2._max.y;
  if (_max.y < _min.y) { _empty = true; return; }

  if (b1._min.z > b2._min.z) _min.z = b1._min.z; else _min.z = b2._min.z;
  if (b1._max.z < b2._max.z) _max.z = b1._max.z; else _max.z = b2._max.z;
  if (_max.z < _min.z) { _empty = true; return; }
}


//-----------------------------------------------------------
void Bounds3d::combine(Bounds3d &b1, Bounds3d &b2)
//-----------------------------------------------------------
{
  if (b1._empty && b2._empty) {
    _empty = true; return;
  }
  else if (b1._empty) {
    _empty = false;
    _min = b2._min; _max = b2._max;
  }
  else if (b2._empty) {
    _empty = false;
    _min = b1._min; _max = b1._max;
  }
  else {
    _empty = false;
    _min = b1._min; _min.min(b2._min);
    _max = b1._max; _max.max(b2._max);
  }  
}


//-----------------------------------------------------------
bool Bounds3d::intersect(Bounds3d &b)
//-----------------------------------------------------------
{
  if (_empty || b._empty) return false;
  if ((b._min.x > _max.x) || (_min.x > b._max.x)) return false;
  if ((b._min.y > _max.y) || (_min.y > b._max.y)) return false;
  if ((b._min.z > _max.z) || (_min.z > b._max.z)) return false;
  return true;
}


//-----------------------------------------------------------
bool Bounds3d::contain(Vector3d &v)
//-----------------------------------------------------------
{
  if (_empty) return false;
  if ((v.x < _min.x) || (v.x > _max.x)) return false;
  if ((v.y < _min.y) || (v.y > _max.y)) return false;
  if ((v.z < _min.z) || (v.z > _max.z)) return false;
  return true;
}

//-----------------------------------------------------------
void Bounds3d::scale(float s, bool centered)
//-----------------------------------------------------------
{
  if (_empty) return;
  if (centered) {
    Vector3d center; center.add(_min, _max);
    center.scale(0.5);
    Vector3d r;
    r.sub(_min, center); r.scale(s); _min.add(center, r);
    r.sub(_max, center); r.scale(s); _max.add(center, r);
  }
  else {
    _min.scale(s);
    _max.scale(s);
  }    
}
