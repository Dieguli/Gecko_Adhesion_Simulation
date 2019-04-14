//---------------------------------------------------------------------------

#ifndef ViewH
#define ViewH

#include "Vector3d.h"
//---------------------------------------------------------------------------


// -------------------------------------------------------
class View {
// -------------------------------------------------------
// handles mouse movements to change camera position

public:
  View() { Reset(); }

  void Reset() {
    mouse_x = mouse_y = 0;
    x_rot = y_rot = 0;
    x_trans = y_trans = z_trans = 0;
    scale = 1.0;
  }

  void LookAt(float world_scale, Vector3d &world_center);

  void Click(int x, int y) {
    mouse_x = x;
    mouse_y = y;
  }

  // =========
  // MODIFIERS
  void DragRotate(int x, int y);
  void DragTranslate(int x, int y);
  void DragScale(int x, int y);

private:

  // mouse positions
  int mouse_x, mouse_y;

  // rotations
  float	x_rot, y_rot;

  // translation
  float x_trans, y_trans, z_trans;

  // scale
  float scale;
};


#endif
