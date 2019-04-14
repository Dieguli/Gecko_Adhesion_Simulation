//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "View.h"
#include <math.h>
#include <gl\gl.h>
#include <gl\glu.h>

//---------------------------------------------------------------------------

//--------------------------------------------------------
void View::DragRotate(int x, int y)
//--------------------------------------------------------
{
  float tmp = (M_PI / 10);
  y_rot += tmp*(x - mouse_x);
  x_rot += tmp*(y - mouse_y);
  if (x_rot < -89.9999) x_rot = -89.9999;
  if (x_rot > 89.9999) x_rot = 89.9999;
  mouse_x = x;
  mouse_y = y;
}

//--------------------------------------------------------
void View::DragTranslate(int x, int y)
//--------------------------------------------------------
{
  int dx = x - mouse_x;
  int dy = y - mouse_y;

  Vector3d v; v.set(dx/5.0,-dy/5.0,0.0);
  v.rotateX(-x_rot*M_PI/180.0);
  v.rotateY(-y_rot*M_PI/180.0);

  x_trans += v.x;
  y_trans += v.y;
  z_trans += v.z;
  mouse_x = x;
  mouse_y = y;
}

//--------------------------------------------------------
void View::DragScale(int x, int y) {
//--------------------------------------------------------
  float tmp = 1 + (y-mouse_y)*0.01;
  scale *= tmp;
  x_trans *= tmp;
  y_trans *= tmp;
  z_trans *= tmp;
  mouse_x = x;
  mouse_y = y;
}

//----------------------------------------------------------
void View::LookAt(float world_scale, Vector3d &world_center)
//----------------------------------------------------------
{
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, 220,   0, 0, 0,   0, 1, 0);

  // rotate
  glRotatef(x_rot, 1.0, 0.0, 0.0);
  glRotatef(y_rot, 0.0, 1.0, 0.0);

  // translate - change what's at the center of the screen
  glTranslatef(x_trans, y_trans, z_trans);

  // scale
  float s = scale * world_scale;
  glScalef(s, s, s);

  // translate - center the world's bounding box
  glTranslatef(-world_center.x, -world_center.y, -world_center.z);
}



