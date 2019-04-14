#include <vcl.h>
#pragma hdrstop

#include <gl/gl.h>
#include <gl/glu.h>
#include <math.h>

#include "Draw.h"
#include "Vector3d.h"
#include "Utils.h"
#include "jpeg.hpp"
#include <dir.h>

#define PLANE_SIZE 1.0  // m
#define NUM_TILES 10

//---------------------------------------------------------------------------

struct TextureRef {
  char     filename[FILENAME_LEN];
  GLuint   texID;
  int      width, height;
};

//---------------------------------------------------------------------------

const GLfloat glColor[13][4] = {
  {1.0, 1.0, 1.0, 1.0},
  {0.316, 0.171, 0.156, 1.0},
  {0.723, 0.324, 0.285, 1.0},
  {1.0, 0.0, 0.0, 1.0},
  {0.7, 0.0, 0.0, 1.0},
  {0.15625, 0.78125, 0.2929, 1.0},
  {0.0, 1.0, 0.0, 1.0},
  {0.0, 0.0, 1.0, 1.0},
  {0.0, 0.0, 0.0, 1.0},
  {0.7, 0.7, 0.7, 1.0},
  {0.5, 0.5, 0.5, 1.0},
  {0.2, 0.2, 0.2, 1.0},
  {1.0, 1.0, 0.0, 1.0}};

const GLfloat shiny_high[1]        = {100};  // a value between 1 & 128
const GLfloat shiny_low[1]         = {10};


//---------------------------------------------------------------------------
Draw3d draw3d;
//---------------------------------------------------------------------------


//---------------------------------------------------------------------------
Draw3d::Draw3d()
//---------------------------------------------------------------------------
{
  drawVertices = false;
  drawTriangles = true;
  smoothTriangles = true;
  drawTetras = true;
  drawSurface = false;
  drawOrientations = false;
  drawGroundPlane = true;
  drawStress = false;
  drawForces = false;
  drawForceSum = true;
}


//---------------------------------------------------------------------------
Draw3d::~Draw3d()
//---------------------------------------------------------------------------
{
}


//---------------------------------------------------------------------------
void Draw3d::initOpenGL()
//---------------------------------------------------------------------------
{
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LINE_SMOOTH);

  glClearColor(0.1, 0.1, 0.4, 1.0);
  glShadeModel(GL_SMOOTH);

  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  glMaterialfv(GL_FRONT, GL_SPECULAR, glColor[DRAW_WHITE]);
  glMaterialfv(GL_FRONT, GL_SHININESS, shiny_high);

  glLightfv(GL_LIGHT0, GL_AMBIENT, glColor[DRAW_MEDIUM_GREY]);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, glColor[DRAW_WHITE]);
  glLightfv(GL_LIGHT0, GL_SPECULAR, glColor[DRAW_WHITE]);
  GLfloat light_direction[]    = {10.0, 10.0, 10.0, 0.0};
  glLightfv(GL_LIGHT0, GL_POSITION, light_direction);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, glColor[DRAW_DARK_GREY]);
}


/* ------------------------------------------------------------------------ */
void Draw3d::smoothTriangle(Vector3d &pos1, Vector3d &pos2, Vector3d &pos3,
                        Vector3d &n1, Vector3d &n2, Vector3d &n3, float *color)
/* ------------------------------------------------------------------------ */
{
  glBegin(GL_TRIANGLES);
    glColor4f(color[0], color[1], color[2], color[3]);
    glNormal3f(n1.x, n1.y, n1.z);
    glVertex3f(pos1.x, pos1.y, pos1.z);
    glNormal3f(n2.x, n2.y, n2.z);
    glVertex3f(pos2.x, pos2.y, pos2.z);
    glNormal3f(n3.x, n3.y, n3.z);
    glVertex3f(pos3.x, pos3.y, pos3.z);
  glEnd();
}


/* ------------------------------------------------------------------------ */
void Draw3d::triangle(Vector3d &pos1, Vector3d &pos2, Vector3d &pos3,
                        Vector3d &normal, float *color)
/* ------------------------------------------------------------------------ */
{
  glBegin(GL_TRIANGLES);
    glColor4f(color[0], color[1], color[2], color[3]);
    glNormal3f(normal.x, normal.y, normal.z);
    glVertex3f(pos1.x, pos1.y, pos1.z);
    glVertex3f(pos2.x, pos2.y, pos2.z);
    glVertex3f(pos3.x, pos3.y, pos3.z);
  glEnd();
}


/* ------------------------------------------------------------------------ */
void Draw3d::arrow(Vector3d &pos, Vector3d &n, float r, int color)
/* ------------------------------------------------------------------------ */
{
  Vector3d v1; v1 = pos;
  Vector3d v2; v2 = pos; v2.add(n);

  Vector3d n1, n2, n3;
  n1 = n;
  n1.normalize();
  n2.set(n1.y, n1.z, n1.x);
  n2.cross(n, n2); n2.normalize();
  n3.cross(n, n2); n3.normalize();

  int i,j;
  int num = 6;
  float phi;

  glBegin(GL_QUADS);
  glColor4fv(glColor[color]);

  phi = 0.0;
  for (i = 0; i <= num; i++) {
    Vector3d v, vOld, r1, r2;
    vOld = v;
    r1 = n2; r1.scale(r * cos(phi));
    r2 = n3; r2.scale(r * sin(phi));
    v.add(r1, r2);
    if (i > 0) {
      glNormal3f(v.x, v.y, v.z);
      glVertex3f(v1.x + v.x, v1.y + v.y, v1.z + v.z);
      glNormal3f(vOld.x, vOld.y, vOld.z);
      glVertex3f(v1.x + vOld.x, v1.y + vOld.y, v1.z + vOld.z);
      glNormal3f(vOld.x, vOld.y, vOld.z);
      glVertex3f(v2.x + vOld.x, v2.y + vOld.y, v2.z + vOld.z);
      glNormal3f(v.x, v.y, v.z);
      glVertex3f(v2.x + v.x, v2.y + v.y, v2.z + v.z);
    }
    phi += 2.0 * M_PI / num;
  }
  glEnd();
}


/* ------------------------------------------------------------------------ */
void Draw3d::cube(Vector3d &center, float size, int color)
/* ------------------------------------------------------------------------ */
{
  float cube[6][4][3] = {{{1,1,1},{1,-1,1},{1,-1,-1},{1,1,-1}},
                        {{1,1,1},{1,1,-1},{-1,1,-1},{-1,1,1}},
                        {{-1,-1,-1},{-1,-1,1},{-1,1,1},{-1,1,-1}},
                        {{-1,-1,-1},{1,-1,-1},{1,-1,1},{-1,-1,1}},
                        {{1,1,1},{-1,1,1},{-1,-1,1},{1,-1,1}},
                        {{-1,-1,-1},{-1,1,-1},{1,1,-1},{1,-1,-1}}};
  float normal[6][3] =  {{1,0,0}, {0,1,0}, {-1,0,0},
                         {0,-1,0}, {0,0,1}, {0,0,-1}};
  glBegin(GL_QUADS);
  glColor4fv(glColor[color]);

  for (int i = 0; i < 6; i++) {
    glNormal3f(normal[i][0], normal[i][1], normal[i][2]);
    for (int j = 0; j < 4; j++) {
      glVertex3f(center.x + size * cube[i][j][0],
                 center.y + size * cube[i][j][1],
                 center.z + size * cube[i][j][2]);
    }
  }
  glEnd();
}


/* ------------------------------------------------------------------------ */
void Draw3d::groundPlane(float y)
/* ------------------------------------------------------------------------ */
{
  if (!drawGroundPlane) return;
  glBegin(GL_QUADS);
    glColor4f(0.4, 0.4, 0.4, 1.0);
    glNormal3f(0.0, 1.0, 0.0);
    glVertex3f(-PLANE_SIZE, y, -PLANE_SIZE);
    glVertex3f(-PLANE_SIZE, y,  PLANE_SIZE);
    glVertex3f( PLANE_SIZE, y,  PLANE_SIZE);
    glVertex3f( PLANE_SIZE, y, -PLANE_SIZE);
  glEnd();
}


// -----------------------------------------------------------------
void Draw3d::project(Vector3d &v, int &xi, int &yi, float &depth)
// -----------------------------------------------------------------
{
  // OpenGL view transformation

  GLint viewPort[4];
  GLdouble modelMatrix[16];
  GLdouble projMatrix[16];
  glGetIntegerv(GL_VIEWPORT, viewPort);
  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);

  GLdouble winX, winY, winZ;
  gluProject((GLdouble) v.x, (GLdouble) v.y, (GLdouble) v.z,
      modelMatrix, projMatrix, viewPort, &winX, &winY, &winZ);

  xi = winX; yi = viewPort[3] - winY - 1; depth = winZ;
}

// ---------------------------------------------------------------------
void Draw3d::unProject(int xi, int yi, float depth, Vector3d &v)
// ---------------------------------------------------------------------
{
  // OpenGL view transformation

  GLint viewPort[4];
  GLdouble modelMatrix[16];
  GLdouble projMatrix[16];
  glGetIntegerv(GL_VIEWPORT, viewPort);
  glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);

  yi = viewPort[3] - yi - 1;
  GLdouble wx, wy, wz;
  gluUnProject((GLdouble) xi, (GLdouble) yi, (GLdouble) depth,
    modelMatrix, projMatrix, viewPort, &wx, &wy, &wz);
  v.x = wx; v.y = wy; v.z = wz;
}

