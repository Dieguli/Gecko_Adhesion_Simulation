//---------------------------------------------------------------------------

#ifndef DrawH
#define DrawH

#include "stdio.h"

class Vector3d;

#define DRAW_WHITE 0
#define DRAW_BROWN 1
#define DRAW_PEACH 2
#define DRAW_RED 3
#define DRAW_BRICKRED 4
#define DRAW_JELLO_GREEN 5
#define DRAW_GREEN 6
#define DRAW_BLUE 7
#define DRAW_BLACK 8
#define DRAW_LT_GREY 9
#define DRAW_MEDIUM_GREY 10
#define DRAW_DARK_GREY 11
#define DRAW_YELLOW 12

//---------------------------------------------------------------------------
#define FILENAME_LEN 80

//---------------------------------------------------------------------------

class Draw3d {
public:
  Draw3d();
  ~Draw3d();

  void initOpenGL();

  // ---- draw procs -----------------------------------
  void arrow(Vector3d &pos, Vector3d &n, float r, int color);
  void cube(Vector3d &center, float size, int color);

  void smoothTriangle(Vector3d &pos1, Vector3d &pos2, Vector3d &pos3,
                        Vector3d &n1, Vector3d &n2, Vector3d &n3, float *color);
  void triangle(Vector3d &pos1, Vector3d &pos2, Vector3d &pos3,
                  Vector3d &normal, float *color);

  void groundPlane(float y);

  void project(Vector3d &v, int &xi, int &yi, float &depth);
  void unProject(int xi, int yi, float depth, Vector3d &v);

  // --- parameters ------------------------
  bool drawVertices;
  bool drawTriangles;
  bool smoothTriangles;
  bool drawTetras;
  bool drawSurface;
  bool drawOrientations;
  bool drawGroundPlane;
  bool drawStress;
  bool drawForces;
  bool drawForceSum;
  int  dummy;          // compiler bug: drawForceSum overlaps texturePath!
};

//---------------------------------------------------------------------------

extern Draw3d draw3d;


#endif
