//---------------------------------------------------------------------------

#ifndef WorldH
#define WorldH
//---------------------------------------------------------------------------

#include <stdio.h>
#include "Utils.h"

class TetraMesh;
class Vector3d;
class Bounds3d;

//---------------------------------------------------------------------------
class World
//--------------------------------------------------------------------
{
public:
  World(char *applicationPath);
  ~World();

  void add(TetraMesh *mesh);
  void animate();
  void draw();

  void mouseDown(int x, int y);
  void mouseMove(int x, int y);
  void mouseUp(int x, int y);
  void getBounds(Bounds3d &bounds);
  int numTetras();
  int numVertices();
  int numTetraMeshes() { return _numTetraMeshes; }

  //Interacción entre vértices
  void interact (TetraMesh *mesh1, TetraMesh *mesh2);
  TetraMesh **_tetraMeshes;
  int _numTetraMeshes;
  int _tetraMeshesSize;

private:
  void setupTestWorld();
  float _groundPlaneY;
  float _groundPlaneTexNr;

  //int _numTetraMeshes;
  //int _tetraMeshesSize;
  //TetraMesh **_tetraMeshes;
};


#endif
