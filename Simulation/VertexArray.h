//---------------------------------------------------------------------------

#ifndef VertexArrayH
#define VertexArrayH
//---------------------------------------------------------------------------

class Vertex;
class Matrix3d;
class Vector3d;

//--------------------------------------------------------------------
class VertexArray {
//--------------------------------------------------------------------
public:
  VertexArray();
  ~VertexArray();

  void   clear() { _numVertices = 0; }
  void   add(Vertex *vertex);
  int    count() { return _numVertices; }

  void   deleteAll();

  Vertex* get(int nr) {
    if (nr < 0 || nr >= _numVertices) return 0;
    return _vertices[nr];
  }

  int find(Vertex *v);
  void remove(Vertex *v);


  Vertex* operator[](int i) { return _vertices[i]; }

  bool contains(Vertex *v) {
    Vertex **vi = _vertices;
    for (int i = 0; i < _numVertices; i++, vi++)
      if (*vi == v) return true;
    return false;
  }

private:
  Vertex **_vertices;
  int   _verticesSize;
  int   _numVertices;
};

#endif
