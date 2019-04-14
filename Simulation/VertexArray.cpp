#include "VertexArray.h"

//---------------------------------------------------------------------------


#define INIT_ARRAY_SIZE 1000

//--------------------------------------------------------------------
VertexArray::VertexArray()
//--------------------------------------------------------------------
{
  _vertices = new Vertex*[INIT_ARRAY_SIZE];
  _numVertices = 0;
  _verticesSize = INIT_ARRAY_SIZE;
}

//--------------------------------------------------------------------
VertexArray::~VertexArray()
//--------------------------------------------------------------------
{
  delete[] _vertices;
}


//--------------------------------------------------------------------
void VertexArray::deleteAll()
//--------------------------------------------------------------------
{
  for (int i = 0; i < _numVertices; i++)
    delete _vertices[i];
  _numVertices = 0;
}


//--------------------------------------------------------------------
void VertexArray::add(Vertex *vertex)
//--------------------------------------------------------------------
{
  if (_numVertices >= _verticesSize) {
    _verticesSize *= 2;
    Vertex **newList = new Vertex*[_verticesSize];
    for (int i = 0; i < _numVertices; i++)
      newList[i] = _vertices[i];
    delete[] _vertices;
    _vertices = newList;
  }
  _vertices[_numVertices] = vertex;
  _numVertices++;
}


//--------------------------------------------------------------------
int VertexArray::find(Vertex *v)
//--------------------------------------------------------------------
{
  int i = 0;
  while ((i < _numVertices) && (_vertices[i] != v)) i++;
  if (i >= _numVertices) return -1;
  return i;
}


//--------------------------------------------------------------------
void VertexArray::remove(Vertex *v)
//--------------------------------------------------------------------
{
  int i = find(v); if (i < 0) return;
  _vertices[i] = _vertices[_numVertices-1];
  _numVertices--;
}
