#include "TetraArray.h"
#include <stdio.h>

//---------------------------------------------------------------------------

#define INIT_ARRAY_SIZE 1000
#include "Tetra.h"

//--------------------------------------------------------------------
TetraArray::TetraArray()
//--------------------------------------------------------------------
{
  _tetras = new Tetra*[INIT_ARRAY_SIZE];
  _numTetras = 0;
  _tetrasSize = INIT_ARRAY_SIZE;
  _currentNr = 0;
}

//--------------------------------------------------------------------
TetraArray::TetraArray(int initSize)
//--------------------------------------------------------------------
{
  _tetras = new Tetra*[initSize];
  _numTetras = 0;
  _tetrasSize = initSize;
}

//--------------------------------------------------------------------
TetraArray::~TetraArray()
//--------------------------------------------------------------------
{
  delete[] _tetras;
}


//--------------------------------------------------------------------
void TetraArray::deleteAll()
//--------------------------------------------------------------------
{
  for (int i = 0; i < _numTetras; i++)
    delete _tetras[i];
  _numTetras = 0;
}


//--------------------------------------------------------------------
void TetraArray::add(Tetra *tetra)
//--------------------------------------------------------------------
{
  if (_numTetras >= _tetrasSize) {
    _tetrasSize *= 2;
    Tetra **newList = new Tetra*[_tetrasSize];
    for (int i = 0; i < _numTetras; i++)
      newList[i] = _tetras[i];
    delete[] _tetras;
    _tetras = newList;
  }
  _tetras[_numTetras] = tetra;
  _numTetras++;
}


//--------------------------------------------------------------------
int TetraArray::find(Tetra *tetra)
//--------------------------------------------------------------------
{
  int i = 0;
  while ((i < _numTetras) && (_tetras[i] != tetra)) i++;
  if (i >= _numTetras) return -1;
  return i;
}


//--------------------------------------------------------------------
void TetraArray::rewind()
//--------------------------------------------------------------------
{
  _currentNr = 0;
}


//--------------------------------------------------------------------
Tetra *TetraArray::getNext()
//--------------------------------------------------------------------
{
  while((_currentNr < _numTetras) && _tetras[_currentNr]->_deleted)
    _currentNr++;
  if (_currentNr >= _numTetras) return NULL;
  Tetra *t = _tetras[_currentNr];
  _currentNr++;
  return t;
}


//--------------------------------------------------------------------
void TetraArray::removeAllDeleted()
//--------------------------------------------------------------------
{
  int i = 0;
  while (i < _numTetras) {
    if (_tetras[i]->_deleted) {
      delete _tetras[i];
      _tetras[i] = _tetras[_numTetras-1];
      _numTetras--;
    }
    else i++;
  }
}

//--------------------------------------------------------------------
void TetraArray::remove(int i)
//--------------------------------------------------------------------
{
  if (i < 0 || i >= _numTetras) return;
  _tetras[i] = _tetras[_numTetras-1];
  _numTetras--;
}
