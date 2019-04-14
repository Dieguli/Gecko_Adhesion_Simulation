//---------------------------------------------------------------------------

#ifndef TetraArrayH
#define TetraArrayH
//---------------------------------------------------------------------------

class Tetra;

//--------------------------------------------------------------------
class TetraArray {
//--------------------------------------------------------------------
public:
  TetraArray();
  TetraArray(int initSize);
  ~TetraArray();

  void   deleteAll();

  void   clear() { _numTetras = 0; }
  void   add(Tetra *tetra);
  int    count() { return _numTetras; }
  int    find(Tetra *tetra);
  void   removeAllDeleted();
  void   remove(int i);

  Tetra* get(int nr) {
    if (nr < 0 || nr >= _numTetras) return 0;
    return _tetras[nr];
  }

  void   rewind();
  Tetra *getNext();

  Tetra* operator[](int i) { return _tetras[i]; }

private:
  Tetra **_tetras;
  int   _tetrasSize;
  int   _numTetras;
  int   _currentNr;
};


#endif
