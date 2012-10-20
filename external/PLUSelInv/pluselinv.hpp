#ifndef _PLUSELINV_HPP
#define _PLUSELINV_HPP

/**********************************************************************
 * LLIN: Common utilities
/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "comobject.hpp"
#include "vec2t.hpp"
#include "vec3t.hpp"
#include "numtns.hpp"
#include "serialize.hpp"

using std::vector;
using std::pair;
using std::map;
using std::set;
using std::cerr;
using std::cout;
using std::ostream;
using std::istream;
using std::istringstream;
using std::ostringstream;
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::setw;
using std::scientific;

/**********************************************************************
 * LLIN: Output LU structure
/**********************************************************************/


class LBlock{
  public:
    int _blkind;
    int _nrows;
    int _ncols;
    IntNumVec _rows;
    DblNumMat _nzval;
  public:
    int blkind() {return _blkind;}
    int nrows()  {return _nrows;}
    int ncols()  { return _ncols;}
    IntNumVec& rows() {return _rows;}
    DblNumMat& nzval() {return _nzval;}
    int rows(int i) {return _rows[i];}
    double nzval(int i, int j) {return _nzval(i,j);}

    LBlock() {_blkind=-1;_nrows=0;_ncols=0;}
    ~LBlock() {}
    LBlock& operator= (const LBlock& LB) {
      _blkind = LB._blkind;
      _nrows  = LB._nrows;
      _ncols  = LB._ncols;
      _rows   = LB._rows;
      _nzval  = LB._nzval;
      return *this;
    }

};

class UBlock{
  public:
    int _blkind;
    int _nrows;
    int _ncols;
    IntNumVec _cols;
    DblNumMat _nzval;
  public:
    int blkind() {return _blkind;}
    int nrows()  {return _nrows;}
    int ncols()  { return _ncols;}
    IntNumVec& cols() {return _cols;}
    DblNumMat& nzval() {return _nzval;}
    int cols(int i) {return _cols[i];}
    double nzval(int i, int j) {return _nzval(i,j);}

    UBlock() {_blkind=-1;_nrows=0;_ncols=0;}
    ~UBlock() {}
    
    UBlock& operator= (const UBlock& UB) {
      _blkind = UB._blkind;
      _nrows  = UB._nrows;
      _ncols  = UB._ncols;
      _cols   = UB._cols;
      _nzval  = UB._nzval;
      return *this;
    }
};

class PMatrix{
  public:
    int _ndof;
    int _nsupers;
    int _nbc;
    int _nbr;
    IntNumVec _xsup;
    IntNumVec _supno;
    IntNumVec _nblkl;
    IntNumVec _nblku;
    vector<vector<LBlock> > _L;
    vector<vector<UBlock> > _U;
  public:
    PMatrix(){;}
    ~PMatrix(){;}
    int ndof(){return _ndof;}
    int nsupers(){return _nsupers;}
    int nbc(){return _nbc;}
    int nbr(){return _nbr;}
    int xsup(int i){return _xsup[i];}
    int supsize(int i){return _xsup[i+1]-_xsup[i];}
    int superno(int i){return _supno[i];}
    int nblkl(int i){return _nblkl[i];}
    int nblku(int i){return _nblku[i];}
    vector<vector<LBlock> >& L(){return _L;}
    vector<LBlock>& L(int i){return _L[i];}
    vector<vector<UBlock> >& U(){return _U;}
    vector<UBlock>& U(int i){return _U[i];}

    int DumpL(string filename, gridinfo_t* grid);
    int DumpU(string filename, gridinfo_t* grid);
    
    int DumpLBlock(int isup, int jsup, string filename, gridinfo_t* grid);
    int DumpUBlock(int isup, int jsup, string filename, gridinfo_t* grid);

    int CondDiagBlock(gridinfo_t* grid);
    int DumpDiagVec(string filename, gridinfo_t *grid);

};

class BlockPtn
{
public:
  IntNumVec _ownerinfo;
public:
  BlockPtn() {;}
  ~BlockPtn() {;}
  IntNumVec& ownerinfo() { return _ownerinfo; }
  int owner(int key) {
    return _ownerinfo(key);
  }
};


#define KITEMS  5
#define ISOURCE 0
#define ITARGET 1
#define ISUPER  2
#define IROWIND 3
#define NROWS   4

int SuperLU2SelInv(int n, LUstruct_t *LUstruct, gridinfo_t *grid,
                   PMatrix& PMloc);

int DiagTri(PMatrix& PMloc, int ksup, gridinfo_t *grid);
int ScaleLinv(PMatrix& PMloc, int ksup, gridinfo_t *grid);
int ScaleUinv(PMatrix& PMloc, int ksup, gridinfo_t *grid);
int DiagInnerProd(PMatrix& PMloc, vector<vector<int> >& localEtree, int ksup, gridinfo_t *grid);
int SinvPL(PMatrix& PMloc, vector<vector<int> >& localEtree, int ksup, gridinfo_t* grid);
int SinvPU(PMatrix& PMloc, vector<vector<int> >& localEtree, int ksup, gridinfo_t* grid);

int PLUSelInv(gridinfo_t *grid, PMatrix& PMloc, vector<vector<int> >& localEtree);
int ConstructLocalEtree(int n, gridinfo_t *grid, PMatrix& PMloc, 
			vector<vector<int> >& localEtree);

int BcastLBlock(LBlock& LB, int mykey, int srckey, MPI_Comm comm);
int BcastUBlock(UBlock& UB, int mykey, int srckey, MPI_Comm comm);

int DumpLocalEtree(vector<vector<int> >& localEtree, gridinfo_t* grid);

//void blockinvert(Lstruct *PMlocL, LUstruct_t *LUstruct, gridinfo_t *grid,
//                 int jsup, int n);
//void ScalePMbyL(Lstruct *PMlocL, gridinfo_t *grid, int ksup, int nsupers);
//void SetPL(Lstruct *PMlocL, gridinfo_t *grid, PLstruct *PL, 
//           int ksup, int nsupers);

#define LBlock_Number 5
enum {
  LBlock_blkind   = 0,
  LBlock_nrows    = 1,
  LBlock_ncols    = 2,
  LBlock_rows     = 3,
  LBlock_nzval    = 4,
};

int inline serialize(LBlock& val, ostream& os, const vector<int>& mask){
  int i = 0;
  if(mask[i]==1) serialize(val._blkind, os, mask); i++;
  if(mask[i]==1) serialize(val._nrows,  os, mask); i++;
  if(mask[i]==1) serialize(val._ncols,  os, mask); i++;
  if(mask[i]==1) serialize(val._rows, os, mask);   i++;
  if(mask[i]==1) serialize(val._nzval, os, mask);  i++;
  iA( i == LBlock_Number );
  return 0;
}

int inline deserialize(LBlock& val, istream& is, const vector<int>& mask){
  int i = 0;
  if(mask[i]==1) deserialize(val._blkind, is, mask); i++;
  if(mask[i]==1) deserialize(val._nrows,  is, mask); i++;
  if(mask[i]==1) deserialize(val._ncols,  is, mask); i++;
  if(mask[i]==1) deserialize(val._rows,   is, mask); i++;
  if(mask[i]==1) deserialize(val._nzval,  is, mask); i++; 
  iA( i == LBlock_Number );
  return 0;
}

#define UBlock_Number 5
enum {
  UBlock_blkind   = 0,
  UBlock_nrows    = 1,
  UBlock_ncols    = 2,
  UBlock_cols     = 3,
  UBlock_nzval    = 4,
};



int inline serialize(UBlock& val, ostream& os, const vector<int>& mask){
  int i = 0;
  if(mask[i]==1) serialize(val._blkind, os, mask); i++;
  if(mask[i]==1) serialize(val._nrows,  os, mask); i++;
  if(mask[i]==1) serialize(val._ncols,  os, mask); i++;
  if(mask[i]==1) serialize(val._cols, os, mask);   i++;
  if(mask[i]==1) serialize(val._nzval, os, mask);  i++;
  iA( i == UBlock_Number );
  return 0;
}

int inline deserialize(UBlock& val, istream& is, const vector<int>& mask){
  int i = 0;
  if(mask[i]==1) deserialize(val._blkind, is, mask); i++;
  if(mask[i]==1) deserialize(val._nrows,  is, mask); i++;
  if(mask[i]==1) deserialize(val._ncols,  is, mask); i++;
  if(mask[i]==1) deserialize(val._cols,   is, mask); i++;
  if(mask[i]==1) deserialize(val._nzval,  is, mask); i++; 
  iA( i == UBlock_Number );
  return 0;
}

inline ostream& operator<<( ostream& os, LBlock& LB)
{
  os<<"LBlock"<< endl;
  os<<"blkind = "<< LB.blkind() << endl;
  os<<"nrows = "<< LB.nrows() << endl;
  os<<"ncols = "<< LB.ncols() << endl;
  os<<"rows = "<< LB.rows() << endl;
  os<<"nzval = "<< LB.nzval() << endl;
  return os;
}

inline ostream& operator<<( ostream& os, UBlock& UB)
{
  os<<"UBlock"<< endl;
  os<<"blkind = "<< UB.blkind() << endl;
  os<<"nrows = "<< UB.nrows() << endl;
  os<<"ncols = "<< UB.ncols() << endl;
  os<<"colss = "<< UB.cols() << endl;
  os<<"nzval = "<< UB.nzval() << endl;
  return os;
}


#endif
