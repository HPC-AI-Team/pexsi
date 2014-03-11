#ifndef _SUPERLU_DIST_COMPLEX_GRIDDATA_HPP_
#define _SUPERLU_DIST_COMPLEX_GRIDDATA_HPP_

#include "superlu_GridData.hpp"

namespace PEXSI{

class ComplexGridInfo;

class ComplexGridData{//: public GridData{
  public:
    ComplexGridInfo * info_;
  public:

    ComplexGridData();
    ~ComplexGridData();
    virtual void GridInit( MPI_Comm comm, Int nprow, Int npcol );
    virtual void GridExit(  );
};


}

#endif
