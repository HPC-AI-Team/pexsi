#ifndef _SUPERLU_DIST_REAL_GRIDDATA_HPP_
#define _SUPERLU_DIST_REAL_GRIDDATA_HPP_

#include "superlu_GridData.hpp"

namespace PEXSI{

class RealGridInfo;

class RealGridData{//: public GridData{
  public:
    RealGridInfo * info_;
  public:

    RealGridData();
    ~RealGridData();
    virtual void GridInit( MPI_Comm comm, Int nprow, Int npcol );
    virtual void GridExit(  );
};


}

#endif
