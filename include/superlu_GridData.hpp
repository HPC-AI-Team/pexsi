#ifndef _SUPERLU_DIST_GRIDDATA_HPP_
#define _SUPERLU_DIST_GRIDDATA_HPP_

#include "environment.hpp"

namespace PEXSI{

class GridData{
  public:
    virtual void GridInit( MPI_Comm comm, Int nprow, Int npcol ) =0;
    virtual void GridExit(  ) =0;
};


}

#endif
