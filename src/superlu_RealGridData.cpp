#include "superlu_RealGridData.hpp"

#include "superlu_ddefs.h"
#include "Cnames.h"


// SuperLUGrid class
namespace PEXSI{


class RealGridInfo{
  friend class RealGridData;
  protected:
	  gridinfo_t          grid;
};

    RealGridData::RealGridData(){
      info_ = new RealGridInfo;
    }

    RealGridData::~RealGridData(){
      delete info_;
    }



    void RealGridData::GridInit( MPI_Comm comm, Int nprow, Int npcol ){
    	superlu_gridinit(comm, nprow, npcol, &info_->grid);
    }

    void RealGridData::GridExit(  ){
    	superlu_gridexit(&info_->grid);
    }

}
