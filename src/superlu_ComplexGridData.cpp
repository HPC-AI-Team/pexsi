#include "superlu_ComplexGridData.hpp"

#include "superlu_zdefs.h"
#include "Cnames.h"
// SuperLUGrid class
namespace PEXSI{


class ComplexGridInfo{
  friend class ComplexGridData;
  protected:
	  gridinfo_t          grid;
};

    ComplexGridData::ComplexGridData(){
      info_ = new ComplexGridInfo;
    }

    ComplexGridData::~ComplexGridData(){
      delete info_;
    }

    void ComplexGridData::GridInit( MPI_Comm comm, Int nprow, Int npcol ){
    	superlu_gridinit(comm, nprow, npcol, &info_->grid);
    }

    void ComplexGridData::GridExit(  ){
    	superlu_gridexit(&info_->grid);
    }

}
