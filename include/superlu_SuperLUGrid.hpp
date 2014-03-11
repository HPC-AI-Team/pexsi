#ifndef _SUPERLU_DIST_SUPERLUGRID_HPP_
#define _SUPERLU_DIST_SUPERLUGRID_HPP_

#include "environment.hpp"

#include "superlu_RealGridData.hpp"
#include "superlu_ComplexGridData.hpp"


namespace PEXSI{

  /// @class SuperLUGrid
  /// @brief A thin interface for the gridinfo_t structure in SuperLU.
  template< typename T > class SuperLUGrid{
    /// @brief SuperLUMatrix can have access to the grid information.
    public:
    void *        ptrData;
    public:
    SuperLUGrid( MPI_Comm comm, int nprow, int npcol ){};
    ~SuperLUGrid(){};
  };


  template< > class SuperLUGrid<Real>{
    /// @brief SuperLUMatrix can have access to the grid information.
    public:
    RealGridData*        ptrData;
    public:
    SuperLUGrid( MPI_Comm comm, int nprow, int npcol );
    ~SuperLUGrid();
  };


  template< > class SuperLUGrid<Complex>{
    /// @brief SuperLUMatrix can have access to the grid information.
    public:
    ComplexGridData*        ptrData;
    public:
    SuperLUGrid( MPI_Comm comm, int nprow, int npcol );
    ~SuperLUGrid();
  };

}


#include "superlu_SuperLUGrid_impl.hpp"

#endif
