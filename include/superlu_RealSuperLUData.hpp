#ifndef _SUPERLU_DIST_REAL_SUPERLUDATA_HPP_
#define _SUPERLU_DIST_REAL_SUPERLUDATA_HPP_

//#include "superlu_SuperLUData.hpp"
#include "environment.hpp"

#include "sparse_matrix_impl.hpp"
#include "nummat_impl.hpp"
#include "numvec_impl.hpp"
#include "superlu_SuperLUGrid.hpp"
#include "superlu_SuperLUOptions.hpp"

namespace PEXSI{

  struct SuperNodeType;
  template<typename T> class PMatrix;

class RealSuperLUData_internal;

class RealSuperLUData{//: public SuperLUData{
  protected:
    RealSuperLUData_internal * ptrData;
  public:
    RealSuperLUData( const SuperLUGrid<Real>& g, const SuperLUOptions& opt );
    ~RealSuperLUData();
  
		virtual Int m() const;
		virtual Int n() const;
		virtual void DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Real>& sparseA ); 
		virtual void DestroyAOnly(); 
		virtual void SymbolicFactorize(); 
		virtual void Distribute(); 
		virtual void NumericalFactorize(); 
    virtual void ConvertNRlocToNC	( RealSuperLUData * aptrData );
		virtual void MultiplyGlobalMultiVector( NumMat<Real>& xGlobal, NumMat<Real>& bGlobal ); 
		virtual void DistributeGlobalMultiVector( NumMat<Real>& xGlobal, NumMat<Real>& xLocal ); 
		virtual void GatherDistributedMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& xLocal ); 
		virtual void SolveDistMultiVector( NumMat<Real>& bLocal, DblNumVec& berr ); 
		virtual void CheckErrorDistMultiVector( NumMat<Real>& xLocal, NumMat<Real>& xTrueLocal ); 
		virtual void LUstructToPMatrix( PMatrix<Real>& PMloc ); 
		virtual void SymbolicToSuperNode( SuperNodeType& super );
};


}

#endif
