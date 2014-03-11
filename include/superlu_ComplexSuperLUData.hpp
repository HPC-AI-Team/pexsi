#ifndef _SUPERLU_DIST_COMPLEX_SUPERLUDATA_HPP_
#define _SUPERLU_DIST_COMPLEX_SUPERLUDATA_HPP_

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

class ComplexSuperLUData_internal;

class ComplexSuperLUData{//: public SuperLUData{
  protected:
    ComplexSuperLUData_internal * ptrData;
  public:
    ComplexSuperLUData( const SuperLUGrid<Complex>& g, const SuperLUOptions& opt );
    ~ComplexSuperLUData();

		virtual Int m() const;
		virtual Int n() const;
		virtual void DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Complex>& sparseA ); 
		virtual void DestroyAOnly(); 
		virtual void SymbolicFactorize(); 
		virtual void Distribute(); 
		virtual void NumericalFactorize(); 
    virtual void ConvertNRlocToNC	( ComplexSuperLUData * aptrData );
		virtual void MultiplyGlobalMultiVector( NumMat<Complex>& xGlobal, NumMat<Complex>& bGlobal ); 
		virtual void DistributeGlobalMultiVector( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal ); 
		virtual void GatherDistributedMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal ); 
		virtual void SolveDistMultiVector( NumMat<Complex>& bLocal, DblNumVec& berr ); 
		virtual void CheckErrorDistMultiVector( NumMat<Complex>& xLocal, NumMat<Complex>& xTrueLocal ); 
		virtual void LUstructToPMatrix( PMatrix<Complex>& PMloc ); 
		virtual void SymbolicToSuperNode( SuperNodeType& super );
};


}

#endif
