#ifndef _SUPERLU_DIST_INTERF_IMPL_HPP_
#define _SUPERLU_DIST_INTERF_IMPL_HPP_


namespace PEXSI{

inline SuperLUMatrix<Real>::SuperLUMatrix	( const SuperLUGrid<Real>& g, const SuperLUOptions& opt )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::SuperLUMatrix");
#endif
  ptrData = new RealSuperLUData(g,opt);
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Real>::SuperLUMatrix  ----- 

inline SuperLUMatrix<Real>::~SuperLUMatrix	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::~SuperLUMatrix");
#endif
  delete ptrData;
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Real>::~SuperLUMatrix  ----- 


inline Int SuperLUMatrix<Real>::m (  ) const	
{
	return ptrData->m();
} 		// -----  end of method SuperLUMatrix<Real>::m  ----- 



inline Int SuperLUMatrix<Real>::n (  ) const	
{
	return ptrData->n();
} 		// -----  end of method SuperLUMatrix<Real>::n  ----- 

inline void
SuperLUMatrix<Real>::DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Real>& sparseA )
{
#ifndef _RELEASE_
	PushCallStack( "SuperLUMatrix<Real>::DistSparseMatrixToSuperMatrixNRloc" );
#endif
  ptrData->DistSparseMatrixToSuperMatrixNRloc(sparseA );
#ifndef _RELEASE_
	PopCallStack();
#endif
	return;
} 		// -----  end of method SuperLUMatrix<Real>::DistSparseMatrixToSuperMatrixNRloc ----- 


inline void
SuperLUMatrix<Real>::DestroyAOnly	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::DestroyAOnly");
#endif
  ptrData->DestroyAOnly();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::DestroyAOnly  ----- 

inline void
SuperLUMatrix<Real>::SymbolicFactorize	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::SymbolicFactorize");
#endif
  ptrData->SymbolicFactorize();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::SymbolicFactorize  ----- 


inline void
SuperLUMatrix<Real>::Distribute	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::Distribute");
#endif
  ptrData->Distribute();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::Distribute  ----- 


inline void
SuperLUMatrix<Real>::NumericalFactorize	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::NumericalFactorize");
#endif
  ptrData->NumericalFactorize();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::NumericalFactorize  ----- 


inline void
SuperLUMatrix<Real>::ConvertNRlocToNC	( SuperLUMatrix& AGlobal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::ConvertNRlocToNC");
#endif
  ptrData->ConvertNRlocToNC(AGlobal.ptrData);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::ConvertNRlocToNC  ----- 

inline void
SuperLUMatrix<Real>::MultiplyGlobalMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& bGlobal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::MultiplyGlobalMultiVector");
#endif
  ptrData->MultiplyGlobalMultiVector(xGlobal, bGlobal);
#ifndef _RELEASE_
	PopCallStack();
#endif
 
	return ;
} 		// -----  end of method SuperLUMatrix<Real>::MultiplyGlobalMultiVector  ----- 


inline void
SuperLUMatrix<Real>::DistributeGlobalMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& xLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::DistributeGlobalMultiVector");
#endif
  ptrData->DistributeGlobalMultiVector(xGlobal, xLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::DistributeGlobalMultiVector  ----- 


inline void SuperLUMatrix<Real>::GatherDistributedMultiVector	( NumMat<Real>& xGlobal, NumMat<Real>& xLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::GatherDistributedMultiVector");
#endif
  ptrData->GatherDistributedMultiVector(xGlobal, xLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::GatherDistributedMultiVector  ----- 


inline void
SuperLUMatrix<Real>::SolveDistMultiVector	( NumMat<Real>& bLocal, DblNumVec& berr )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::SolveDistMultiVector");
#endif
  ptrData->SolveDistMultiVector(bLocal, berr );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::SolveDistMultiVector  ----- 


inline void
SuperLUMatrix<Real>::CheckErrorDistMultiVector	( NumMat<Real>& xLocal, NumMat<Real>& xTrueLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::CheckErrorDistMultiVector");
#endif
  ptrData->CheckErrorDistMultiVector(xLocal, xTrueLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::CheckErrorDistMultiVector  ----- 


inline void
SuperLUMatrix<Real>::LUstructToPMatrix	( PMatrix<Real>& PMloc )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::LUstructToPMatrix");
#endif
  ptrData->LUstructToPMatrix(PMloc);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::LUstructToPMatrix  ----- 



inline void
SuperLUMatrix<Real>::SymbolicToSuperNode	( SuperNodeType& super )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Real>::SymbolicToSuperNode");
#endif
  ptrData->SymbolicToSuperNode(super);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Real>::SymbolicToSuperNode  ----- 

}

namespace PEXSI{

inline SuperLUMatrix<Complex>::SuperLUMatrix	( const SuperLUGrid<Complex>& g, const SuperLUOptions& opt )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::SuperLUMatrix");
#endif
  ptrData = new ComplexSuperLUData(g,opt);
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Complex>::SuperLUMatrix  ----- 

inline SuperLUMatrix<Complex>::~SuperLUMatrix	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::~SuperLUMatrix");
#endif
  delete ptrData;
#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix<Complex>::~SuperLUMatrix  ----- 

inline Int SuperLUMatrix<Complex>::m (  ) const	
{
return ptrData->m();
}		// -----  end of method SuperLUMatrix<Complex>::m  ----- 

inline Int SuperLUMatrix<Complex>::n (  ) const	
{
	return ptrData->n();
} 		// -----  end of method SuperLUMatrix<Complex>::n  ----- 

inline void
SuperLUMatrix<Complex>::DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Complex>& sparseA )
{
#ifndef _RELEASE_
	PushCallStack( "SuperLUMatrix<Complex>::DistSparseMatrixToSuperMatrixNRloc" );
#endif
  ptrData->DistSparseMatrixToSuperMatrixNRloc(sparseA );
#ifndef _RELEASE_
	PopCallStack();
#endif
	return;
} 		// -----  end of method SuperLUMatrix<Complex>::DistSparseMatrixToSuperMatrixNRloc ----- 

inline void
SuperLUMatrix<Complex>::DestroyAOnly	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::DestroyAOnly");
#endif
  ptrData->DestroyAOnly();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::DestroyAOnly  ----- 

inline void
SuperLUMatrix<Complex>::SymbolicFactorize	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::SymbolicFactorize");
#endif
  ptrData->SymbolicFactorize();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::SymbolicFactorize  ----- 

inline void
SuperLUMatrix<Complex>::Distribute	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::Distribute");
#endif
  ptrData->Distribute();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::Distribute  ----- 

inline void
SuperLUMatrix<Complex>::NumericalFactorize	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::NumericalFactorize");
#endif
  ptrData->NumericalFactorize();
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::NumericalFactorize  ----- 

inline void
SuperLUMatrix<Complex>::ConvertNRlocToNC	( SuperLUMatrix& AGlobal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::ConvertNRlocToNC");
#endif
  ptrData->ConvertNRlocToNC(AGlobal.ptrData);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::ConvertNRlocToNC  ----- 

inline void
SuperLUMatrix<Complex>::MultiplyGlobalMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& bGlobal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::MultiplyGlobalMultiVector");
#endif
  ptrData->MultiplyGlobalMultiVector(xGlobal, bGlobal);
#ifndef _RELEASE_
	PopCallStack();
#endif
 
	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::MultiplyGlobalMultiVector  ----- 

inline void
SuperLUMatrix<Complex>::DistributeGlobalMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::DistributeGlobalMultiVector");
#endif
  ptrData->DistributeGlobalMultiVector(xGlobal, xLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::DistributeGlobalMultiVector  ----- 

inline void SuperLUMatrix<Complex>::GatherDistributedMultiVector	( NumMat<Complex>& xGlobal, NumMat<Complex>& xLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::GatherDistributedMultiVector");
#endif
  ptrData->GatherDistributedMultiVector(xGlobal, xLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::GatherDistributedMultiVector  ----- 

inline void
SuperLUMatrix<Complex>::SolveDistMultiVector	( NumMat<Complex>& bLocal, DblNumVec& berr )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::SolveDistMultiVector");
#endif
  ptrData->SolveDistMultiVector(bLocal, berr );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::SolveDistMultiVector  ----- 

inline void
SuperLUMatrix<Complex>::CheckErrorDistMultiVector	( NumMat<Complex>& xLocal, NumMat<Complex>& xTrueLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::CheckErrorDistMultiVector");
#endif
  ptrData->CheckErrorDistMultiVector(xLocal, xTrueLocal );
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::CheckErrorDistMultiVector  ----- 

inline void
SuperLUMatrix<Complex>::LUstructToPMatrix	( PMatrix<Complex>& PMloc )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::LUstructToPMatrix");
#endif
  ptrData->LUstructToPMatrix(PMloc);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::LUstructToPMatrix  ----- 

inline void
SuperLUMatrix<Complex>::SymbolicToSuperNode	( SuperNodeType& super )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix<Complex>::SymbolicToSuperNode");
#endif
  ptrData->SymbolicToSuperNode(super);
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix<Complex>::SymbolicToSuperNode  ----- 

}





#endif
