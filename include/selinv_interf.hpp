/// @file selinv_interf.hpp
/// @brief Interface with the sequential version of SelInv.
/// @author Lin Lin
/// @date 2012-11-28
#ifndef _SELINV_INTERF_HPP_
#define _SELINV_INTERF_HPP_

#include "sparse_matrix_impl.hpp"

namespace PEXSI{

class SelInvInterface{
private:
	bool   isSelInvInitialized_;
	Int    token_;

public:

	SelInvInterface( );

	~SelInvInterface();

	void SymbolicFactorize( SparseMatrix<Scalar>& A, Int order, Int* perm, Int& Lnnz );

	void NumericalFactorize( SparseMatrix<Scalar>& A );

	void Solve( Scalar* x, Scalar* rhs );

	void SelInv( SparseMatrix<Scalar>& Ainv );
};

}

#endif  // _SELINV_INTERF_HPP_
