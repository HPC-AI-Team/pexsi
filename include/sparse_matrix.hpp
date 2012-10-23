#ifndef _SPARSE_MATRIX_HPP_
#define _SPARSE_MATRIX_HPP_

#include "environment_impl.hpp"
#include "numvec_impl.hpp"

namespace  PEXSI{

// *********************************************************************
// Sparse matrix in the compressed sparse column format (CSC)
// *********************************************************************

template <class F> struct SparseMatrix{
	Int          size;                            // Matrix dimension
	Int          nnz;                             // Number of nonzeros
	IntNumVec    colptr;                          // Column index pointer
	IntNumVec    rowind;                          // Starting row index pointer
	NumVec<F>    nzval;                           // Nonzero values for the sparse matrix
};

// Commonly used
typedef NumVec<Real>       DblSparseMatrix;
typedef NumVec<Complex>    CpxSparseMatrix;

} // namespace PEXSI

#endif // _SPARSE_MATRIX_HPP_
