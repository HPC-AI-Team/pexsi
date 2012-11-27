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
typedef SparseMatrix<Real>       DblSparseMatrix;
typedef SparseMatrix<Complex>    CpxSparseMatrix;

// *********************************************************************
// Sparse matrix in the compressed sparse column format (CSC) and
// distributed with column major partition. (This can be trivially
// transformed to the row major partition for SYMMETRIC matrices.
//
// TODO add class Type so that Type can be NR (row major distributed)
// and NCloc (column major distributed)
// TODO Comment on the size of indices
// *********************************************************************

template <class F> struct DistSparseMatrix{
	Int          size;                            // Matrix dimension
	Int          nnz;                             // Number of nonzeros
	Int          nnzLocal;                        // Number of local nonzeros
	IntNumVec    colptrLocal;                     // Local column index pointer
	IntNumVec    rowindLocal;                     // Local starting row index pointer
	NumVec<F>    nzvalLocal;                      // Local nonzero values for the sparse matrix
	MPI_Comm     comm;                            // MPI Communicator
};

// Commonly used
typedef DistSparseMatrix<Real>       DblDistSparseMatrix;
typedef DistSparseMatrix<Complex>    CpxDistSparseMatrix;



} // namespace PEXSI

#endif // _SPARSE_MATRIX_HPP_
