/// @file sparse_matrix.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @author Lin Lin
/// @date 2012-11-28
#ifndef _SPARSE_MATRIX_HPP_
#define _SPARSE_MATRIX_HPP_

#include "environment_impl.hpp"
#include "numvec_impl.hpp"

namespace  PEXSI{

/// @struct SparseMatrix
/// 
/// @brief SparseMatrix describes a sequential sparse matrix saved in
/// compressed sparse column format.
///
/// Note
/// ----
///
/// Since in PEXSI and PPEXSI only symmetric matrix is considered, the
/// compressed sparse row format will also be represented by the
/// compressed sparse column format.
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

/// @struct DistSparseMatrix
///
/// @brief DistSparseMatrix describes a Sparse matrix in the compressed
/// sparse column format (CSC) and distributed with column major partition. 
///
/// Note
/// ----
/// 
/// Since in PEXSI and PPEXSI only symmetric matrix is considered, the
/// compressed sparse row format will also be represented by the
/// compressed sparse column format.
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
