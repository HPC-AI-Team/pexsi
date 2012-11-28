/// @file sparse_matrix_decl.hpp
/// @brief Sparse matrix and Distributed sparse matrix in compressed
/// column format.
/// @author Lin Lin
/// @date 2012-11-10
#ifndef _SPARSE_MATRIX_DECL_HPP_
#define _SPARSE_MATRIX_DECL_HPP_

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
	/// @brief Matrix dimension.
	Int          size;         

	/// @brief Total number of nonzeros elements.
	Int          nnz;                             

	/// @brief Local number of local nonzeros elements on this processor.
	Int          nnzLocal;                        

	/// @brief Dimension numColLocal + 1, storing the pointers to the
	/// nonzero row indices and nonzero values in rowptrLocal and
	/// nzvalLocal, respectively.  numColLocal is the number
	/// of local columns saved on this processor. The indices are 1-based
	/// (FORTRAN-convention), i.e.  colptrLocal[0] = 1. 
	IntNumVec    colptrLocal;                     

	/// @brief Dimension nnzLocal, storing the nonzero indices.
	/// The indices are 1-based (FORTRAN-convention), i.e. the first row
	/// index is 1. 
	IntNumVec    rowindLocal;                    
	
	/// @brief Dimension nnzLocal, storing the nonzero values.
	NumVec<F>    nzvalLocal;                      

	/// @brief MPI communicator
	MPI_Comm     comm;        
};

// Commonly used
typedef DistSparseMatrix<Real>       DblDistSparseMatrix;
typedef DistSparseMatrix<Complex>    CpxDistSparseMatrix;



} // namespace PEXSI

#endif // _SPARSE_MATRIX_DECL_HPP_
