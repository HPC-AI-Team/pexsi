/*
	 Copyright (c) 2012 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Author: Lin Lin
	 
   This file is part of PEXSI. All rights reserved.

	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions are met:

	 (1) Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.
	 (3) Neither the name of the University of California, Lawrence Berkeley
	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
	 be used to endorse or promote products derived from this software without
	 specific prior written permission.

	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	 You are under no obligation whatsoever to provide any bug fixes, patches, or
	 upgrades to the features, functionality or performance of the source code
	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
	 available either publicly, or directly to Lawrence Berkeley National
	 Laboratory, without imposing a separate written license agreement for such
	 Enhancements, then you hereby grant the following license: a non-exclusive,
	 royalty-free perpetual license to install, use, modify, prepare derivative
	 works, incorporate into other computer software, distribute, and sublicense
	 such enhancements or derivative works thereof, in binary and source code form.
*/
/// @file superlu_dist_interf.hpp
/// @brief Inteface with SuperLU_Dist (version 3.0 and later)
/// @date 2012-11-12
#ifndef _SUPERLU_DIST_INTERF_HPP_
#define _SUPERLU_DIST_INTERF_HPP_

// Interface with PSelInv
#include "pselinv.hpp"

// Interface with sparse matrix (CSC format)
#include  "sparse_matrix_impl.hpp"

// Interface with LAPACK
#include  "lapack.hpp"

namespace PEXSI{
  struct SuperNode;
  class PMatrix;

/// @class SuperLUGrid
/// @brief A thin interface for the gridinfo_t strucutre in SuperLU.
class SuperLUGrid{
	/// @brief SuperLUMatrix can have access to the grid information.
	friend class SuperLUMatrix;
private:
	struct           GridData;
	GridData*        ptrData;
public:
	SuperLUGrid( MPI_Comm comm, int nprow, int npcol );

	~SuperLUGrid();

};


/// @struct SuperLUOptions
/// @brief A thin interface for setting the SuperLU options from
/// outside.  
///
struct SuperLUOptions{
	Int              numProcSymbFact;
  Int              maxPipelineDepth; 

	std::string      ColPerm;

	// Member functions to setup the default value
	SuperLUOptions(): numProcSymbFact(0), maxPipelineDepth(-1), ColPerm("MMD_AT_PLUS_A") {}
};

/// @class SuperLUMatrix
/// @brief An thin interface to keep the main code insulated from
/// the source code of SuperLU. 
///
///
/// Procedure for factorization
/// ---------------------------
///
/// - Construct a distributed compressed sparse column format matrix
/// (complex arithmetic for example).
///  
///   \code{.cpp}
///   DistSparseMatrix<Complex>  AMat;
///   ...(Construct AMat)
///   \endcode
///
/// - Create a SuperMatrix from a DistSparseMatrix.
///   
///   \code{.cpp}
///   SuperLUMatrix luMat( g );
///   luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
///   \endcode
///
/// - Symbolic factorization.
///
///   \code{.cpp}
///   luMat.SymbolicFactorize();
///   luMat.DestroyAOnly();
///		\endcode
///
///		<b> NOTE: </b> Destroy the SuperMatrix by DestroyAOnly is
///		important if the symbolic factorization is to be reused for
///		factorizing matrices of the same pattern.
///
/// - Redistribute and factorize a matrix of the same sparsity pattern
/// as AMat.
///
///   \code{.cpp}
///   DistSparseMatrix<Complex>  BMat;
///   ...(Construct BMat)
///   luMat.DistSparseMatrixToSuperMatrixNRloc( BMat );
///   luMat.Distribute();
///		\endcode
///
/// - Numerical factorization
///
///   \code{.cpp}
///   luMat.NumericalFactorize();
///   \endcode
///
/// - (Optional) Solve multivectors (for consistency check).
///
///   Construct a global matrix.
///
///   \code{.cpp}
///   SuperLUMatrix A1( g ), GA( g );
///   A1.DistSparseMatrixToSuperMatrixNRloc( AMat );
///   A1.ConvertNRlocToNC( GA );
///   \endcode
///
///   Construct the distributed right hand sides and the exact
///   solution.
///
///   \code{.cpp}
///   CpxNumMat xTrueGlobal(n, nrhs), bGlobal(n, nrhs);
///   CpxNumMat xTrueLocal, bLocal;
///   UniformRandom( xTrueGlobal );
///   GA.MultiplyGlobalMultiVector( xTrueGlobal, bGlobal );
///   A1.DistributeGlobalMultiVector( xTrueGlobal, xTrueLocal );
///   A1.DistributeGlobalMultiVector( bGlobal,     bLocal );
///   \endcode
///   
///   Solve and check the error.
///
///   \code{.cpp}
///   luMat.SolveDistMultiVector( bLocal, berr );
///   luMat.CheckErrorDistMultiVector( bLocal, xTrueLocal );
///   \endcode
///
///
/// Note
/// ----
///
/// - The memory allocation etc does not need to be managed explicitly
///	using the SuperLUMatrix class.
///
/// - For parallel selected inversion, see PMatrix.
///
class SuperLUMatrix{
private:
	/// @struct SuperLUData
	///
	/// @brief  Data of a matrix in the SuperLU format.
	///
	/// SuperLUData is only used to define ptrData which is a private member
	/// SuperLUMatrix, and is only defined in superlu_dist_interf.cpp.
	///
	struct        SuperLUData;
	/// @brief Actual pointer to save the data related to SuperLU.
	SuperLUData*  ptrData;
public:
	
	SuperLUMatrix( const SuperLUGrid& g, const SuperLUOptions& opt = SuperLUOptions() );
	~SuperLUMatrix();

	Int m() const;

	Int n() const;

	/// @brief DistSparseMatrixToSuperMatrixNRloc converts a distributed
	/// sparse matrix in compressed sparse column format into the SuperLU
	/// compressed row format. 
	///
	/// TODO real arithmetic
	/// FIXME: IntNumVec convention.  Assumes a symmetric matrix
	void DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Scalar>& sparseA );

	/// @brief DestroyAOnly releases the data in A but keeps other data,
	/// such as LUstruct. This allows one to perform factorization of
	/// matrices of the same pattern, such as the option
	/// fact = SamePattern_SameRowPerm in SuperLU_DIST.
	void DestroyAOnly();

	/// @brief SymbolicFactorize factorizes the superlu matrix symbolically.
	void SymbolicFactorize();

	/// @brief Distribute redistrbutes the SuperMatrix in parallel so that
	/// it is ready for the numerical factorization.
	///
	/// TODO Real arithmetic
	void Distribute();

	/// @brief NumericalFactorize performs LU factorization numerically.  
	///
	/// The matrix should have been permuted and distributed.  
	///
	/// TODO Real arithmetic. 
	void NumericalFactorize();

	/// @brief ConvertNRlocToNC converts a distributed compressed sparse
	/// row matrix to a global compressed sparse column matrix.
	///
	/// @param[out] AGlobal
	void ConvertNRlocToNC( SuperLUMatrix& AGlobal );

	/// @brief MultiplyGlobalMultiVector computes b = A * x.
	///
	/// Both xGlobal and bGlobal are saved globally.
	///
	/// @param[in] xGlobal
	/// @param[out] bGlobal
	void MultiplyGlobalMultiVector( NumMat<Scalar>& xGlobal, NumMat<Scalar>& bGlobal );

	/// @brief DistributeGlobalMultiVector distributes a global multivector into a
	/// local multivector according to the compressed row format of A.
	///
	/// @param[in] xGlobal
	/// @param[out] xLocal
	void DistributeGlobalMultiVector( NumMat<Scalar>& xGlobal, NumMat<Scalar>& xLocal );

  void GatherDistributedMultiVector	( NumMat<Scalar>& xGlobal, NumMat<Scalar>& xLocal );

	/// @brief SolveDistMultiVector A x = b with b overwritten by x for
	/// distributed multivector.
	///
	/// @param[in,out] bLocal Right hand side savd in the distributed format.
	/// @param[out] berr The componentwise relative backward error of each solution vector.
	void SolveDistMultiVector( NumMat<Scalar>& bLocal, DblNumVec& berr );


	/// @brief CheckErrorDistMultiVector prints out the error by direct
	/// comparison with the true solution in distributed format.
	///
	/// @param[in] xLocal The computed solution.
	/// @param[in] xTrueLocal The true solution.
	void CheckErrorDistMultiVector( NumMat<Scalar>& xLocal, NumMat<Scalar>& xTrueLocal );

	/// @brief LUstructToPMatrix converts the data in LUstruct to PMatrix.
	///
	/// @paramp[out] PMloc
	void LUstructToPMatrix( PMatrix& PMloc );

	/// @brief SymbolicToSuperNode converts the symbolic information to
	/// SuperNode structure in SelInv.
	///
	/// @param[out] super
	void SymbolicToSuperNode( SuperNode& super );
};





} // namespace PEXSI

#endif // _SUPERLU_DIST_INTERF_HPP_

