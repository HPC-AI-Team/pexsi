#ifndef _SUPERLU_DIST_INTERF_HPP_
#define _SUPERLU_DIST_INTERF_HPP_

// Interface with PSelInv
#include "pselinv.hpp"

// Interface with sparse matrix (CSC format)
#include  "sparse_matrix.hpp"

// Interface with LAPACK
#include  "lapack.hpp"

namespace PEXSI{

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

/// @class SuperLUMatrix
/// @brief An thin interface to keep the main code insulated from
/// the source code of SuperLU. 
class SuperLUMatrix{
private:
	/// @struct SuperLUData
	///
	/// @brief  Data of a matrix in the SuperLU format.
	///
	/// SuperLUData is only used to define ptrData which is a private member
	/// SuperLUMatrix.
	struct        SuperLUData;
	/// @brief Actual pointer to save the data related to SuperLU.
	SuperLUData*  ptrData;
public:

	SuperLUMatrix( const SuperLUGrid& g );

	~SuperLUMatrix();

	int m() const;

	int n() const;

	/// @brief DistSparseMatrixToSuperMatrixNRloc converts a distributed
	/// sparse matrix in compressed sparse column format into the SuperLU
	/// compressed row format. 
	///
	/// TODO real arithmetic
	/// FIXME: IntNumVec convention.  Assumes a symmetric matrix
	void DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Scalar>& sparseA );

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





//void LUstructToPMatrix( const SuperLUData& luData, Grid &grid, PMatrix& PMloc );
//
//void ConstructCommunicationPattern( CommunicationPattern& commPattern, 
//		const Grid& grid, const PMatrix& PMloc);
//



} // namespace PEXSI

#endif // _SUPERLU_DIST_INTERF_HPP_

