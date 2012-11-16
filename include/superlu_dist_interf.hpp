#ifndef _SUPERLU_DIST_INTERF_HPP_
#define _SUPERLU_DIST_INTERF_HPP_

// Interface with PSelInv
//#include "pselinv.hpp"

// Interface with sparse matrix (CSC format)
#include "sparse_matrix.hpp"


namespace PEXSI{

/// @brief A thin interface for the gridinfo_t strucutre in SuperLU.
class SuperLUGrid{
	friend class SuperLUMatrix;
private:
	struct           GridData;
	GridData         *ptrData;
public:
	SuperLUGrid( MPI_Comm comm, int nprow, int npcol );

	~SuperLUGrid();
};

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
	SuperLUData   *ptrData;
public:

	SuperLUMatrix();

	~SuperLUMatrix();

	void Initialize( MPI_Comm comm, int nprow, int npcol  );

	void Finalize();

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

//	void DistMultiplyGlobalVector( Scalar*);

};


//void LUstructToPMatrix( const SuperLUData& luData, Grid &grid, PMatrix& PMloc );
//
//void ConstructCommunicationPattern( CommunicationPattern& commPattern, 
//		const Grid& grid, const PMatrix& PMloc);
//
//void SuperLUSymbolicFactorization( SuperLUData& luData );
//
//void SuperLURedistribute(
//
//
//void DistSparseMatrixToSuperMatrixNRloc( SuperLUData& luData, DistSparseMatrix<Complex>& A );



} // namespace PEXSI

#endif // _SUPERLU_DIST_INTERF_HPP_

