/// @file pselinv.hpp
/// @brief Header file for PSelInv
/// @author Lin Lin
/// @version 0.1
/// @date 2012-11-14
#ifndef _PSELINV_HPP_
#define _PSELINV_HPP_

// *********************************************************************
//  Common utilities
// *********************************************************************

#include  "environment_impl.hpp"
#include	"numvec_impl.hpp"
#include	"nummat_impl.hpp" 
#include  "sparse_matrix.hpp"
#include  "mpi_interf.hpp"
#include	"utility.hpp"
#include	"blas.hpp"
#include	"lapack.hpp"

namespace PEXSI{

//const char UPPER = 'U';
//const char LOWER = 'L';
//const char FROMRIGHT = 'R';
//const char FROMLEFT  = 'L';
//const char TRAN = 'T';
//const char NOTRAN = 'N';
//const char ALL   = 'A';
//const char UDIAG = 'U';
//const char NOUDIAG = 'N';


/**********************************************************************
 * Basic PSelInv data structure
 **********************************************************************/

/// @struct Grid 
///
/// @brief Grid is the PSelInv way of defining the grid.  
///
/// Grid should be consistent with the grid used by SuperLU.
///
/// NOTE: It is your responsibility to make sure that the SuperLUGrid
/// and Grid used for SelInv are the same.
struct Grid{
	// Data
  MPI_Comm    comm;
  MPI_Comm    rowComm;
	MPI_Comm    colComm;
	Int         mpirank;
  Int         mpisize; 
	Int         numProcRow;
	Int         numProcCol;

	// Member function
	Grid( MPI_Comm Bcomm, int nprow, int npcol );
	~Grid();
};

/// @struct SuperNode
///
/// @brief SuperNode describes mapping between supernode and column, the
/// permutation information, and potentially the elimination tree (not
/// implemented here).
/// 
/// superIdx[i] is the supernode index to which column i belongs. 
/// This is the same as supno[i] in SuperLU.
///
/// superPtr[s] is the leading column of the s-th supernode (as in
/// colptr).  This is the same as xsup[s] in SuperLU.
///
///	e.g.   superIdx  0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
///	       superPtr  0 1 2 4 7 12
///
/// This is allocated during symbolic factorization SYMBFACT.
///
/// perm is the permutation vector.  Symmetric permutation is assumed.
/// perm is the same as ScalePermstruct -> perm_c.
///
struct SuperNode{
  IntNumVec   perm;              
	IntNumVec   superIdx;
	IntNumVec   superPtr;
};

/// @struct LBlock
///
/// @brief LBlock stores a lower triangular block, including the
/// diagonal block in PSelInv.
struct LBlock{
	// Variables
	Int               blockIdx;
	Int               numRow;
	Int               numCol;
	IntNumVec         rows;
	NumMat<Scalar>    nzval;

	// Member functions;
	LBlock() {blockIdx = -1; numRow = 0; numCol =0;}
	~LBlock() {}
	LBlock& operator = (const LBlock& LB) {
		blockIdx    = LB.blockIdx;
		numRow      = LB.numRow;
		numCol      = LB.numCol;
		rows        = LB.rows;
		nzval       = LB.nzval;
		return *this;
	}
};

/// @struct UBlock
///
/// @brief UBlock stores a upper triangular block, excluding the
/// diagonal block in PSelInv.
struct UBlock{
	// Variables
	Int               blockIdx;
	Int               numRow;
	Int               numCol;
	IntNumVec         cols;
	NumMat<Scalar>    nzval;

	// Member functions;
	UBlock() {blockIdx = -1; numRow = 0; numCol =0;}
	~UBlock() {}
	UBlock& operator = (const UBlock& UB) {
		blockIdx    = UB.blockIdx;
		numRow      = UB.numRow;
		numCol      = UB.numCol;
		cols        = UB.cols;
		nzval       = UB.nzval;
		return *this;
	}
};

// *********************************************************************
// SuperLU style utility functions
// 
// The define MACROS are removed so that the code is more portable.
// *********************************************************************

inline Int MYPROC( const Grid* g )
{ return g->mpirank; }

inline Int MYROW( const Grid* g )
{ return g->mpirank / g->numProcCol; }

inline Int MYCOL( const Grid* g )
{ return g->mpirank % g->numProcCol; }

inline Int PROW( Int bnum, const Grid* g ) 
{ return bnum % g->numProcRow; }

inline Int PCOL( Int bnum, const Grid* g ) 
{ return bnum % g->numProcCol; }

inline Int PNUM( Int i, Int j, const Grid* g )
{ return i * g->numProcCol + j; }

inline Int LBi( Int bnum, const Grid* g )
{ return bnum / g->numProcRow; }

inline Int LBj( Int bnum, const Grid* g)
{ return bnum / g->numProcCol; }

inline Int GBi( Int iLocal, const Grid* g )
{ return iLocal * g->numProcRow + MYROW( g ); }

inline Int GBj( Int jLocal, const Grid* g )
{ return jLocal * g->numProcCol + MYCOL( g ); }

inline Int CEILING( Int a, Int b )
{ return (a%b) ? ( a/b + 1 ) : ( a/b ); }

inline Int BlockIdx( Int i, const SuperNode *s )
{ return s->superIdx[i]; }

inline Int FirstBlockCol( Int bnum, const SuperNode *s )
{ return s->superPtr[bnum]; }	

inline Int SuperSize( Int bnum, const SuperNode *s )
{ return s->superPtr[bnum+1] - s->superPtr[bnum]; } 

inline Int NumSuper( const SuperNode *s )
{ return s->superPtr.m() - 1; }

inline Int NumCol( const SuperNode *s )
{ return s->superIdx.m(); }


// *********************************************************************
// Serialize / Deserialize
// *********************************************************************

// L part

namespace LBlockMask{
enum {
	BLOCKIDX,
	NUMROW,
	NUMCOL,
	ROWS,
	NZVAL,
	TOTAL_NUMBER
};
}

Int inline serialize(LBlock& val, std::ostream& os, const std::vector<Int>& mask){
  Int i = 0;
  if(mask[i]==1) serialize(val.blockIdx, os, mask); i++;
  if(mask[i]==1) serialize(val.numRow,  os, mask); i++;
  if(mask[i]==1) serialize(val.numCol,  os, mask); i++;
  if(mask[i]==1) serialize(val.rows, os, mask);   i++;
  if(mask[i]==1) serialize(val.nzval, os, mask);  i++;
  return 0;
}

Int inline deserialize(LBlock& val, std::istream& is, const std::vector<Int>& mask){
  Int i = 0;
  if(mask[i]==1) deserialize(val.blockIdx, is, mask); i++;
  if(mask[i]==1) deserialize(val.numRow,  is, mask); i++;
  if(mask[i]==1) deserialize(val.numCol,  is, mask); i++;
  if(mask[i]==1) deserialize(val.rows,   is, mask); i++;
  if(mask[i]==1) deserialize(val.nzval,  is, mask); i++; 
  return 0;
}

// U part
namespace UBlockMask{
enum {
	BLOCKIDX,
	NUMROW,
	NUMCOL,
	COLS,
	NZVAL,
	TOTAL_NUMBER
};
}

Int inline serialize(UBlock& val, std::ostream& os, const std::vector<Int>& mask){
  Int i = 0;
  if(mask[i]==1) serialize(val.blockIdx, os, mask); i++;
  if(mask[i]==1) serialize(val.numRow,  os, mask); i++;
  if(mask[i]==1) serialize(val.numCol,  os, mask); i++;
  if(mask[i]==1) serialize(val.cols, os, mask);   i++;
  if(mask[i]==1) serialize(val.nzval, os, mask);  i++;
  return 0;
}

Int inline deserialize(UBlock& val, std::istream& is, const std::vector<Int>& mask){
  Int i = 0;
  if(mask[i]==1) deserialize(val.blockIdx, is, mask); i++;
  if(mask[i]==1) deserialize(val.numRow,  is, mask); i++;
  if(mask[i]==1) deserialize(val.numCol,  is, mask); i++;
  if(mask[i]==1) deserialize(val.cols,   is, mask); i++;
  if(mask[i]==1) deserialize(val.nzval,  is, mask); i++; 
  return 0;
}


/**********************************************************************
 * Main data structure in PSelInv: PMatrix
 **********************************************************************/

/// @class PMatrix
///
/// @brief PMatrix contains the main data structure and the
/// computational routine for the parallel selected inversion.  
///
/// NOTE: In the current version of PMatrix, square grid is assumed.
/// This assumption is only used when sending the information to
/// cross-diagonals, i.e. from L_{ik} to U_{ki}.  This assumption can be
/// relaxed later.
/// 

class PMatrix{

private:
	// *********************************************************************
	// Variables
	// *********************************************************************
	// Data variables
	const Grid*           grid_;
	const SuperNode*      super_;
	std::vector<std::vector<LBlock> > L_;
	std::vector<std::vector<UBlock> > U_;

	// Communication variables
	NumMat<bool>                       isSendToDown_;
	NumMat<bool>                       isSendToRight_;
	NumVec<bool>                       isSendToCrossDiagonal_;

	NumVec<bool>                       isRecvFromUp_;
	NumVec<bool>                       isRecvFromLeft_;
	NumVec<bool>                       isRecvFromCrossDiagonal_;

private:
	// *********************************************************************
	// Private member functions 
	// *********************************************************************
  
	/// @brief PreSelInv computes the inverse of the diagonal blocks, 
	/// scale the off-diagonal L blocks by L_{kk}^{-1}, and fill the
	/// U_{ki} blocks by L_{ik}.
	///
	/// PreSelInv assumes that ConstructCommunicationPattern has been
	/// executed.
	void PreSelInv( );

public:
	// *********************************************************************
	// Public member functions 
	// *********************************************************************

	PMatrix( const PEXSI::Grid* g, const PEXSI::SuperNode* s );

	~PMatrix();

	Int NumCol() const { return super_ -> superIdx.m(); }

	Int NumSuper() const { return super_ ->superPtr.m() - 1; }

	Int NumLocalBlockCol() const { return CEILING( this->NumSuper(), grid_->numProcCol ); }

	Int NumLocalBlockRow() const { return CEILING( this->NumSuper(), grid_->numProcRow); }

	/// @brief NumBlockL returns the number of nonzero L blocks for the
	/// local block column jLocal.
	Int NumBlockL( Int jLocal ) const { return L_[jLocal].size(); }

	/// @brief NumBlockU returns the number of nonzero U blocks for the
	/// local block row iLocal.
  Int NumBlockU( Int iLocal ) const { return U_[iLocal].size(); }

	const PEXSI::Grid* Grid() const { return grid_; }
  
	const PEXSI::SuperNode* SuperNode() const { return super_; }	

	/// @brief L returns the vector of nonzero L blocks for the local
	/// block column jLocal.
  std::vector<LBlock>& L( Int jLocal ) { return L_[jLocal]; } 	

	/// @brief U returns the vector of nonzero U blocks for the local
	/// block row iLocal.
	std::vector<UBlock>& U( Int iLocal ) { return U_[iLocal]; }

	/// @brief ConstructCommunicationPattern constructs the communication
	/// pattern to be used later in the selected inversion stage.
	void ConstructCommunicationPattern( );

	/// @brief SelInv is the main function for the selected inversion.
	///
	/// Algorithm description for PSelInv
	/// ---------------------------------
	///
	/// PSelInv is a right-looking based parallel selected inversion
	/// subroutine for sparse symmetric matrices.  Static pivoting is
	/// used in this version.
	///
	/// Although the matrix is symmetric, the key idea of the current
	/// implementation of PSelInv is that the upper-triangular matrix is
	/// saved (in the form of UBlock).  Such redundant information is
	/// effective for reducing the complexity of the communication
	/// pattern.  
	///
	/// At each supernode ksup, three subroutines are executed:
	///
	/// - UpdateL.
	/// 
	/// The blocks in the processor row of ksup first sends the nonzero
	/// blocks in U(ksup, jsup) to the Schur complements Ainv(isup, jsup).
	/// At the same time the blocks in the processor column of ksup sends
	/// the nonzero blocks (only row indices) to the Schur complement
	/// Ainv(isup, jsup).  Then 
	///
	/// \sum_{jsup} A^{-1}(isup, jsup) U^{T}(ksup, jsup)
	///
	/// is performed.  In this procedure, only processors with
	/// isRecvFromUp[ksup] == true && isRecvFromLeft[ksup] == true
	/// participate in the computation.
	///
	///
	/// The result is reduced to the processor column ksup
	/// within the same processor row. 
	///
	/// - UpdateD.
	///
	/// The diagonal block (ksup, ksup) is simply updated by a reduce
	/// procedure within the column processor group of ksup.
	///
	/// - UpdateU.
	///
	/// The Ainv(ksup, isup) blocks can be simply updated via
	/// the update from the cross diagonal processors.
	///
	/// The cross diagonal processor is only well defined for square
	/// grids.  For a P x P square grid, (ip, jp) is the cross diagonal
	/// processor of (jp, ip) if ip != jp.  The current version of SelInv
	/// only works for square processor grids.
	///
	///
	/// Communication pattern
	/// ---------------------
	///
	/// The communication is controlled by 3 sending varaibles and 3
	/// receiving variables. The first dimension of all the sending and
	/// receiving variables are numSuper.  There is redundant information
	/// saved, but this increases the readability of the output by only
	/// increasing a small amount of memory cost for indexing.  This set of
	/// sending / receiving mechanism avoids the double indexing of the
	/// supernodes and can scale to matrices of large size.
	///
	/// - isSendToDown:  
	///
	///   Dimension: numSuper x numProcRow
	///
	///   Role     : At supernode ksup, if isSendToDown(ksup, ip) == true, send
	///   all local blocks {(ksup, jsup) | jsup > ksup} to the processor row ip.
	///
	/// - isRecvFromUp:
	///
	///   Dimension: numSuper
	///
	///   Role     : At supernode ksup, if isRecvFromUp(ksup) == true,
	///   receive blocks from the processor owning the block row of ksup
	///   within the same column processor group.
	///
	/// - isSendToRight:
	///
	///   Dimension: numSuper x numProcCol
	///
	///   Role     : At supernode ksup, if isSendToRight(jp, ksup) == true, send
	///   all local blocks {(isup, ksup) | isup > ksup} to the processor
	///   column jp.
	///
	/// - isRecvFromLeft:
	///   
	///   Dimension: numSuper
	///
	///   Role     : At supernode ksup, if isRecvFromLeft(ksup) == true,
	///   receive blocks from the processor owning the block column of
	///   ksup within the same row processor group.
	///
	/// - isSendToCrossDiagonal:
	///
	///   Dimension: numSuper
	///
	///   Role     : At supernode ksup, if isSendToCrossDiagonal(ksup) ==
	///   true, send all local blocks {(isup, ksup) | isup > ksup} to the
	///   cross-diagonal processor.  NOTE: This requires a square processor
	///   grid.
	///
	/// - isRecvCrossDiagonal:
	///
	///   Dimension: numSuper
	///
	///   Role     : At supernode ksup, if isRecvFromCrossDiagonal(ksup) ==
	///   true, receive from the cross-diagonal processor.  NOTE: This
	///   requires a square processor grid.
	///   
	///
	/// Future work
	/// -----------
	///
	/// - Pipelining:
	///
	///   Approximately 50% of speedup.
	///   
	///   Reference: S. Li and J. Demmel, SuperLU_DIST : A Scalable
	///   Distributed-Memory Sparse Direct Solver for Unsymmetric Linear
	///   Systems, ACM TOMS, 2003
	///
	/// - Look-ahead and static scheduling:
	///
	///   Approximately 2.5~3.5 times speedup.
	///
	///   Reference: I. Yamazaki and S. Li, New Scheduling Strategies and Hybrid Programming for
	///   a Parallel Right-looking Sparse LU Factorization Algorithm on
	///   Multicore Cluster Systems, IPDPS 2012
	///
	///
	void SelInv( );

	/// @brief PMatrixToDistSparseMatrix converts the PMatrix into a
	/// distributed compressed sparse column matrix format, according
	/// block row partition.
	///
	/// @param[in] A Input sparse matrix to enforce the sparsity
	/// pattern.
	/// @param[out] B Output sparse matrix with the same sparsity pattern
	/// as A.
	void PMatrixToDistSparseMatrix( 
			const DistSparseMatrix<Scalar>& A,
			DistSparseMatrix<Scalar>& B );

	/// @brief Diagonal extracts the diagonal elements of the PMatrix.
	///
	/// 1) diagNaturalOrder is permuted back to the natural order
	///
	/// 2) diagNaturalOrder is shared by all processors in grid_->comm through a
	/// Allreduce procedure.
	void Diagonal( NumVec<Scalar>& diagNaturalOrder );

};



} // namespace PEXSI

#endif // _PSELINV_HPP_
