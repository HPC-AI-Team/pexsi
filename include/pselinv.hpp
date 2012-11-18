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

inline Int LBi( Int bnum, const Grid* g)
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
	std::vector<std::vector<bool> >    isSendToDown;
	std::vector<std::vector<bool> >    isSendToRight;
	std::vector<std::vector<bool> >    isSendToCrossDiagonal;

	std::vector<bool>                  isRecvFromUp;
	std::vector<bool>                  isRecvFromLeft;
	std::vector<bool>                  isRecvFromCrossDiagonal;

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

	void UpdateL( Int ksup );

	void UpdateD( Int ksup );

	void UpdateU( Int ksup );

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

//	Int DumpL(string filename, gridinfo_t* grid);
//	Int DumpU(string filename, gridinfo_t* grid);
//
//	Int DumpLBlock(Int isup, Int jsup, string filename, gridinfo_t* grid);
//	Int DumpUBlock(Int isup, Int jsup, string filename, gridinfo_t* grid);
//
//	//    Int CondDiagBlock(gridinfo_t* grid);
//	Int DumpDiagVec(NumVec<Scalar>& globalDiagVec, 
//			string filename, gridinfo_t *grid);
//
};

//
//// Convert PMatrix structure to a CSC matrix in order to perform trace
//// operations. 
//void PMatrixToCSCMatrix(Int n, gridinfo_t *grid, PMatrix& PMloc);
//
//
//
//
//
//Int DiagTri(PMatrix& PMloc, Int ksup, gridinfo_t *grid);
//Int ScaleLinv(PMatrix& PMloc, Int ksup, gridinfo_t *grid);
//Int ScaleUinv(PMatrix& PMloc, Int ksup, gridinfo_t *grid);
//Int DiagInnerProd(PMatrix& PMloc, std::vector<std::vector<Int> >& localEtree, Int ksup, gridinfo_t *grid);
//Int SinvPL(PMatrix& PMloc, std::vector<std::vector<Int> >& localEtree, Int ksup, gridinfo_t* grid);
//Int SinvPU(PMatrix& PMloc, std::vector<std::vector<Int> >& localEtree, Int ksup, gridinfo_t* grid);
//
//Int PLUSelInv(gridinfo_t *grid, PMatrix& PMloc, std::vector<std::vector<Int> >& localEtree);
//Int ConstructLocalEtree(Int n, gridinfo_t *grid, PMatrix& PMloc, 
//		std::vector<std::vector<Int> >& localEtree);
//
//Int BcastLBlock(LBlock& LB, Int mykey, Int srckey, MPI_Comm comm);
//Int BcastUBlock(UBlock& UB, Int mykey, Int srckey, MPI_Comm comm);
//
//Int DumpLocalEtree(std::vector<std::vector<Int> >& localEtree, gridinfo_t* grid);

//void blockinvert(Lstruct *PMlocL, LUstruct_t *LUstruct, gridinfo_t *grid,
//                 Int jsup, Int n);
//void ScalePMbyL(Lstruct *PMlocL, gridinfo_t *grid, Int ksup, Int nsupers);
//void SetPL(Lstruct *PMlocL, gridinfo_t *grid, PLstruct *PL, 
//           Int ksup, Int nsupers);



} // namespace PEXSI

#endif // _PSELINV_HPP_
