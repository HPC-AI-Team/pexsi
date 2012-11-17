#ifndef _PSELINV_HPP_
#define _PSELINV_HPP_

// *********************************************************************
//  Common utilities
// *********************************************************************

#include  "environment_impl.hpp"
#include	"numvec_impl.hpp"
#include	"nummat_impl.hpp" 

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
 * Output LU structure
 **********************************************************************/

/// @struct Grid 
///
/// @brief Grid is the PSelInv way of defining the grid.  
///
/// Grid should be consistent with the grid used by SuperLU.
struct Grid{
	// Data
  MPI_Comm    comm;
  MPI_Comm    rowComm;
	MPI_Comm    colComm;
  Int         mpisize; 
	Int         numProcRow;
	Int         numProcCol;

	// Member function
	Grid( MPI_Comm Bcomm, int nprow, int npcol );
	~Grid();
};

/// @struct SuperNode
///
/// @brief SuperNode describes mapping between supernode and column.
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
struct SuperNode{
	IntNumVec   superPtr;
	IntNumVec   superIdx;
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

class PMatrix{
private:
	Grid*           grid;
//private:
//
//	Int _ndof;
//	Int _nsupers;
//	Int _nbc;
//	Int _nbr;
//	IntNumVec _xsup;
//	IntNumVec _supno;
//	IntNumVec _nblkl;
//	IntNumVec _nblku;
//	std::vector<std::vector<LBlock> > _L;
//	std::vector<std::vector<UBlock> > _U;
//
//public:
//	Int _ndof;
//	Int _nsupers;
//	Int _nbc;
//	Int _nbr;
//	IntNumVec _xsup;
//	IntNumVec _supno;
//	IntNumVec _nblkl;
//	IntNumVec _nblku;
//	std::vector<std::vector<LBlock> > _L;
//	std::vector<std::vector<UBlock> > _U;
public:
	PMatrix(){;}
	~PMatrix(){;}
//	Int ndof(){return _ndof;}
//	Int nsupers(){return _nsupers;}
//	Int nbc(){return _nbc;}
//	Int nbr(){return _nbr;}
//	Int xsup(Int i){return _xsup[i];}
//	Int supsize(Int i){return _xsup[i+1]-_xsup[i];}
//	Int superno(Int i){return _supno[i];}
//	Int nblkl(Int i){return _nblkl[i];}
//	Int nblku(Int i){return _nblku[i];}
//	std::vector<std::vector<LBlock> >& L(){return _L;}
//	std::vector<LBlock>& L(Int i){return _L[i];}
//	std::vector<std::vector<UBlock> >& U(){return _U;}
//	std::vector<UBlock>& U(Int i){return _U[i];}
//
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

//class BlockPtn
//{
//public:
//	IntNumVec _ownerinfo;
//public:
//	BlockPtn() {;}
//	~BlockPtn() {;}
//	IntNumVec& ownerinfo() { return _ownerinfo; }
//	Int owner(Int key) {
//		return _ownerinfo(key);
//	}
//};

// *********************************************************************
// Main functions for PLUSelInv
// *********************************************************************

//#define KITEMS  5
//#define ISOURCE 0
//#define ITARGET 1
//#define ISUPER  2
//#define IROWIND 3
//#define NROWS   4

//Int SuperLU2SelInv(Int n, LUstruct_t *LUstruct, gridinfo_t *grid,
//		PMatrix& PMloc);
//
//// Convert a CSC Matrix to PMatrix structure in order to perform trace
//// operations. 
////void CSCMatrixToPMatrix();
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


// *********************************************************************
// SuperLU style utility functions
// 
// The define MACROS are removed so that the code is more portable.
// *********************************************************************

inline Int MYROW( const Int mpirank, const Grid* g )
{ return mpirank / g->numProcCol; }

inline Int MYCOL( const Int mpirank, const Grid* g )
{ return mpirank % g->numProcCol; }

inline Int PROW( const Int bnum, const Grid* g ) 
{ return bnum % g->numProcRow; }

inline Int PCOL( const Int bnum, const Grid* g ) 
{ return bnum % g->numProcCol; }

inline Int PNUM( const Int i, const Int j, const Grid* g )
{ return i * g->numProcCol + j; }

inline Int LBi( const Int bnum, const Grid* g)
{ return bnum / g->numProcRow; }

inline Int LBj( const Int bnum, const Grid* g)
{ return bnum / g->numProcCol; }

inline Int CEILING( const Int a, const Int b )
{ return (a%b) ? ( a/b + 1 ) : ( a/b ); }

inline Int BlockIdx( const Int i, const SuperNode *s )
{ return s->superIdx[i]; }

inline Int FirstBlockCol( const Int bnum, const SuperNode *s )
{ return s->superPtr[bnum]; }	

inline Int SuperSize( const Int bnum, const SuperNode *s )
{ return s->superPtr[bnum+1] - s->superPtr[bnum]; } 

} // namespace PEXSI

#endif // _PSELINV_HPP_
