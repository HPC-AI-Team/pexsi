/// @file superlu_dist_interf.cpp
/// @brief Interface to SuperLU_Dist (v3.2).
/// @author Lin Lin
/// @date 2012-11-15
#include "superlu_dist_interf.hpp"
#ifdef _USE_COMPLEX_
#include "superlu_zdefs.h"
#include "Cnames.h"
extern "C"{
void
pzsymbfact(superlu_options_t *options, SuperMatrix *A, 
		ScalePermstruct_t *ScalePermstruct, gridinfo_t *grid,
		LUstruct_t *LUstruct, SuperLUStat_t *stat, int *info);
}
#else
// TODO pdsymbfact
#include "superlu_ddefs.h"
#include "Cnames.h"
extern "C"{ void
pdsymbfact(superlu_options_t *options, SuperMatrix *A, 
		ScalePermstruct_t *ScalePermstruct, gridinfo_t *grid,
		LUstruct_t *LUstruct, SuperLUStat_t *stat, int *info);
}
#endif

// SuperLUGrid class
namespace PEXSI{

struct SuperLUGrid::GridData{
	gridinfo_t          grid;
};


SuperLUGrid::SuperLUGrid	( MPI_Comm comm, Int nprow, Int npcol )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUGrid::SuperLUGrid");
#endif
	ptrData = new GridData;
	if( ptrData == NULL )
		throw std::runtime_error( "SuperLUGrid cannot be allocated." );

	superlu_gridinit(comm, nprow, npcol, &ptrData->grid);


#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUGrid::SuperLUGrid  ----- 


SuperLUGrid::~SuperLUGrid	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUGrid::~SuperLUGrid");
#endif
	// NOTE: Do not call 
	//
	// superlu_gridexit(&ptrData->grid);
	//
	// Since this subroutine frees grid->comm.  This is not what you want
	// to do when grid->comm is to be used by other subroutines.

	// The following lines are taken from
	// superlu_gridexit:SuperLU_DIST/SRC/superlu_grid.c 
	MPI_Comm_free( &ptrData->grid.rscp.comm );
	MPI_Comm_free( &ptrData->grid.cscp.comm );
	// SuperLU defines this new  type in superlu_grid.c
#ifdef _USE_COMPLEX_
	if ( SuperLU_MPI_DOUBLE_COMPLEX != MPI_DATATYPE_NULL ) {
		MPI_Type_free( &SuperLU_MPI_DOUBLE_COMPLEX );
	}
#endif
	
	delete ptrData;

#ifndef _RELEASE_
	PopCallStack();
#endif
	return ;
} 		// -----  end of method SuperLUGrid::~SuperLUGrid  ----- 

} // namespace PEXSI

// SuperLUMatrix class
namespace PEXSI{

struct SuperLUMatrix::SuperLUData{
	/// @brief SuperLU matrix. 
	SuperMatrix         A;                        

	/// @brief SuperLU options. 
	///
	/// Note
	/// ----
	///
	/// It is important to have 
	///
	/// options.RowPerm           = NOROWPERM;
	/// 
	/// to make sure that symmetric permutation is used.
	///
	superlu_options_t   options;                  

	/// @brief Saves the permutation vectors.  Only perm_c (permutation of
	/// column as well as rows due to the symmetric permutation) will be used.
	ScalePermstruct_t   ScalePermstruct;          

	/// @brief SuperLU grid structure.
	gridinfo_t*         grid;

	/// @brief Saves the supernodal partition as well as the numerical
	/// values and structures of the L and U structure.
	LUstruct_t          LUstruct;

	/// @brief Used for solve for multivectors.
	SOLVEstruct_t       SOLVEstruct;

	/// @brief SuperLU statistics
	SuperLUStat_t       stat;

	/// @brief SuperLU information
	Int                 info;

	// The following are for consistency checks
	bool                isSuperMatrixAllocated;
	bool                isSuperMatrixFactorized;
	bool                isScalePermstructAllocated;
	bool                isLUstructAllocated;
};

SuperLUMatrix::SuperLUMatrix	( const SuperLUGrid& g, const SuperLUOptions& opt )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::SuperLUMatrix");
#endif
	ptrData = new SuperLUData;
	if( ptrData == NULL )
		throw std::runtime_error( "SuperLUMatrix cannot be allocated." );
	
	ptrData->isSuperMatrixAllocated     = false;
	ptrData->isScalePermstructAllocated = false;
	ptrData->isLUstructAllocated        = false;

	// Options
	superlu_options_t& options = ptrData->options;
	
	set_default_options_dist(&options);

	// The default value of ColPerm uses the default value from SuperLUOptions
	options.Fact              = DOFACT;
	options.RowPerm           = NOROWPERM; // IMPORTANT for symmetric matrices
	options.IterRefine        = NOREFINE;
	options.ParSymbFact       = NO;
	options.Equil             = NO; 
	options.ReplaceTinyPivot  = NO;
	options.PrintStat         = YES;
	options.SolveInitialized  = NO;

	if ( opt.ColPerm == "NATURAL" ){
		options.ColPerm = NATURAL;
	} 
	else if( opt.ColPerm == "MMD_AT_PLUS_A" ){
		options.ColPerm = MMD_AT_PLUS_A;
	}
	else if( opt.ColPerm == "METIS_AT_PLUS_A" ){
		options.ColPerm = METIS_AT_PLUS_A;
	}
	else if( opt.ColPerm == "PARMETIS" ){
		options.ColPerm           = PARMETIS;
		options.ParSymbFact       = YES;
//		options.ParSymbFact       = NO;
	}
	else{
		std::ostringstream msg;
		msg << opt.ColPerm << " is not a supported ColPerm type. Try (case sensitive) " << std::endl
			  << "NATURAL | MMD_AT_PLUS_A | METIS_AT_PLUS_A | PARMETIS" << std::endl;
		throw std::runtime_error( msg.str().c_str() );
	}

	// Setup grids
  ptrData->grid = &(g.ptrData->grid);

#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix::SuperLUMatrix  ----- 


SuperLUMatrix::~SuperLUMatrix	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::~SuperLUMatrix");
#endif
	SuperMatrix& A = ptrData->A;
	if( ptrData->isLUstructAllocated ){
		Destroy_LU(A.ncol, ptrData->grid, &ptrData->LUstruct);
		LUstructFree(&ptrData->LUstruct); 
	}
	if( ptrData->isScalePermstructAllocated ){
		ScalePermstructFree(&ptrData->ScalePermstruct);
	}
	if( ptrData->options.SolveInitialized ){
		// TODO real arithmetic
		zSolveFinalize(&ptrData->options, &ptrData->SOLVEstruct);
	}
	if( ptrData->isSuperMatrixAllocated ){
		this->DestroyAOnly();
	}
	
	delete ptrData;

#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SuperLUMatrix::~SuperLUMatrix  ----- 


Int SuperLUMatrix::m (  ) const	
{
	return ptrData->A.nrow;
} 		// -----  end of method SuperLUMatrix::m  ----- 


Int SuperLUMatrix::n (  ) const	
{
	return ptrData->A.ncol;
} 		// -----  end of method SuperLUMatrix::n  ----- 

void SuperLUMatrix::DistSparseMatrixToSuperMatrixNRloc( DistSparseMatrix<Scalar>& sparseA ){
#ifndef _RELEASE_
	PushCallStack( "SuperLUMatrix::DistSparseMatrixToSuperMatrixNRloc" );
#endif
	if( ptrData->isSuperMatrixAllocated == true ){
		throw std::logic_error( "SuperMatrix is already allocated." );
	}
	gridinfo_t* grid = ptrData->grid;

	Int mpirank = grid->iam;
	Int mpisize = grid->nprow * grid->npcol;

	Int numRowLocal = sparseA.colptrLocal.m() - 1;
	Int numRowLocalFirst = sparseA.size / mpisize;
	Int firstRow = mpirank * numRowLocalFirst;

  int_t *colindLocal, *rowptrLocal;
	doublecomplex *nzvalLocal;
	rowptrLocal = (int_t*)intMalloc_dist(numRowLocal+1);
	colindLocal = (int_t*)intMalloc_dist(sparseA.nnzLocal); 
	nzvalLocal  = (doublecomplex*)doublecomplexMalloc_dist(sparseA.nnzLocal);
  
	std::copy( sparseA.colptrLocal.Data(), sparseA.colptrLocal.Data() + sparseA.colptrLocal.m(),
			rowptrLocal );
	std::copy( sparseA.rowindLocal.Data(), sparseA.rowindLocal.Data() + sparseA.rowindLocal.m(),
			colindLocal );
	std::copy( sparseA.nzvalLocal.Data(), sparseA.nzvalLocal.Data() + sparseA.nzvalLocal.m(),
			(Complex*)nzvalLocal );



	// Important to adjust from FORTRAN convention (1 based) to C convention (0 based) indices
	for(Int i = 0; i < sparseA.rowindLocal.m(); i++){
		colindLocal[i]--;
	}

	for(Int i = 0; i < sparseA.colptrLocal.m(); i++){
		rowptrLocal[i]--;
	}

	// Construct the distributed matrix according to the SuperLU_DIST format
	zCreate_CompRowLoc_Matrix_dist(&ptrData->A, sparseA.size, sparseA.size, sparseA.nnzLocal, 
			numRowLocal, firstRow,
			nzvalLocal, colindLocal, rowptrLocal,
			SLU_NR_loc, SLU_Z, SLU_GE);

	ptrData->isSuperMatrixAllocated = true;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return;

} 		// -----  end of method SuperLUMatrix::DistSparseMatrixToSuperMatrixNRloc ----- 


void
SuperLUMatrix::DestroyAOnly	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::DestroyAOnly");
#endif
	if( ptrData->isSuperMatrixAllocated == false ){
		throw std::logic_error( "SuperMatrix has not been allocated." );
	}
	switch ( ptrData->A.Stype ){
		case SLU_NC:
			Destroy_CompCol_Matrix_dist(&ptrData->A);
			break;
		case SLU_NR_loc:
			Destroy_CompRowLoc_Matrix_dist(&ptrData->A);
			break;
		default:
			std::ostringstream msg;
			msg << "Type " << SLU_NR_loc << " is to be destroyed" << std::endl
				<< "This is an unsupported SuperMatrix format to be destroyed." << std::endl;
			throw std::runtime_error( msg.str().c_str() );
	}
	ptrData->isSuperMatrixAllocated = false;
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix::DestroyAOnly  ----- 

void
SuperLUMatrix::SymbolicFactorize	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::SymbolicFactorize");
#endif
	if( ptrData->isScalePermstructAllocated ){
		throw std::logic_error( "ScalePermstruct is already allocated." );
	}
	if( ptrData->isLUstructAllocated){
		throw std::logic_error( "LUstruct is already allocated." );
	}
	if( ptrData->options.RowPerm != NOROWPERM ){
		throw std::logic_error( "For PEXSI there must be no row permutation." );
	}
	
	SuperMatrix&  A = ptrData->A;

	ScalePermstructInit(A.nrow, A.ncol, &ptrData->ScalePermstruct);
	LUstructInit(A.nrow, A.ncol, &ptrData->LUstruct);

	PStatInit(&ptrData->stat);
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "Before symbfact subroutine." << std::endl;
#endif
	pzsymbfact(&ptrData->options, &A, &ptrData->ScalePermstruct, ptrData->grid, 
			&ptrData->LUstruct, &ptrData->stat, &ptrData->info);
	PStatFree(&ptrData->stat);

	ptrData->isScalePermstructAllocated = true;
  ptrData->isLUstructAllocated        = true; 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix::SymbolicFactorize  ----- 


void
SuperLUMatrix::Distribute	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::Distribute");
#endif
	if( ptrData->isScalePermstructAllocated == false ){
		throw std::logic_error( "ScalePermstruct has not been allocated by SymbolicFactorize." );
	}	
	if( ptrData->isLUstructAllocated == false ){
		throw std::logic_error( "LUstruct has not been allocated by SymbolicFactorize." );
	}	
	if( ptrData->isSuperMatrixAllocated == false ){
		throw std::logic_error( "SuperMatrix has not been allocated." );
	}	

	Int* perm_c = ptrData->ScalePermstruct.perm_c;
	NRformat_loc *Astore = (NRformat_loc *) ptrData->A.Store;
	Int* colind = Astore->colind;
	Int  nnzLocal = Astore->nnz_loc;
	// Apply column permutation to the original distributed A
	for(Int j = 0; j < nnzLocal; j++)
		colind[j] = perm_c[colind[j]];
	// Distribute Pc*Pr*diag(R)*A*diag(C)*Pc' into L and U storage.  
  // NOTE: the row permutation Pc*Pr is applied internally in the
  // distribution routine. 
	float dist_mem_use = pzdistribute(SamePattern_SameRowPerm, ptrData->A.nrow, 
			&ptrData->A, &ptrData->ScalePermstruct, NULL, &ptrData->LUstruct, ptrData->grid);

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix::Distribute  ----- 


void
SuperLUMatrix::NumericalFactorize	(  )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::NumericalFactorize");
#endif
	if( !ptrData->isLUstructAllocated ){
		throw std::logic_error( "LUstruct has not been allocated." );
	}
	// Estimate the 1-norm
	char norm[1]; *norm = '1';
	double anorm = pzlangs( norm, &ptrData->A, ptrData->grid );
  
	PStatInit(&ptrData->stat);
	pzgstrf(&ptrData->options, ptrData->A.nrow, ptrData->A.ncol, 
			anorm, &ptrData->LUstruct, ptrData->grid, &ptrData->stat, &ptrData->info); 
	PStatFree(&ptrData->stat);
	if( ptrData->info ){
		std::ostringstream msg;
		msg << "Numerical factorization error, info =  " << ptrData->info << std::endl;
		throw std::runtime_error( msg.str().c_str() );
	}
	
	// Prepare for Solve.
	ptrData->options.Fact = FACTORED;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix::NumericalFactorize  ----- 


void
SuperLUMatrix::ConvertNRlocToNC	( SuperLUMatrix& AGlobal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::ConvertNRlocToNC");
#endif
  if( !ptrData->isSuperMatrixAllocated ){
		throw std::runtime_error( "The local SuperMatrix has not been allocated." );
	}
  if( AGlobal.ptrData->isSuperMatrixAllocated ){
		throw std::runtime_error( "The global SuperMatrix has been allocated." );
	}
	// TODO make sure the two grids are the same
  
	// TODO real arithmetic
	const Int NEED_VALUE = 1;
	pzCompRow_loc_to_CompCol_global(NEED_VALUE, &ptrData->A, ptrData->grid, 
			&AGlobal.ptrData->A);

	AGlobal.ptrData->isSuperMatrixAllocated = true;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix::ConvertNRlocToNC  ----- 
void
SuperLUMatrix::MultiplyGlobalMultiVector	( NumMat<Scalar>& xGlobal, NumMat<Scalar>& bGlobal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::MultiplyGlobalMultiVector");
#endif
	char trans[1]; *trans = 'N';
	Int m = xGlobal.m();
	Int nrhs = xGlobal.n();
	bGlobal.Resize( m, nrhs );

  if( ptrData->A.Stype != SLU_NC ){
		std::ostringstream msg;
		msg << "MultiplyGlobalMultiVector requires SLU_NC matrix with type " << SLU_NC << std::endl
			<< "The matrix is of type " << ptrData->A.Stype << std::endl
			<< "Consider using ConvertNRlocToNC subroutine" << std::endl;
		throw std::runtime_error( msg.str().c_str() );
	}	
  zFillRHS_dist(trans, nrhs, (doublecomplex*)xGlobal.Data(), m, 
			&ptrData->A, (doublecomplex*) bGlobal.Data(), m);
#ifndef _RELEASE_
	PopCallStack();
#endif
 
	return ;
} 		// -----  end of method SuperLUMatrix::MultiplyGlobalMultiVector  ----- 


void
SuperLUMatrix::DistributeGlobalMultiVector	( NumMat<Scalar>& xGlobal, NumMat<Scalar>& xLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::DistributeGlobalMultiVector");
#endif
	SuperMatrix& A = ptrData->A;
  if( ptrData->A.Stype != SLU_NR_loc ){
		std::ostringstream msg;
		msg << "DistributeGlobalMultiVector requires SLU_NR_loc matrix with type " << SLU_NR_loc << std::endl
			<< "The matrix is of type " << ptrData->A.Stype << std::endl;
		throw std::runtime_error( msg.str().c_str() );
	}	

	NRformat_loc *Astore = (NRformat_loc *) ptrData->A.Store;
	
	Int numRowLocal = Astore->m_loc;
	Int firstRow    = Astore->fst_row;
	Int nrhs = xGlobal.n();
	
	xLocal.Resize( numRowLocal, nrhs );
	SetValue( xLocal, SCALAR_ZERO );
	for( Int j = 0; j < nrhs; j++ ){
		std::copy( xGlobal.VecData(j)+firstRow, xGlobal.VecData(j)+firstRow+numRowLocal,
				xLocal.VecData(j) );
	}

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix::DistributeGlobalMultiVector  ----- 


void
SuperLUMatrix::SolveDistMultiVector	( NumMat<Scalar>& bLocal, DblNumVec& berr )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::SolveDistMultiVector");
#endif
	Int nrhs = bLocal.n();
	NRformat_loc *Astore = (NRformat_loc *) ptrData->A.Store;
	Int numRowLocal = Astore->m_loc;
	Int firstRow    = Astore->fst_row;
	
	berr.Resize( nrhs );

	// TODO Real arithmetic

	PStatInit(&ptrData->stat);
	pzgssvx(&ptrData->options, &ptrData->A, &ptrData->ScalePermstruct, 
			(doublecomplex*)bLocal.Data(), numRowLocal, nrhs, ptrData->grid,
			&ptrData->LUstruct, &ptrData->SOLVEstruct, berr.Data(), 
			&ptrData->stat, &ptrData->info);
	PStatFree(&ptrData->stat); 

	if( ptrData->info ){
		std::ostringstream msg;
		msg << "Numerical solve error, info =  " << ptrData->info << std::endl;
		throw std::runtime_error( msg.str().c_str() );
	}


#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix::SolveDistMultiVector  ----- 


void
SuperLUMatrix::CheckErrorDistMultiVector	( NumMat<Scalar>& xLocal, NumMat<Scalar>& xTrueLocal )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::CheckErrorDistMultiVector");
#endif
	Int nrhs = xLocal.n();
	NRformat_loc *Astore = (NRformat_loc *) ptrData->A.Store;
	Int numRowLocal = Astore->m_loc;

	pzinf_norm_error(ptrData->grid->iam, numRowLocal, nrhs, 
			(doublecomplex*)xLocal.Data(), numRowLocal, 
			(doublecomplex*)xTrueLocal.Data(), numRowLocal, ptrData->grid);

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix::CheckErrorDistMultiVector  ----- 


void
SuperLUMatrix::LUstructToPMatrix	( PMatrix& PMloc )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::LUstructToPMatrix");
#endif
	const LocalLU_t* Llu   = ptrData->LUstruct.Llu;
	const Grid* grid       = PMloc.Grid();
	const SuperNode* super = PMloc.SuperNode();
	Int numSuper = PMloc.NumSuper();

  // L part   
#ifndef _RELEASE_
	PushCallStack("L part");
#endif
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "LUstructToPMatrix::L part" << std::endl;
#endif
	for( Int jb = 0; jb < PMloc.NumLocalBlockCol(); jb++ ){
		Int bnum = GBj( jb, grid );
		if( bnum >= numSuper ) continue;

		Int cnt = 0;                                // Count for the index in LUstruct
		Int cntval = 0;                             // Count for the nonzero values
		const Int* index = Llu->Lrowind_bc_ptr[jb];
		if( index ){ 
			// Not an empty column, start a new column then.
			std::vector<LBlock>& Lcol = PMloc.L(jb);
			Lcol.resize( index[cnt++] );
			Int lda = index[cnt++];

			for( Int iblk = 0; iblk < Lcol.size(); iblk++ ){
				LBlock& LB     = Lcol[iblk];
				LB.blockIdx    = index[cnt++];
				LB.numRow      = index[cnt++];
				LB.numCol      = super->superPtr[bnum+1] - super->superPtr[bnum];
				LB.rows        = IntNumVec( LB.numRow, true, const_cast<Int*>(&index[cnt]) );
				LB.nzval.Resize( LB.numRow, LB.numCol );   
				SetValue( LB.nzval, SCALAR_ZERO ); 
				cnt += LB.numRow;
				
				lapack::Lacpy( 'A', LB.numRow, LB.numCol, 
						(Scalar*)(Llu->Lnzval_bc_ptr[jb]+cntval), lda, 
						LB.nzval.Data(), LB.numRow );

				cntval += LB.numRow;


#if ( _DEBUGlevel_ >= 1 )
				statusOFS 
					<< "L part: bnum = " << bnum << ", numBlock = " << Lcol.size()
					<< ", blockIdx = " << LB.blockIdx
					<< ", numRow = " << LB.numRow 
					<< ", numCol = " << LB.numCol << std::endl;
#endif 

			} // for(iblk)
		}  // if(index)
	} // for(jb)
#ifndef _RELEASE_
	PopCallStack();
#endif

	// U part
#ifndef _RELEASE_
	PushCallStack("U part");
#endif
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "LUstructToPMatrix::U part" << std::endl;
#endif
	for( Int ib = 0; ib < PMloc.NumLocalBlockRow(); ib++ ){
		Int bnum = GBi( ib, grid );
		if( bnum >= numSuper ) continue;

		Int cnt = 0;                                // Count for the index in LUstruct
		Int cntval = 0;                             // Count for the nonzero values
		const Int*    index = Llu->Ufstnz_br_ptr[ib]; 
		const Scalar* pval  = reinterpret_cast<const Scalar*>(Llu->Unzval_br_ptr[ib]);
		if( index ){ 
			// Not an empty row
			// Compute the number of nonzero columns 
			std::vector<UBlock>& Urow = PMloc.U(ib);
			Urow.resize( index[cnt++] );
			cnt = BR_HEADER;

			std::vector<Int> cols;                    //Save the nonzero columns in the current block
			for(Int jblk = 0; jblk < Urow.size(); jblk++ ){
				cols.clear();
				UBlock& UB = Urow[jblk];
				UB.blockIdx = index[cnt];
				UB.numRow = super->superPtr[bnum+1] - super->superPtr[bnum];
				cnt += UB_DESCRIPTOR;
				for( Int j = FirstBlockCol( UB.blockIdx, super ); 
						 j < FirstBlockCol( UB.blockIdx+1, super ); j++ ){
					Int firstRow = index[cnt++];
					if( firstRow != FirstBlockCol( bnum+1, super ) )
						cols.push_back(j);
				}
				// Rewind the index
				cnt -= super->superPtr[UB.blockIdx+1] - super->superPtr[UB.blockIdx];

				UB.numCol = cols.size();
				UB.cols   = IntNumVec( cols.size(), true, &cols[0] );
				UB.nzval.Resize( UB.numRow, UB.numCol );
				SetValue( UB.nzval, SCALAR_ZERO );

				Int cntcol = 0;
				for( Int j = 0; 
					   j < super->superPtr[UB.blockIdx+1] - super->superPtr[UB.blockIdx]; j++ ){
					Int firstRow = index[cnt++];
					if( firstRow != FirstBlockCol( bnum+1, super ) ){
						Int tnrow = FirstBlockCol( bnum+1, super ) - firstRow;
						lapack::Lacpy( 'A', tnrow, 1, &pval[cntval], tnrow,
								&UB.nzval(firstRow - FirstBlockCol(bnum, super), cntcol),
								UB.numRow );
						cntcol ++;
						cntval += tnrow;
					}
				} // for( j )

#if ( _DEBUGlevel_ >= 1 )
				statusOFS 
					<< "U part: bnum = " << bnum << ", numBlock = " << Urow.size()
					<< ", blockIdx = " << UB.blockIdx
					<< ", numRow = " << UB.numRow 
					<< ", numCol = " << UB.numCol << std::endl;
#endif 

			} // for (jblk)



		} // if( index )

	} // for(ib)
#ifndef _RELEASE_
	PopCallStack();
#endif

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix::LUstructToPMatrix  ----- 



void
SuperLUMatrix::SymbolicToSuperNode	( SuperNode& super )
{
#ifndef _RELEASE_
	PushCallStack("SuperLUMatrix::SymbolicToSuperNode");
#endif
	Int n = ptrData->A.ncol;
	// Permutation vector
	Int *perm_c = ptrData->ScalePermstruct.perm_c;
	super.perm.Resize( n );
	std::copy( perm_c, perm_c + n, super.perm.Data() );

	// Construct the inverse permutation vector
	super.permInv.Resize( n );
	for( Int i = 0; i < n; i++ ){
		super.permInv(i) = i;
	}

	std::sort( super.permInv.Data(), super.permInv.Data() + n,
			IndexComp<IntNumVec&>(super.perm) );


  // Supernodal information

	Int *xsup    = ptrData -> LUstruct.Glu_persist -> xsup;
	Int *superno = ptrData -> LUstruct.Glu_persist -> supno;
  Int *etree = ptrData -> LUstruct.etree;

	Int numSuper = superno[ n-1 ] + 1;
  super.superPtr.Resize( numSuper + 1 );
	std::copy( xsup, xsup + numSuper + 1, super.superPtr.Data() );
	super.superIdx.Resize( n );
	std::copy( superno, superno + n, super.superIdx.Data() );
  super.etree.Resize( n );
	std::copy( etree, etree + n, super.etree.Data() );
  
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SuperLUMatrix::SymbolicToSuperNode  ----- 

} // namespace PEXSI
