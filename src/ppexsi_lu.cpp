/// @file ppexsi_lu.cpp
/// @brief Inertia subroutine for PPEXSI that directly uses SuperLU.  
///
/// This file is listed separately from other PEXSI routines because the
/// inertia count uses real arithmetic and may not be compatible with
/// the complex superlu_dist.
///
/// @author Lin Lin
/// @date 2013-07-18

#include "ppexsi.hpp"
// NOTE: IMPORTANT: Since some macros in pselinv.hpp and superlu_ddefs.h
// share the same name (such as MYROW, MYCOL), superlu_ddefs.h MUST be
// included AFTER ppexsi.hpp
#include "superlu_ddefs.h"
#include "Cnames.h"

extern "C"{
void
pzsymbfact(superlu_options_t *options, SuperMatrix *A, 
		ScalePermstruct_t *ScalePermstruct, gridinfo_t *grid,
		LUstruct_t *LUstruct, SuperLUStat_t *stat, int *numProcSymbFact,
		int *info);

void
pdsymbfact(superlu_options_t *options, SuperMatrix *A, 
		ScalePermstruct_t *ScalePermstruct, gridinfo_t *grid,
		LUstruct_t *LUstruct, SuperLUStat_t *stat, int *numProcSymbFact,
		int *info);
}


namespace PEXSI{

 


void PPEXSIData::CalculateNegativeInertiaReal( 
		const std::vector<Real>&       shiftVec, 
		std::vector<Real>&             inertiaVec,
		const DistSparseMatrix<Real>&  HMat,
		const DistSparseMatrix<Real>&  SMat,
		std::string                    ColPerm,
		Int                            numProcSymbFact
		){
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateNegativeInertiaReal");
#endif

	// *********************************************************************
	// Initialize
	// *********************************************************************

	// SuperLU Data structure

	SuperMatrix         A;
	superlu_options_t   options;
	ScalePermstruct_t   ScalePermstruct;          
	gridinfo_t*         grid;
	LUstruct_t          LUstruct;
	SOLVEstruct_t       SOLVEstruct;
	SuperLUStat_t       stat;
	Int                 info;

	// Initialize the SuperLU grid
	// Note: cannot use gridSuperLU_ because of the pImpl treatment.
	// Here it assumes that gridSuperLU_ and gridSelInv_ follow the same
	// distribution.

	grid = new gridinfo_t;
	superlu_gridinit( gridSelInv_->comm, gridSelInv_->numProcRow, 
			gridSelInv_->numProcCol, grid );

	// Get the diagonal indices for H and save it n diagIdxLocal_
	{
		Int numColLocal      = HMat.colptrLocal.m() - 1;
		Int numColLocalFirst = HMat.size / gridSelInv_->mpisize;
		Int firstCol         = gridSelInv_->mpirank * numColLocalFirst;
		
		diagIdxLocal_.clear();

		for( Int j = 0; j < numColLocal; j++ ){
			Int jcol = firstCol + j + 1;
			for( Int i = HMat.colptrLocal(j)-1; 
				 	 i < HMat.colptrLocal(j+1)-1; i++ ){
				Int irow = HMat.rowindLocal(i);
				if( irow == jcol ){
					diagIdxLocal_.push_back( i );
				}
			}
		} // for (j)
	}


	// Set SuperLU options
	set_default_options_dist(&options);

	// The default value of ColPerm uses the default value from SuperLUOptions
	options.Fact              = DOFACT;
	options.RowPerm           = NOROWPERM; // IMPORTANT for symmetric matrices
	options.IterRefine        = NOREFINE;
	options.ParSymbFact       = NO;
	options.Equil             = NO; 
	options.ReplaceTinyPivot  = NO;
	// For output information such as # of nonzeros in L and U
	// and the memory cost, set PrintStat = YES
	options.PrintStat         = NO;
	options.SolveInitialized  = NO;
	
	if ( ColPerm == "NATURAL" ){
		options.ColPerm = NATURAL;
	} 
	else if( ColPerm == "MMD_AT_PLUS_A" ){
		options.ColPerm = MMD_AT_PLUS_A;
	}
	else if( ColPerm == "METIS_AT_PLUS_A" ){
		options.ColPerm = METIS_AT_PLUS_A;
	}
	else if( ColPerm == "PARMETIS" ){
		options.ColPerm           = PARMETIS;
		options.ParSymbFact       = YES;
	}

	
	// *********************************************************************
	// Symbolic factorization.  
	// Each numPoleGroup perform independently
	// *********************************************************************
	Real timeSta, timeEnd;
	GetTime( timeSta );
	// Generate the data pattern for the SuperMatrix A
	{
		Int numRowLocal = HMat.colptrLocal.m() - 1;
		Int numRowLocalFirst = HMat.size / gridSelInv_->mpisize;
		Int firstRow = gridSelInv_->mpirank * numRowLocalFirst;

		int_t *colindLocal, *rowptrLocal;
		double *nzvalLocal;
		rowptrLocal = (int_t*)intMalloc_dist(numRowLocal+1);
		colindLocal = (int_t*)intMalloc_dist(HMat.nnzLocal); 
		nzvalLocal  = (double*)doubleMalloc_dist(HMat.nnzLocal);

		std::copy( HMat.colptrLocal.Data(), HMat.colptrLocal.Data() + HMat.colptrLocal.m(),
				rowptrLocal );
		std::copy( HMat.rowindLocal.Data(), HMat.rowindLocal.Data() + HMat.rowindLocal.m(),
				colindLocal );

		// The value is not used here
		std::copy( HMat.nzvalLocal.Data(), HMat.nzvalLocal.Data() + HMat.nzvalLocal.m(),
				nzvalLocal );

		// Important to adjust from FORTRAN convention (1 based) to C convention (0 based) indices
		for(Int i = 0; i < HMat.rowindLocal.m(); i++){
			colindLocal[i]--;
		}

		for(Int i = 0; i < HMat.colptrLocal.m(); i++){
			rowptrLocal[i]--;
		}

		// Construct the distributed matrix according to the SuperLU_DIST format
		dCreate_CompRowLoc_Matrix_dist(&A, HMat.size, HMat.size, HMat.nnzLocal, 
				numRowLocal, firstRow,
				nzvalLocal, colindLocal, rowptrLocal,
				SLU_NR_loc, SLU_D, SLU_GE);
	}
	GetTime( timeEnd );
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << "SuperMatrix is constructed." << std::endl;
	statusOFS << "Time for SuperMatrix construction is " <<
		timeEnd - timeSta << " [s]" << std::endl << std::endl;
#endif
	GetTime( timeSta );
	// Symbolic factorization
	{
		ScalePermstructInit(A.nrow, A.ncol, &ScalePermstruct);
		LUstructInit(A.nrow, A.ncol, &LUstruct);

		PStatInit(&stat);
#if ( _DEBUGlevel_ >= 1 )
		statusOFS << "Before symbfact subroutine." << std::endl;
#endif
		pdsymbfact(&options, &A, &ScalePermstruct, grid, 
				&LUstruct, &stat, &numProcSymbFact, &info);
		PStatFree(&stat);
	}
	GetTime( timeEnd );
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << "Symbolic factorization is finished." << std::endl;
#endif
	statusOFS << "Time for symbolic factorization is " <<
		timeEnd - timeSta << " [s]" << std::endl << std::endl;
	
	// NOTE: A does not need to be destroyed.  The values of A can
	// be reused in later factorization.

	Real timeShiftSta, timeShiftEnd;
	
	Int numShift = shiftVec.size();
	std::vector<Real>  inertiaVecLocal(numShift);
	inertiaVec.resize(numShift);
	for(Int l = 0; l < numShift; l++){
		inertiaVecLocal[l] = 0.0;
		inertiaVec[l]      = 0.0;
	}

	for(Int l = 0; l < numShift; l++){
		if( gridPole_->mpirank / gridPole_->numProcCol == 
				l % gridPole_->numProcRow ){
			//  MYROW( gridPole_ ) == PROW( l, gridPole_ )

			GetTime( timeShiftSta );

			statusOFS << "Shift " << l << " = " << shiftVec[l] 
				<< " processing..." << std::endl;

			// *********************************************************************
			// Data conversion
			// *********************************************************************
			
			{
				NRformat_loc *Astore = (NRformat_loc *) A.Store;
				Astore = (NRformat_loc *) A.Store;
				double* AnzvalLocal  = (double*)(Astore->nzval);
				if( SMat.size != 0 ){
					// S is not an identity matrix
					for( Int i = 0; i < HMat.nnzLocal; i++ ){
						AnzvalLocal[i] = HMat.nzvalLocal[i] - shiftVec[l] * SMat.nzvalLocal[i];
					}
				}
				else{
					// S is an identity matrix
					for( Int i = 0; i < HMat.nnzLocal; i++ ){
						AnzvalLocal[i] = HMat.nzvalLocal[i];
					}

					for( Int i = 0; i < diagIdxLocal_.size(); i++ ){
						AnzvalLocal[ diagIdxLocal_[i] ] -= shiftVec[l];
					}
				} // if (SMat.size != 0 )
			}



			// *********************************************************************
			// Factorization
			// *********************************************************************

			Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

			GetTime( timeTotalFactorizationSta );

			// Data redistribution
#if ( _DEBUGlevel_ >= 1 )
			statusOFS << "Before Distribute." << std::endl;
#endif
			{
				Int* perm_c = ScalePermstruct.perm_c;
				NRformat_loc *Astore = (NRformat_loc *) A.Store;
				Int* colind   = Astore->colind;
				Int  nnzLocal = Astore->nnz_loc;

				// Important: recompute the permuted colind
				std::copy( HMat.rowindLocal.Data(), HMat.rowindLocal.Data() +
						HMat.rowindLocal.m(), colind );

				// Important to adjust from FORTRAN convention (1 based) to C convention (0 based) indices
				for(Int i = 0; i < HMat.rowindLocal.m(); i++){
					colind[i]--;
				}

				// Apply column permutation to the original distributed A
				for(Int j = 0; j < nnzLocal; j++)
					colind[j] = perm_c[colind[j]];
				// Distribute Pc*Pr*diag(R)*A*diag(C)*Pc' into L and U storage.  
				// NOTE: the row permutation Pc*Pr is applied internally in the
				// distribution routine. 
				float dist_mem_use = pddistribute(SamePattern_SameRowPerm, A.nrow, 
						&A, &ScalePermstruct, NULL, &LUstruct, grid);
				statusOFS << "Memory usage = " << dist_mem_use << std::endl;
			}
#if ( _DEBUGlevel_ >= 1 )
			statusOFS << "After Distribute." << std::endl;
#endif

			// Numerical factorization
#if ( _DEBUGlevel_ >= 1 )
			statusOFS << "Before NumericalFactorize." << std::endl;
#endif
			{
				// Estimate the 1-norm
				char norm[1]; *norm = '1';
				double anorm = pdlangs( norm, &A, grid );

				PStatInit(&stat);
				pdgstrf(&options, A.nrow, A.ncol, 
						anorm, &LUstruct, grid, &stat, &info); 
				PStatFree(&stat);
				if( info ){
					std::ostringstream msg;
					msg << "Numerical factorization error, info =  " << info << std::endl;
					throw std::runtime_error( msg.str().c_str() );
				}
			}
#if ( _DEBUGlevel_ >= 1 )
			statusOFS << "After NumericalFactorize." << std::endl;
#endif

			GetTime( timeTotalFactorizationEnd );

			statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 

			// *********************************************************************
			// Compute inertia
			// *********************************************************************
			Real timeInertiaSta, timeInertiaEnd;
			GetTime( timeInertiaSta );

			// Compute the negative inertia of the matrix, and save it to
			// inertiaVecLocal[l]
			Real inertiaLocal = 0.0;
			{
				Int  n = A.ncol;
				Int  numSuper = LUstruct.Glu_persist -> supno[ n-1 ] + 1;
				Int *xsup = LUstruct.Glu_persist -> xsup;
				LocalLU_t* Llu   = LUstruct.Llu;
				for( Int jb = 0; jb < CEILING( numSuper, gridSelInv_->numProcCol ); jb++ ){
					Int bnum = GBj( jb, gridSelInv_ );
					if( bnum >= numSuper ) continue;

					Int cnt = 0;                                // Count for the index in LUstruct
					Int cntval = 0;                             // Count for the nonzero values
					const Int* index = Llu->Lrowind_bc_ptr[jb];
					// We only need the information from the diagonal block
					if( index ){ 
						Int  numBlock = index[cnt++];
						Int  lda      = index[cnt++];
						Int  blockIdx = index[cnt++];

						if( blockIdx == bnum ){
							// Diagonal block occurs on this processor.
							// By the convention in SuperLU, the diagonal block occurs
							// at the first block in L.
							Int  numRow   = index[cnt++];
							Int  numCol   = xsup[bnum+1] - xsup[bnum];
							// Check numRow == numCol
							if( numRow != numCol ){
								std::ostringstream msg;
								msg << "This is not a diagonal block." << std::endl 
									<< "blockIdx  = " << blockIdx << std::endl
									<< "numRow    = " << numRow 
									<< ", numCol = " << numCol << std::endl;
								throw std::runtime_error( msg.str().c_str() );
							}
							NumMat<Real>  nzval( numRow, numCol );

							lapack::Lacpy( 'A', numRow, numCol, 
									(Real*)(Llu->Lnzval_bc_ptr[jb]), lda, 
									nzval.Data(), numRow );

							for( Int i = 0; i < numRow; i++ ){
								if( nzval(i, i) < 0 )
									inertiaLocal++;
							}
						}

					}  // if(index)
				} // for(jb)
			} // Compute the inertia

			mpi::Allreduce( &inertiaLocal, &inertiaVecLocal[l], 1, 
					MPI_SUM, gridPole_->rowComm );
		} // if I am in charge of this shift

	} // for(l)

	// Collect all the negative inertia together
	mpi::Allreduce( &inertiaVecLocal[0], &inertiaVec[0], numShift, 
			MPI_SUM, gridPole_->colComm );

	// Free the data
	Destroy_LU(A.ncol, grid, &LUstruct);
	LUstructFree(&LUstruct); 
	ScalePermstructFree(&ScalePermstruct);
	Destroy_CompRowLoc_Matrix_dist(&A);
	
	superlu_gridexit(grid);
	delete grid;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::CalculateNegativeInertiaReal ----- 


} //  namespace PEXSI

