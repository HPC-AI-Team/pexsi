#include "pselinv.hpp"

namespace PEXSI{

Grid::Grid	( MPI_Comm Bcomm, int nprow, int npcol )
{
#ifndef _RELEASE_
	PushCallStack("Grid::Grid");
#endif
	Int info;
	MPI_Initialized( &info );
	if( !info ){
		throw std::logic_error( "MPI has not been initialized." );
	}
	comm = Bcomm;
  MPI_Comm_rank( comm, &mpirank );
	MPI_Comm_size( comm, &mpisize );
	if( mpisize != nprow * npcol ){
		throw std::logic_error( "mpisize != nprow * npcol." ); 
	}
	if( nprow != npcol ){
		throw std::runtime_error( "The current version of SelInv only works for square processor grids." ); }

  numProcRow = nprow;
  numProcCol = npcol;
  	
	Int myrow = mpirank / npcol;
	Int mycol = mpirank % npcol;

	MPI_Comm_split( comm, myrow, mycol, &rowComm );
	MPI_Comm_split( comm, mycol, myrow, &colComm );


#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method Grid::Grid  ----- 


Grid::~Grid	(  )
{
#ifndef _RELEASE_
	PushCallStack("Grid::~Grid");
#endif
	if( comm != MPI_COMM_NULL && comm != MPI_COMM_WORLD ){
		MPI_Comm_free( &rowComm );
		MPI_Comm_free( &colComm ); 
		MPI_Comm_free( &comm ); 
	}

#ifndef _RELEASE_
	PopCallStack();
#endif
	return ;
} 		// -----  end of method Grid::~Grid  ----- 

} // namespace PEXSI


namespace PEXSI{


PMatrix::PMatrix ( const PEXSI::Grid* g, const PEXSI::SuperNode* s ):grid_(g), super_(s)
{
#ifndef _RELEASE_
	PushCallStack("PMatrix::PMatrix");
#endif
  L_.resize( this->NumLocalBlockCol() );
	U_.resize( this->NumLocalBlockRow() );

//#if ( _DEBUGlevel_ >= 1 )
//	statusOFS << "mpirank = " << MYPROC(grid_) << std::endl;
//	statusOFS << "myrow   = " << MYROW(grid_) << std::endl; 
//	statusOFS << "mycol   = " << MYCOL(grid_) << std::endl; 
//#endif

#ifndef _RELEASE_
	PopCallStack();
#endif
	return ;
} 		// -----  end of method PMatrix::PMatrix  ----- 

PMatrix::~PMatrix() {}	


void
PMatrix::ConstructCommunicationPattern	(  )
{
#ifndef _RELEASE_
	PushCallStack("PMatrix::ConstructCommunicationPattern");
#endif
  Int numSuper = this->NumSuper();
#ifndef _RELEASE_
	PushCallStack( "Initialize the communicatino pattern" );
#endif
	isSendToDown_.Resize( numSuper, grid_->numProcRow );
	isSendToRight_.Resize( numSuper, grid_->numProcCol );
	isSendToCrossDiagonal_.Resize( numSuper );
	SetValue( isSendToDown_, false );
	SetValue( isSendToRight_, false );
	SetValue( isSendToCrossDiagonal_, false );

	isRecvFromUp_.Resize( numSuper );
	isRecvFromLeft_.Resize( numSuper );
	isRecvFromCrossDiagonal_.Resize( numSuper );
  SetValue( isRecvFromUp_, false );
	SetValue( isRecvFromLeft_, false );
	SetValue( isRecvFromCrossDiagonal_, false );

#ifndef _RELEASE_
	PopCallStack();
#endif


#ifndef _RELEASE_
	PushCallStack( "Local column communication" );
#endif
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "Local column communication" << std::endl;
#endif
	// localColBlockRowIdx stores the nonzero block indices for each local block column.
	// The nonzero block indices including contribution from both L and U.
	// Dimension: numLocalBlockCol x numNonzeroBlock
	std::vector<std::set<Int> >   localColBlockRowIdx;

	localColBlockRowIdx.resize( this->NumLocalBlockCol() );

	for( Int ksup = 0; ksup < numSuper; ksup++ ){
		// All block columns perform independently
		if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
			std::vector<Int>  tBlockRowIdx;
			tBlockRowIdx.clear();

			// L part
			std::vector<LBlock>& Lcol = this->L( LBj(ksup, grid_) );
			for( Int ib = 0; ib < Lcol.size(); ib++ ){
				tBlockRowIdx.push_back( Lcol[ib].blockIdx );
			}

			// U part
			for( Int ib = 0; ib < this->NumLocalBlockRow(); ib++ ){
				std::vector<UBlock>& Urow = this->U(ib);
				for( Int jb = 0; jb < Urow.size(); jb++ ){
					if( Urow[jb].blockIdx == ksup ){
						tBlockRowIdx.push_back( GBi( ib, grid_ ) );
					}
				}
			}

			// Communication
      std::vector<Int> tAllBlockRowIdx;
			mpi::Allgatherv( tBlockRowIdx, tAllBlockRowIdx, grid_->colComm );

			localColBlockRowIdx[LBj( ksup, grid_ )].insert(
					tAllBlockRowIdx.begin(), tAllBlockRowIdx.end() );
			
#if ( _DEBUGlevel_ >= 1 )
      statusOFS 
				<< " Column block " << ksup 
				<< " has the following nonzero block rows" << std::endl;
			for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, grid_ )].begin();
					si != localColBlockRowIdx[LBj( ksup, grid_ )].end();
					si++ ){
				statusOFS << *si << "  ";
			}
			statusOFS << std::endl; 
#endif

		} // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )
	} // for(ksup)


#ifndef _RELEASE_
	PopCallStack();
#endif


#ifndef _RELEASE_
	PushCallStack( "Local row communication" );
#endif
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "Local row communication" << std::endl;
#endif
	// localRowBlockColIdx stores the nonzero block indices for each local block row.
	// The nonzero block indices including contribution from both L and U.
	// Dimension: numLocalBlockRow x numNonzeroBlock
	std::vector<std::set<Int> >   localRowBlockColIdx;

	localRowBlockColIdx.resize( this->NumLocalBlockRow() );

	for( Int ksup = 0; ksup < numSuper; ksup++ ){
		// All block columns perform independently
		if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
			std::vector<Int>  tBlockColIdx;
			tBlockColIdx.clear();

			// U part
			std::vector<UBlock>& Urow = this->U( LBi(ksup, grid_) );
			for( Int jb = 0; jb < Urow.size(); jb++ ){
				tBlockColIdx.push_back( Urow[jb].blockIdx );
			}

			// L part
			for( Int jb = 0; jb < this->NumLocalBlockCol(); jb++ ){
				std::vector<LBlock>& Lcol = this->L(jb);
				for( Int ib = 0; ib < Lcol.size(); ib++ ){
					if( Lcol[ib].blockIdx == ksup ){
						tBlockColIdx.push_back( GBj( jb, grid_ ) );
					}
				}
			}

			// Communication
      std::vector<Int> tAllBlockColIdx;
			mpi::Allgatherv( tBlockColIdx, tAllBlockColIdx, grid_->rowComm );

			localRowBlockColIdx[LBi( ksup, grid_ )].insert(
					tAllBlockColIdx.begin(), tAllBlockColIdx.end() );
			
#if ( _DEBUGlevel_ >= 1 )
      statusOFS 
				<< " Row block " << ksup 
				<< " has the following nonzero block columns" << std::endl;
			for( std::set<Int>::iterator si = localRowBlockColIdx[LBi( ksup, grid_ )].begin();
					si != localRowBlockColIdx[LBi( ksup, grid_ )].end();
					si++ ){
				statusOFS << *si << "  ";
			}
			statusOFS << std::endl; 
#endif

		} // if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
	} // for(ksup)

#ifndef _RELEASE_
	PopCallStack();
#endif


#ifndef _RELEASE_
	PushCallStack("SendToDown / RecvFromUp");
#endif
	for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
		// Loop over all the supernodes to the right of ksup
		for( Int jsup = ksup + 1; jsup < numSuper - 1; jsup++ ){
			Int jsupLocalBlockCol = LBj( jsup, grid_ );
			Int jsupProcCol = PCOL( jsup, grid_ );
			if( MYCOL( grid_ ) == jsupProcCol ){
				// SendToDown / RecvFromUp only if (ksup, jsup) is nonzero.
				if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 ){
          for( std::set<Int>::iterator si = localColBlockRowIdx[jsupLocalBlockCol].begin();
						   si != localColBlockRowIdx[jsupLocalBlockCol].end(); si++	 ){
            Int isup = *si;
						Int isupProcRow = PROW( isup, grid_ );
						if( isup > ksup ){
							if( MYROW( grid_ ) == isupProcRow ){
								isRecvFromUp_(ksup) = true;
							}
							if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
								isSendToDown_( ksup, isupProcRow ) = true;
							}
						} // if( isup > ksup )
					} // for (si)
				} // if( localColBlockRowIdx[jsupLocalBlockCol].count( ksup ) > 0 )
			} // if( MYCOL( grid_ ) == PCOL( jsup, grid_ ) )
					
		} // for(jsup)
	} // for(ksup)
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "isSendToDown:" << isSendToDown_ << std::endl;
	statusOFS << std::endl << "isRecvFromUp:" << isRecvFromUp_ << std::endl;
#endif
 
#ifndef _RELEASE_
	PopCallStack();
#endif


#ifndef _RELEASE_
	PushCallStack("SendToRight / RecvFromLeft");
#endif
  for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
		// Loop over all the supernodes below ksup
		for( Int isup = ksup + 1; isup < numSuper - 1; isup++ ){
			Int isupLocalBlockRow = LBi( isup, grid_ );
			Int isupProcRow       = PROW( isup, grid_ );
			if( MYROW( grid_ ) == isupProcRow ){
				// SendToRight / RecvFromLeft only if (isup, ksup) is nonzero.
				if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 ){
					for( std::set<Int>::iterator si = localRowBlockColIdx[isupLocalBlockRow].begin();
							 si != localRowBlockColIdx[isupLocalBlockRow].end(); si++ ){
						Int jsup = *si;
						Int jsupProcCol = PCOL( jsup, grid_ );
						if( jsup > ksup ){
							if( MYCOL( grid_ ) == jsupProcCol ){
								isRecvFromLeft_(ksup) = true;
							}
              if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
								isSendToRight_( ksup, jsupProcCol ) = true;
							}
						}
					} // for (si)
				} // if( localRowBlockColIdx[isupLocalBlockRow].count( ksup ) > 0 )
			} // if( MYROW( grid_ ) == isupProcRow )
		} // for (isup)
  }	 // for (ksup)
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "isSendToRight:" << isSendToRight_ << std::endl;
	statusOFS << std::endl << "isRecvFromLeft:" << isRecvFromLeft_ << std::endl;
#endif

#ifndef _RELEASE_
	PopCallStack();
#endif

#ifndef _RELEASE_
	PushCallStack("SendToCrossDiagonal / RecvFromCrossDiagonal");
#endif
	for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
		if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
			for( std::set<Int>::iterator si = localColBlockRowIdx[LBj( ksup, grid_ )].begin();
					 si != localColBlockRowIdx[LBj( ksup, grid_ )].end(); si++ ){
				Int isup = *si;
				Int isupProcRow = PROW( isup, grid_ );
				if( isup > ksup && MYROW( grid_ ) == isupProcRow ){
					isSendToCrossDiagonal_( ksup ) = true;
				}
			} // for (si)
		} // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )
	} // for (ksup)
  
	for( Int ksup = 0; ksup < numSuper-1; ksup++ ){
		if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
			for( std::set<Int>::iterator si = localRowBlockColIdx[ LBi(ksup, grid_) ].begin();
					 si != localRowBlockColIdx[ LBi(ksup, grid_) ].end(); si++ ){
				Int jsup = *si;
				Int jsupProcCol = PCOL( jsup, grid_ );
				if( jsup > ksup && MYCOL(grid_) == jsupProcCol ){
					isRecvFromCrossDiagonal_[ksup] = true;
				}
			} // for (si)
		} // if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
	} // for (ksup)
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "isSendToCrossDiagonal:" << isSendToCrossDiagonal_ << std::endl;
	statusOFS << std::endl << "isRecvFromCrossDiagonal:" << isRecvFromCrossDiagonal_ << std::endl;
#endif

#ifndef _RELEASE_
	PopCallStack();
#endif


#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PMatrix::ConstructCommunicationPattern  ----- 



void
PMatrix::SelInv	(  )
{
#ifndef _RELEASE_
	PushCallStack("PMatrix::PSelInv");
#endif
  Int numSuper = this->NumSuper(); 

	this->PreSelInv();

	for( Int ksup = numSuper - 2; ksup >= 0; ksup-- ){
		this->UpdateL( ksup );

		this->UpdateD( ksup );

		this->UpdateU( ksup );
	}

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PMatrix::PSelInv  ----- 


void
PMatrix::PreSelInv	(  )
{
#ifndef _RELEASE_
	PushCallStack("PMatrix::PreSelInv");
#endif
  
  Int numSuper = this->NumSuper(); 
  	
#ifndef _RELEASE_
	PushCallStack("L(i,k) <- L(i,k) * L(k,k)^{-1}");
#endif
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "L(i,k) <- L(i,k) * L(k,k)^{-1}" << std::endl << std::endl; 
#endif
	for( Int ksup = 0; ksup < numSuper; ksup++ ){
		if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
			// Broadcast the diagonal L block
			NumMat<Scalar> nzvalLDiag;
			std::vector<LBlock>& Lcol = this->L( LBj( ksup, grid_ ) );
			if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
				nzvalLDiag = Lcol[0].nzval;
				if( nzvalLDiag.m() != SuperSize(ksup, super_) ||
						nzvalLDiag.n() != SuperSize(ksup, super_) ){
					throw std::runtime_error( "The size of the diagonal block of L is wrong." );
				}
			} // Owns the diagonal block
			else
			{
				nzvalLDiag.Resize( SuperSize(ksup, super_), SuperSize(ksup, super_) );
			}
			MPI_Bcast( (void*)nzvalLDiag.Data(), 
					sizeof(Scalar) * SuperSize(ksup, super_) * SuperSize(ksup, super_),
					MPI_BYTE, PROW( ksup, grid_ ), grid_->colComm );

			// Triangular solve
			for( Int ib = 0; ib < Lcol.size(); ib++ ){
				LBlock& LB = Lcol[ib];
				if( LB.blockIdx > ksup ){
#if ( _DEBUGlevel_ >= 1 )
					// Check the correctness of the triangular solve for the first local column
					if( LBj( ksup, grid_ ) == 0 ){
						statusOFS << "Diag   L(" << ksup << ", " << ksup << "): " << nzvalLDiag << std::endl;
						statusOFS << "Before solve L(" << LB.blockIdx << ", " << ksup << "): " << LB.nzval << std::endl;
					}
#endif
					blas::Trsm( 'R', 'L', 'N', 'U', LB.numRow, LB.numCol, SCALAR_ONE,
							nzvalLDiag.Data(), LB.numCol, LB.nzval.Data(), LB.numRow );
#if ( _DEBUGlevel_ >= 1 )
					// Check the correctness of the triangular solve for the first local column
					if( LBj( ksup, grid_ ) == 0 ){
						statusOFS << "After solve  L(" << LB.blockIdx << ", " << ksup << "): " << LB.nzval << std::endl;
					}
#endif
				}
			}
		} // if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) )
	} // for (ksup)
  

#ifndef _RELEASE_
	PopCallStack();
#endif


#ifndef _RELEASE_
	PushCallStack("U(k,i) <- L(i,k)");
#endif
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "U(k,i) <- L(i,k)" << std::endl << std::endl; 
#endif

	for( Int ksup = 0; ksup < numSuper; ksup++ ){
		Int ksupProcRow = PROW( ksup, grid_ );
		Int ksupProcCol = PCOL( ksup, grid_ );
		MPI_Status mpistatSize, mpistatContent;

		// Sender
		if( isSendToCrossDiagonal_[ksup]  &&
			  MYPROC( grid_ ) !=  PNUM( ksupProcRow, MYROW( grid_ ), grid_ ) ){
      // Pack L data
			std::stringstream sstm;
			std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
			std::vector<LBlock>&  Lcol = this->L( LBj(ksup, grid_) );
			serialize( (Int)Lcol.size(), sstm, NO_MASK );
			for( Int ib = 0; ib < Lcol.size(); ib++ ){
				if( Lcol[ib].blockIdx > ksup ){
					serialize( Lcol[ib], sstm, mask );
				}
			}
			// Send/Recv is possible here due to the one to one correspondence
			// in the case of square processor grid
			mpi::Send( sstm, PNUM( ksupProcRow, MYROW( grid_ ), grid_ ), 
					ksup + numSuper, ksup, grid_->comm );
		} // if I am a sender

		// Receiver
		if( isRecvFromCrossDiagonal_[ksup] ){
			if( PNUM( MYCOL( grid_ ), ksupProcCol, grid_ ) != MYPROC( grid_ ) ){
				std::stringstream sstm;
				mpi::Recv( sstm, PNUM( MYCOL( grid_ ), ksupProcCol, grid_ ), 
						ksup + numSuper, ksup, grid_->comm, mpistatSize, mpistatContent );
				
				// Unpack L data.  
				Int numLBlock;
				std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
				deserialize( numLBlock, sstm, NO_MASK );
				std::vector<LBlock>  Lcol(numLBlock);
				for( Int ib = 0; ib < numLBlock; ib++ ){
					deserialize( Lcol[ib], sstm, mask );
				}
				
				// Update U
				// Make sure that the size of L and the corresponding U blocks match.
				std::vector<UBlock>& Urow = this->U( LBi( ksup, grid_ ) );
				for( Int ib = 0; ib < Lcol.size(); ib++ ){
					LBlock& LB = Lcol[ib];
					if( LB.blockIdx > ksup ){
						bool isUBFound = false;
						for( Int jb = 0; jb < Urow.size(); jb++ ){
							UBlock&  UB = Urow[jb];
							if( LB.blockIdx == UB.blockIdx ){
								// Compare size
								if( LB.numRow != UB.numCol || LB.numCol != UB.numRow ){
									std::ostringstream msg;
									msg << "LB(" << LB.blockIdx << ", " << ksup << ") and UB(" 
										<< ksup << ", " << UB.blockIdx << ")	do not share the same size." << std::endl
										<< "LB: " << LB.numRow << " x " << LB.numCol << std::endl
										<< "UB: " << UB.numRow << " x " << UB.numCol << std::endl;
									throw std::runtime_error( msg.str().c_str() );
								}

								// Note that the order of the column indices of the U
								// block may not follow the order of the row indices,
								// overwrite the information in U.
								UB.cols = LB.rows;
								Transpose( LB.nzval, UB.nzval );
							
								isUBFound = true;
								break;
							} // if( LB.blockIdx == UB.blockIdx )
						} // for (jb)
						// Did not find a matching block
						if( !isUBFound ){
							std::ostringstream msg;
							msg << "LB(" << LB.blockIdx << ", " << ksup << ") did not find a matching block in U." << std::endl;
							throw std::runtime_error( msg.str().c_str() );
						}
					} // if( LB.blockIdx > ksup )
				} // for (ib)
			} // sender is not the same as receiver
			else{
				// L is obtained locally.
				std::vector<LBlock>& Lcol = this->L( LBj( ksup, grid_ ) );
				
				// Update U
				// Make sure that the size of L and the corresponding U blocks match.
				std::vector<UBlock>& Urow = this->U( LBi( ksup, grid_ ) );
				for( Int ib = 0; ib < Lcol.size(); ib++ ){
					LBlock& LB = Lcol[ib];
					if( LB.blockIdx > ksup ){
						bool isUBFound = false;
						for( Int jb = 0; jb < Urow.size(); jb++ ){
							UBlock&  UB = Urow[jb];
							if( LB.blockIdx == UB.blockIdx ){
								// Compare size
								if( LB.numRow != UB.numCol || LB.numCol != UB.numRow ){
									std::ostringstream msg;
									msg << "LB(" << LB.blockIdx << ", " << ksup << ") and UB(" 
										<< ksup << ", " << UB.blockIdx << ")	do not share the same size." << std::endl
										<< "LB: " << LB.numRow << " x " << LB.numCol << std::endl
										<< "UB: " << UB.numRow << " x " << UB.numCol << std::endl;
									throw std::runtime_error( msg.str().c_str() );
								}

								// Note that the order of the column indices of the U
								// block may not follow the order of the row indices,
								// overwrite the information in U.
								UB.cols = LB.rows;
								Transpose( LB.nzval, UB.nzval );
							
								isUBFound = true;
								break;
							} // if( LB.blockIdx == UB.blockIdx )
						} // for (jb)
						// Did not find a matching block
						if( !isUBFound ){
							std::ostringstream msg;
							msg << "LB(" << LB.blockIdx << ", " << ksup << ") did not find a matching block in U." << std::endl;
							throw std::runtime_error( msg.str().c_str() );
						}
					} // if( LB.blockIdx > ksup )
				} // for (ib)
			} // sender is the same as receiver
		} // if I am a receiver
	} // for (ksup)
  
#ifndef _RELEASE_
	PopCallStack();
#endif

#ifndef _RELEASE_
	PushCallStack("L(i,i) <- [L(k,k) * U(k,k)]^{-1} ");
#endif
#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "L(i,i) <- [L(k,k) * U(k,k)]^{-1}" << std::endl << std::endl; 
#endif

	for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
		if( MYROW( grid_ ) == PROW( ksup, grid_ ) &&
		    MYCOL( grid_ ) == PCOL( ksup, grid_ )	){
			IntNumVec ipiv( SuperSize( ksup, super_ ) );
			// Note that the pivoting vector ipiv should follow the FORTRAN
			// notation by adding the +1
			for(Int i = 0; i < SuperSize( ksup, super_ ); i++){
				ipiv[i] = i + 1;
			}
			LBlock& LB = this->L( LBj( ksup, grid_ ) )[0];
#if ( _DEBUGlevel_ >= 1 )
			// Check the correctness of the matrix inversion for the first local column
			if( LBj( ksup, grid_ ) == 0 ){
				statusOFS << "Factorized A (" << ksup << ", " << ksup << "): " << LB.nzval << std::endl;
			}
#endif
			lapack::Getri( SuperSize( ksup, super_ ), LB.nzval.Data(), 
					SuperSize( ksup, super_ ), ipiv.Data() );
#if ( _DEBUGlevel_ >= 1 )
			// Check the correctness of the matrix inversion for the first local column
			if( LBj( ksup, grid_ ) == 0 ){
				statusOFS << "Inversed   A (" << ksup << ", " << ksup << "): " << LB.nzval << std::endl;
			}
#endif
		} // if I need to inverse the diagonal block
	} // for (ksup)
  


#ifndef _RELEASE_
	PopCallStack();
#endif



#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PMatrix::PreSelInv  ----- 


void
PMatrix::UpdateL	( Int ksup )
{
#ifndef _RELEASE_
	PushCallStack("PMatrix::UpdateL");
#endif

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PMatrix::UpdateL  ----- 


void
PMatrix::UpdateD	( Int ksup )
{
#ifndef _RELEASE_
	PushCallStack("PMatrix::UpdateD");
#endif

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PMatrix::UpdateD  ----- 


void
PMatrix::UpdateU	( Int ksup  )
{
#ifndef _RELEASE_
	PushCallStack("PMatrix::UpdateU");
#endif

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PMatrix::UpdateU  ----- 

} // namespace PEXSI


//
///**********************************************************************
// * Convert SuperLU data structure to Lstruct and Ustruct required
// * by parallel selected inversion
// **********************************************************************/
//int SuperLU2SelInv(int n, LUstruct_t *LUstruct, gridinfo_t *grid, 
//		PMatrix& PMloc)
//{
//	int i, j, k, lb, ub, row, ib, ipnl, nrows, blkidx;
//	int nsupers, iam, Pc, Pr, myrow, mycol, nbc, nbr, supsize,lda;
//	int nbl, ldlnz, cnt, cntval, bnum, nnzval, nindex;
//	int nb, c, jb, nsupc, npnl, flg, ncoltot, icol, r, fstrow;
//	int *index;
//	Scalar *pval;
//	int tsize;
//	char All = 'A';
//
//
//	Glu_persist_t *glu_persist = LUstruct -> Glu_persist;
//	LocalLU_t* Llu = LUstruct-> Llu;
//
//	iam = grid->iam;
//	Pc = grid->npcol;
//	Pr = grid->nprow;
//	myrow = MYROW( iam, grid );
//	mycol = MYCOL( iam, grid );
//
//
//	PMloc._supno   = IntNumVec(n, true, glu_persist->supno);
//	PMloc._nsupers = glu_persist->supno[n-1] + 1;
//	PMloc._xsup = IntNumVec(PMloc.nsupers()+1,true,glu_persist->xsup);
//	//LLIN: The size of _xsup is nsupers + 1
//
//	PMloc._ndof = n;
//	PMloc._nbc = CEILING(PMloc.nsupers(), grid->npcol);
//	PMloc._nbr = CEILING(PMloc.nsupers(), grid->nprow);
//	PMloc._nblkl.resize(PMloc.nbc());   setvalue(PMloc._nblkl, 0);
//	PMloc._nblku.resize(PMloc.nbr());   setvalue(PMloc._nblku, 0); 
//	PMloc.L().resize(PMloc.nbc());    
//	PMloc.U().resize(PMloc.nbr());
//
//	/* L part */
//
//	for( ib = 0; ib < PMloc.nbc(); ib++ ){
//		bnum = ( (ib) * grid->npcol ) + mycol;
//		if( bnum >= PMloc.nsupers() ) continue;
//
//		cnt = 0;
//		cntval = 0;
//		index = Llu->Lrowind_bc_ptr[ib];
//		if( index ){ 
//			// Not an empty column, start a new column then.
//			vector<LBlock>& Lcol = PMloc.L(ib);
//
//			PMloc._nblkl[ib] = index[cnt++];
//			Lcol.resize(PMloc.nblkl(ib));
//			lda = index[cnt++];
//
//			for( int iblk= 0; iblk < PMloc.nblkl(ib); iblk++ ){
//				LBlock& LB = Lcol[iblk];
//				LB._blkind = index[cnt++];
//				LB._nrows = index[cnt++];
//				LB._ncols = PMloc.supsize(bnum);
//				LB._rows = IntNumVec(LB.nrows(), true, &index[cnt]);
//				LB._nzval.resize(LB.nrows(), LB.ncols());   
//				setvalue(LB._nzval, SCALAR_ZERO); 
//				cnt += LB.nrows();
//
//#if ( __DEBUGlevel >= 2 )
//				if(iam == 0){
//					fprintf(stderr,"colib = %d, nblock = %d, nrows = %d, ncols = %d\n",
//							ib, LB.blkind(), LB.nrows(), LB.ncols());
//					cerr << LB.rows() << endl;
//				}
//#endif 
//
//#ifdef _USE_COMPLEX_
//				zlacpy_(&All, &LB._nrows, &LB._ncols, 
//						(Complex*)(Llu->Lnzval_bc_ptr[ib]+cntval), 
//						&lda, LB.nzval().data(), &LB._nrows);
//#else
//				dlacpy_(&All, &LB._nrows, &LB._ncols, 
//						Llu->Lnzval_bc_ptr[ib]+cntval, 
//						&lda, LB.nzval().data(), &LB._nrows);
//#endif
//				cntval += LB.nrows();
//
//#if ( __DEBUGlevel >= 2 )
//				if(iam == 0){
//					cerr << "LB nzval follows" << endl;
//					cerr << LB.nzval() << endl;
//				}
//#endif 
//
//			} // for(iblk)
//
//
//		}  // if(index)
//	} // for(ib)
//
//
//	/* U part */
//	for( ib = 0; ib < PMloc.nbr(); ib++ ){
//		bnum = ( (ib) * grid->nprow ) + myrow;
//		if( bnum >= PMloc.nsupers() ) continue;
//
//		index = Llu->Ufstnz_br_ptr[ib];
//		pval = reinterpret_cast<Scalar*>(Llu->Unzval_br_ptr[ib]);
//		if( index ){ 
//			// Not an empty row
//			// Compute the number of nonzero columns 
//			vector<UBlock>& Urow = PMloc.U(ib);
//
//			cnt = 0;
//			cntval = 0;
//
//			PMloc._nblku[ib] = index[cnt++];
//			Urow.resize(PMloc.nblku(ib));
//			cnt = BR_HEADER;
//
//			vector<int> cols;  // save the nonzero columns in the current block
//			for(int iblk = 0; iblk < PMloc.nblku(ib); iblk++){
//				cols.clear();
//				UBlock& UB = Urow[iblk];
//				UB._blkind = index[cnt];
//				UB._nrows = PMloc.supsize(bnum);
//				cnt += UB_DESCRIPTOR;
//				for( j = PMloc.xsup(UB.blkind()); j < PMloc.xsup(UB.blkind()+1); j++ ){
//					fstrow = index[cnt++];
//					if( fstrow != PMloc.xsup(bnum+1) )
//						cols.push_back(j);
//				}
//				cnt -= PMloc.supsize(UB.blkind());  // rewind the index
//
//				UB._ncols = cols.size();
//				UB._cols = IntNumVec(cols.size(), true, &cols[0]);
//				UB._nzval.resize(UB.nrows(), UB.ncols());  
//				setvalue(UB._nzval, SCALAR_ZERO); // LLIN: IMPORTANT
//
//#if ( __DEBUGlevel >= 2 )
//				if(iam == 0){
//					fprintf(stderr,"rowib = %d, nblock = %d, nrows = %d, ncols = %d\n",
//							ib, UB.blkind(), UB.nrows(), UB.ncols());
//					cerr << UB.cols() << endl;
//				}
//#endif 
//
//
//				int tnrow, tncol;
//				int cntcol = 0;
//				for( j = 0; j < PMloc.supsize(UB.blkind()); j++ ){
//					fstrow = index[cnt++];
//					if( fstrow != PMloc.xsup(bnum+1) ){
//						tnrow = (PMloc.xsup(bnum+1)-fstrow);
//						tncol = 1;
//#ifdef _USE_COMPLEX_
//						zlacpy_(&All, &tnrow, &tncol, &pval[cntval], &tnrow, 
//								&UB._nzval(fstrow-PMloc.xsup(bnum),cntcol),
//								&UB._nrows);
//#else
//						dlacpy_(&All, &tnrow, &tncol, &pval[cntval], &tnrow, 
//								&UB._nzval(fstrow-PMloc.xsup(bnum),cntcol),
//								&UB._nrows);
//#endif
//						cntcol ++;
//						cntval += tnrow;
//					}
//				}
//
//#if ( __DEBUGlevel >= 2 )
//				if(iam == 0){
//					cerr << "UB nzval follows" << endl;
//					cerr << UB.nzval() << endl;
//				}
//#endif 
//
//			} // for(iblk) 
//		} // if(index)
//	} // for (ib)
//
//
//	return 0;
//}
//
//
///**********************************************************************
// * Construct the local elimination tree for the communication in
// * PLUSELINV 
// **********************************************************************/
//int ConstructLocalEtree(int n, gridinfo_t *grid, PMatrix& PMloc, 
//		vector<vector<int> >& localEtree)
//{
//	int mpirank;  mpirank = grid->iam; 
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	int nprow = grid->nprow,   npcol = grid->npcol;
//	int nsupers = PMloc.nsupers();
//	MPI_Comm comm = grid->comm;
//	MPI_Comm colcomm = grid->cscp.comm, rowcomm = grid->rscp.comm;
//	const int TWO = 2;
//
//	// information of localEtree to be sent to each other, each pair takes
//	// the form (isup, ksup) where isup >= ksup, assuming a structurally
//	// symmetric matrix. 
//	vector<vector<int> > blockInfoToSend(mpisize);
//
//	// Loop over all the supernodes ksup, and the column group processor
//	// locally communicate. The block info to be sent to other processors
//	// are then updated in blockInfoToSend.
//	for(int ksup = 0; ksup < nsupers; ksup++){
//		int pkrow = PROW( ksup, grid ), pkcol = PCOL(ksup, grid);
//		if( mycol == pkcol ){ 
//			vector<int> localBlockInfo;
//			vector<int> clmBlockInfo;
//
//			// Get local block information
//			int kb = LBj(ksup, grid);
//			for(int ib = 0; ib < PMloc.nblkl(kb); ib++){
//				LBlock& LB = PMloc.L(kb)[ib];
//				int isup = LB.blkind();
//				localBlockInfo.push_back(isup);
//			}
//			// All processors communicate within local column group for
//			// block information
//			int localSize = localBlockInfo.size();
//			vector<int> localSizeVec(nprow);
//			vector<int> localSizeDispls(nprow);
//			MPI_Allgather(&localSize, 1, MPI_INT, 
//					&localSizeVec[0], 1, MPI_INT, colcomm);
//			localSizeDispls[0] = 0;
//			for(int ip = 1; ip < nprow; ip++){
//				localSizeDispls[ip] = localSizeDispls[ip-1] + localSizeVec[ip-1];
//			}
//			int totalSize = localSizeDispls[nprow-1] + localSizeVec[nprow-1];
//			clmBlockInfo.resize(totalSize);
//			MPI_Allgatherv(&localBlockInfo[0], localSize, MPI_INT,
//					&clmBlockInfo[0], &localSizeVec[0], 
//					&localSizeDispls[0], MPI_INT, colcomm);
//
//#if (__DEBUGlevel >= 2)
//			if( mpirank == 0 ){
//				fprintf(stderr, "Proc = %3d, ksup = %5d, localSize = %5d, totalSize = %5d\n ", 
//						mpirank, ksup, localSize, totalSize);
//				cerr << "localSizeVec    = " << IntNumVec(nprow, false, &localSizeVec[0]) << endl;
//				cerr << "localSizeDispls = " << IntNumVec(nprow, false, &localSizeDispls[0]) << endl;
//				cerr << "clmBlockInfo    = " << IntNumVec(totalSize, false, &clmBlockInfo[0]) << endl;
//			}
//#endif
//			// Each processor owning (isup, ksup) send the all the pairs
//			// (jsup, ksup) to processors in the same row occupying (isup,
//			// jsup) 
//			for(int i = 0; i < localSize; i++){
//				int isup = localBlockInfo[i];
//				int pirow = PROW(isup, grid);
//				for(int j = 0; j < totalSize; j++){
//					int jsup = clmBlockInfo[j];
//					int pjcol = PCOL(jsup, grid);
//					int procID = PNUM(pirow, pjcol, grid);
//					for(int l = 0; l < totalSize; l++){
//						blockInfoToSend[procID].push_back(clmBlockInfo[l]); // The first element  is isup
//						blockInfoToSend[procID].push_back(ksup);            // The second element is ksup
//					}
//				} // for(j)
//			} // for(i)
//		} // if(mycol == pkcol)
//	} // for(ksup)
//
//	// Pack the sending data.
//	vector<int> sendBuf;
//	sendBuf.clear();
//	for(int ip = 0; ip < mpisize; ip++){
//		sendBuf.insert(sendBuf.end(), blockInfoToSend[ip].begin(), 
//				blockInfoToSend[ip].end());
//	}
//
//	// All processors perform All-to-All communication and receive the
//	// blockInfo outside this column group.
//	// All processors communicate with the size.
//	vector<int> sendSize(mpisize);
//	vector<int> recvSize(mpisize);
//	for(int ip = 0; ip < mpisize; ip++){
//		sendSize[ip] = blockInfoToSend[ip].size();
//	}
//	MPI_Alltoall(&sendSize[0], 1, MPI_INT,
//			&recvSize[0], 1, MPI_INT, comm);
//
//	blockInfoToSend.clear(); // save memory
//
//#if (__DEBUGlevel >= 2)
//	if( mpirank == 0 ){
//		fprintf(stderr, "Proc = %3d\n", mpirank);
//		cerr << "sendSize : " << IntNumVec(mpisize, false, &sendSize[0]) << endl;
//		cerr << "recvSize : " << IntNumVec(mpisize, false, &recvSize[0]) << endl;
//	}
//
//
//#endif
//
//
//	vector<int>  sendSizeDispls(mpisize);
//	vector<int>  recvSizeDispls(mpisize);
//	int totalSendSize, totalRecvSize;
//	sendSizeDispls[0] = 0;
//	recvSizeDispls[0] = 0;
//	for(int ip = 1; ip < mpisize; ip++){
//		sendSizeDispls[ip] = sendSizeDispls[ip-1] + sendSize[ip-1];
//		recvSizeDispls[ip] = recvSizeDispls[ip-1] + recvSize[ip-1];
//	}
//	totalSendSize = sendSizeDispls[mpisize-1] + sendSize[mpisize-1];
//	totalRecvSize = recvSizeDispls[mpisize-1] + recvSize[mpisize-1];
//
//	vector<int> recvBuf(totalRecvSize);
//
//	// Alltoall communication of the information of the elimination tree.
//	MPI_Alltoallv(&sendBuf[0], &sendSize[0], &sendSizeDispls[0], MPI_INT,
//			&recvBuf[0], &recvSize[0], &recvSizeDispls[0], MPI_INT,
//			comm);
//
//	// After receiving the information, each processor fill its own
//	// localEtree.
//	vector<set<int> > localEtreePrep;
//
//	localEtreePrep.clear();
//	localEtreePrep.resize(nsupers);
//	int isup, ksup;
//	for(int i = 0; i < totalRecvSize; i+=TWO){ // 2 comes from pair
//		isup = recvBuf[i];
//		ksup = recvBuf[i+1];
//		localEtreePrep[ksup].insert(isup); 
//	}
//
//	localEtree.clear();
//	localEtree.resize(nsupers);
//	for(int ksup = 0; ksup < nsupers; ksup++){
//		localEtree[ksup].clear();
//		localEtree[ksup].insert(localEtree[ksup].end(), 
//				localEtreePrep[ksup].begin(),
//				localEtreePrep[ksup].end());
//	}
//	return 0;
//}
//
//
///**********************************************************************
// * Dump the L matrices in MATLAB format.  All processors dump to 
// * filename_mpirank_mpisize.  
// *
// * NOTE: The output is 1-based to be the same as the MATLAB format 
// **********************************************************************/
//int PMatrix::DumpL(string filename, gridinfo_t *grid){
//	MPI_Comm comm = grid->comm;
//
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	iC( MPI_Barrier(comm) );
//	char sepfile[100];
//	sprintf(sepfile, "%s_%d_%d", filename.c_str(), mpirank, mpisize);
//
//
//
//	ofstream fid;
//	fid.open(sepfile, ios::out | ios::trunc);
//	fid.precision(15);
//
//	for(int ib = 0; ib < this->nbc(); ib++){
//		int bnum = ( (ib) * grid->npcol ) + mycol;
//		if( bnum >= this->nsupers() ) continue;
//		if( this->nblkl(ib) ){ // Not an empty column
//			for(int iblk = 0; iblk < this->nblkl(ib); iblk++){
//				LBlock& LB = this->L(ib)[iblk];
//				if(bnum == LB.blkind()){
//					// Diagonal block print the lower off diagonal
//					for(int j = 0; j < LB.ncols(); j++){
//						for(int i = 0; i < LB.nrows(); i++){
//							// NOTE: all-one diagonal is NOT outputed
//							if( i > j ){ 
//								fid << setw(12) << LB.rows(i) + 1 << setw(12) << 
//									j + this->xsup(bnum) + 1 << " " << setw(25) << scientific << 
//									LB.nzval(i,j) << endl;
//							}
//						}
//					} // for (j)
//				} 
//				else{
//					// Off-diagonal blocks
//					for(int j = 0; j < LB.ncols(); j++){
//						for(int i = 0; i < LB.nrows(); i++){
//							fid << setw(12) << LB.rows(i) + 1 << setw(12) << 
//								j + this->xsup(bnum) + 1 << " " << setw(25) << scientific << 
//								LB.nzval(i,j) << endl;
//						}
//					} // for (j)
//				} // if(bnum == LB.blkind())	
//			} // for(iblk) 
//		}  
//	} 
//	// for(ib)
//	fid.close();
//	MPI_Barrier(comm);
//	return 0;
//}
//
//
//
///**********************************************************************
// * Dump the U matrices in MATLAB format.  All processors dump to 
// * filename_mpirank_mpisize.  
// *
// * NOTE: The output is 1-based to be the same as the MATLAB format 
// **********************************************************************/
//int PMatrix::DumpU(string filename, gridinfo_t *grid){
//	MPI_Comm comm = grid->comm;
//
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	iC( MPI_Barrier(comm) );
//	char sepfile[100];
//	sprintf(sepfile, "%s_%d_%d", filename.c_str(), mpirank, mpisize);
//
//	ofstream fid;
//	fid.open(sepfile, ios::out | ios::trunc);
//	fid.precision(15);
//
//	// First output the diagonal blocks saved in L
//	for(int ib = 0; ib < this->nbc(); ib++){
//		int bnum = ( (ib) * grid->npcol ) + mycol;
//		if( bnum >= this->nsupers() ) continue;
//		if( this->nblkl(ib) ){ // Not an empty column
//			for(int iblk = 0; iblk < this->nblkl(ib); iblk++){
//				LBlock& LB = this->L(ib)[iblk];
//				// Only output for the upper-diagonal blocks for the diagonal blocks
//				if(bnum == LB.blkind()){
//					for(int j = 0; j < LB.ncols(); j++){
//						for(int i = 0; i < LB.nrows(); i++){
//							if( i <= j ){ 
//								fid << setw(12) << LB.rows(i) + 1 << setw(12) << 
//									j + this->xsup(bnum) + 1 << " " << setw(25) << scientific << 
//									LB.nzval(i,j) << endl;
//							}
//						}
//					} // for (j)
//				} 
//			} // for(iblk) 
//		}  
//	} // for(ib)
//
//	// Then output the off diagonal blocks saved in U
//	for(int ib = 0; ib < this->nbr(); ib++){
//		int bnum = ( (ib) * grid->nprow ) + myrow;
//		if( bnum >= this->nsupers() ) continue;
//		if( this->nblku(ib) ){ // Not an empty row
//			for(int iblk = 0; iblk < this->nblku(ib); iblk++){
//				UBlock& UB = this->U(ib)[iblk];
//				for(int j = 0; j < UB.ncols(); j++){
//					for(int i = 0; i < UB.nrows(); i++){
//						fid << setw(12) << i + this->xsup(bnum) + 1 << setw(12) << 
//							UB.cols(j) + 1 << " " << setw(25) << scientific << 
//							UB.nzval(i,j) << endl;
//					}
//				} // for (j)
//			} // for(iblk) 
//		}  
//	} 
//	// for(ib)
//	fid.close();
//	MPI_Barrier(comm);
//	return 0;
//}
//
///**********************************************************************
// * Dump one block of L in matlab format
// *
// * NOTE: 
// *
// * o The diagonal block should be dumped by this subroutine
// *
// * o The index is the local index of this block rather than global index
// **********************************************************************/
//int PMatrix::DumpLBlock(int isup, int jsup, string filename, gridinfo_t *grid){
//	MPI_Comm comm = grid->comm;
//
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	iC( MPI_Barrier(comm) );
//
//	iA(isup >= 0 && isup < this->nsupers());
//	iA(jsup >= 0 && jsup < this->nsupers());
//
//
//	for(int ib = 0; ib < this->nbc(); ib++){
//		int bnum = ( (ib) * grid->npcol ) + mycol;
//		if( bnum >= this->nsupers() ) continue;
//		if( bnum != jsup ) continue;
//
//		if( this->nblkl(ib) ){ // Not an empty column
//			int iblk;
//			for(iblk = 0; iblk < this->nblkl(ib); iblk++){
//				LBlock& LB = this->L(ib)[iblk];
//				if( LB.blkind() != isup ) continue;
//
//				// Output the (isup, jsup) block
//				ofstream fid;
//				fid.open(filename.c_str(), ios::out | ios::trunc);
//				fid.precision(15);
//
//				cerr << "processor (" << Index2(myrow,mycol) << 
//					") outputs block (" << Index2(isup,jsup) << ") in MATLAB format" << endl;
//				cerr << "Row Range: " << this->xsup(LB.blkind()) + 1 << " -- " <<
//					this->xsup(LB.blkind()+1) << " , size = " << this->supsize(LB.blkind()) << endl;
//				cerr << "Col Range: " << this->xsup(bnum) + 1<< " -- " <<
//					this->xsup(bnum+1) << " , size = " << this->supsize(bnum) << endl;
//				int rowshift = this->xsup(LB.blkind());
//				for(int j = 0; j < LB.ncols(); j++){
//					for(int i = 0; i < LB.nrows(); i++){
//						fid << setw(12) << LB.rows(i)-rowshift+1 
//							<< setw(12) << j+1  << setw(25) << scientific << 
//							LB.nzval(i,j) << endl;
//					}
//				} // for (j)
//				fid.close();
//			} // for(iblk) 
//		}
//	} // for(ib)
//	MPI_Barrier(comm);
//	return 0;
//}
//
///**********************************************************************
// * Dump one block of U in matlab format
// *
// * NOTE: 
// *
// * o The diagonal block is NOT dumped by this subroutine
// *
// * o The index is the local index of this block rather than global index
// **********************************************************************/
//int PMatrix::DumpUBlock(int isup, int jsup, string filename, gridinfo_t *grid){
//	MPI_Comm comm = grid->comm;
//
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	iC( MPI_Barrier(comm) );
//
//	iA(isup >= 0 && isup < this->nsupers());
//	iA(jsup >= 0 && jsup < this->nsupers());
//
//	int prow = PROW(isup, grid), pcol = PCOL(jsup, grid);
//
//	if( myrow == prow && mycol == pcol ){
//		int ib = LBi(isup, grid);
//		if(this->nblku(ib)){ // Not an empty row
//			for(int jb = 0; jb < this->nblku(ib); jb++){
//				UBlock& UB = this->U(ib)[jb];
//				if( UB.blkind() != jsup ) continue;
//				// Output the (isup, jsup) block
//				ofstream fid;
//				fid.open(filename.c_str(), ios::out | ios::trunc);
//				fid.precision(15);
//
//				cerr << "processor (" << Index2(myrow,mycol) << 
//					") outputs block (" << Index2(isup,jsup) << ") in MATLAB format" << endl;
//				cerr << "Row Range: " << this->xsup(isup) + 1<< " -- " <<
//					this->xsup(isup+1) << " , size = " << this->supsize(isup) << endl;
//				cerr << "Col Range: " << this->xsup(jsup) + 1 << " -- " <<
//					this->xsup(jsup+1) << " , size = " << this->supsize(jsup) << endl;
//				int colshift = this->xsup(jsup);
//				for(int j = 0; j < UB.ncols(); j++){
//					for(int i = 0; i < UB.nrows(); i++){
//						fid << setw(12) << i+1  << setw(12) << UB.cols(j)-colshift+1  
//							<< setw(25) << scientific << 
//							UB.nzval(i,j) << endl;
//					}
//				} // for (j)
//				fid.close();
//			} // for (jb)
//		}
//	} // if( myrow == prow && mycol == pcol )
//
//	MPI_Barrier(comm);
//	return 0;
//}
//
//
//
///**********************************************************************
// * Main function for carrying out parallel selected inversion 
// **********************************************************************/
//int PLUSelInv(gridinfo_t *grid, PMatrix& PMloc, vector<vector<int> >& localEtree)
//{
//	//FIXME
//	timettl1 = 0;  timettl2 = 0;  timettl3 = 0; timettl4 = 0; 
//
//	int n = PMloc.ndof();
//	int nsupers = PMloc.nsupers();
//	double t0, t1;
//
//	//FIXME
//	for(int ksup = nsupers-1; ksup >= 0; ksup--){
//		t0 = MPI_Wtime();
//
//#if (__PRNTlevel >= 2 )
//		if(grid->iam == 0){
//			cout << "ksup = " << ksup << endl;
//		}
//#endif
//		/***********************************
//		 * Only Inversion of the L and U factors in the diagonal block
//		 ************************************/
//		if( ksup == nsupers-1 ){
//			iC(DiagTri(PMloc, ksup, grid));
//
//#if (__PRNTlevel >= 2 )
//			MPI_Barrier(grid->comm);
//			if(grid->iam == 0){
//				cout << "All processors passed DiagTri" << endl;
//			}
//#endif
//			continue;
//		}
//
//
//		/***********************************
//		 * Scale PM(i,k), i>k by -L^{-1}_{k,k}
//		 * Scale PM(k,j), j>k by -U^{-1}_{k,k}
//		 ************************************/
//		iC(ScaleLinv(PMloc, ksup, grid));
//
//#if (__PRNTlevel >= 2 )
//		MPI_Barrier(grid->comm);
//		if(grid->iam == 0){
//			cout << "All processors passed ScaleLinv" << endl;
//		}
//#endif
//
//		iC(ScaleUinv(PMloc, ksup, grid));
//
//#if (__PRNTlevel >= 2 )
//		MPI_Barrier(grid->comm);
//		if(grid->iam == 0){
//			cout << "All processors passed ScaleUinv" << endl;
//		}
//#endif
//
//
//#if ( __DEBUGlevel >= 2 )
//		iC(PMloc.DumpLBlock(ksup+1, ksup, "L", grid));
//#endif
//
//#if ( __DEBUGlevel >= 2 )
//		iC(PMloc.DumpUBlock(ksup, ksup+1, "U", grid));
//#endif
//
//
//		/***********************************
//		 * PL(j) = -L(i,k) L^-1_kk
//		 * PLT(i) = S^-1(i,j) PL(j)
//		 ************************************/
//
//		iC(SinvPL(PMloc, localEtree, ksup, grid));
//
//#if (__PRNTlevel >= 2 )
//		MPI_Barrier(grid->comm);
//		if(grid->iam == 0){
//			cout << "All processors passed SinvPL" << endl;
//		}
//#endif
//
//
//#if ( __DEBUGlevel >= 2 )
//		iC(PMloc.DumpLBlock(ksup+1, ksup, "L", grid));
//#endif
//
//		/***********************************
//		 * Invert diagonal block
//		 ************************************/
//		iC(DiagTri(PMloc, ksup, grid));
//
//#if (__PRNTlevel >= 2 )
//		MPI_Barrier(grid->comm);
//		if(grid->iam == 0){
//			cout << "All processors passed DiagTri" << endl;
//		}
//#endif
//
//
//
//		/***********************************
//		 * Update the diagonal block by
//		 *   U_kk^-1 U_ik S^-1_ij L_jk L^-1_kk
//		 ************************************/
//
//		iC(DiagInnerProd(PMloc, localEtree, ksup, grid));
//
//#if (__PRNTlevel >= 2 )
//		MPI_Barrier(grid->comm);
//		if(grid->iam == 0){
//			cout << "All processors passed DiagInnerProd" << endl;
//		}
//#endif
//
//
//#if ( __DEBUGlevel >= 1 )
//		if(ksup == 0 ) iC(PMloc.DumpLBlock(ksup, ksup, "L", grid));
//#endif
//
//
//
//		/***********************************
//		 * PU(i) = -U^-1(k,k) U(i,k) 
//		 * PLU(j) = sum_i PU(i) S^-1(i,j) 
//		 ************************************/
//		iC(SinvPU(PMloc, localEtree, ksup, grid));
//
//#if (__PRNTlevel >= 2 )
//		MPI_Barrier(grid->comm);
//		if(grid->iam == 0){
//			cout << "All processors passed SinvPU" << endl;
//		}
//#endif
//
//
//#if ( __DEBUGlevel >= 2 )
//		iC(PMloc.DumpUBlock(ksup, ksup+1, "U", grid));
//#endif
//
//
//		t1 = MPI_Wtime();
//		// LL: Not timing this first
//#if (__PRNTlevel >= 2)
//		if(grid->iam == 0) {
//			fprintf(stderr, "Time for solving ksup [%6d] is %10.3f s\n", ksup,
//					t1-t0);
//		}
//#endif
//
//	} // for (ksup)
//
//
//	//FIXME
//	//  std::cout << "timettl1 = " << timettl1 << endl;
//	//  std::cout << "timettl2 = " << timettl2 << endl;
//	//  std::cout << "timettl3 = " << timettl3 << endl;
//	//  std::cout << "timettl4 = " << timettl4 << endl;
//	return 0;
//}
//
///**********************************************************************
// * Triangular inversion of the diagonal blocks (ksup, ksup)
// **********************************************************************/
//int DiagTri(PMatrix& PMloc, int ksup, gridinfo_t *grid){
//	MPI_Comm comm = grid->comm;
//
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	int pkrow = PROW( ksup, grid ), pkcol = PCOL(ksup, grid);
//
//	iC( MPI_Barrier(comm) );
//
//	iA(ksup >= 0 && ksup < PMloc.nsupers());
//
//	if( myrow == pkrow && mycol == pkcol ){
//		int kb = LBj(ksup,grid);
//		iA( PMloc.nblkl(kb) ); // Not an empty column
//		LBlock& LB = PMloc.L(kb)[0]; // Diagonal block
//
//		{
//			int info;
//			int tnrow = PMloc.supsize(LB.blkind());
//			int lwork = -1;
//			NumVec<Scalar> work(1);  
//			IntNumVec ipiv(tnrow);
//
//			for(int i = 0; i < tnrow; i++) ipiv[i] = i+1; //LLIN: IMPORTANT
//
//#ifdef _USE_COMPLEX_
//			zgetri_(&tnrow, LB.nzval().data(), &tnrow, ipiv.data(),
//					work.data(), &lwork, &info);
//			lwork = (int)work[0].real();
//			work.resize(lwork);
//			zgetri_(&tnrow, LB.nzval().data(), &tnrow, ipiv.data(),
//					work.data(), &lwork, &info);
//#else
//			dgetri_(&tnrow, LB.nzval().data(), &tnrow, ipiv.data(),
//					work.data(), &lwork, &info);
//			lwork = (int)work[0];
//			work.resize(lwork);
//			dgetri_(&tnrow, LB.nzval().data(), &tnrow, ipiv.data(),
//					work.data(), &lwork, &info);
//#endif
//			iC(info);
//
//		} 
//	}
//
//	iC(MPI_Barrier(comm));
//
//	return 0;
//}
//
//
///**********************************************************************
// * Scale PM(i,k), i>k by -L^{-1}_{k,k}
// **********************************************************************/
//int ScaleLinv(PMatrix& PMloc, int ksup, gridinfo_t *grid)
//{
//
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	iC( MPI_Barrier(grid->comm) );
//
//	iA(ksup >= 0 && ksup < PMloc.nsupers()-1);
//
//	MPI_Comm colcomm = grid->cscp.comm, rowcomm = grid->rscp.comm;
//	int prow = PROW( ksup, grid ), pcol = PCOL(ksup, grid);
//
//	if( mycol == pcol ){
//		// Each processor in the column of pcol receive L(ksup,ksup) block
//		// from (prow,pcol)
//
//		LBlock LBbuf; 
//		int kb = LBj(ksup, grid);
//
//		if( myrow == prow ){
//			LBbuf = PMloc._L[kb][0];
//		}
//		iC(BcastLBlock(LBbuf, myrow, prow, colcomm));
//
//		int nrows, ncols;
//		for(int ib = 0; ib < PMloc.nblkl(kb); ib++){
//			LBlock& LB = PMloc.L(kb)[ib];
//			int isup = LB.blkind();
//			if( isup > ksup ){
//				nrows = LB.nrows();
//				ncols = LB.ncols();
//#ifdef _USE_COMPLEX_
//				Complex Z_MONE = -Z_ONE;
//				ztrsm_(&FROMRIGHT, &LOWER, &NOTRAN, &UDIAG, &nrows,&ncols,
//						&Z_MONE, LBbuf.nzval().data(), &ncols, LB.nzval().data(), 
//						&nrows);
//#else
//				Real   D_MONE = -D_ONE;
//				dtrsm_(&FROMRIGHT, &LOWER, &NOTRAN, &UDIAG, &nrows,&ncols,
//						&D_MONE, LBbuf.nzval().data(), &ncols, LB.nzval().data(), 
//						&nrows);
//#endif
//			}
//		}
//	} // if( mycol == pcol )
//
//	iC( MPI_Barrier(grid->comm) );
//
//	return 0;
//}
//
///**********************************************************************
// * Scale PM(k,j), j>k by -U^{-1}_{k,k}
// **********************************************************************/
//int ScaleUinv(PMatrix& PMloc, int ksup, gridinfo_t *grid)
//{
//
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	iC( MPI_Barrier(grid->comm) );
//
//	iA(ksup >= 0 && ksup < PMloc.nsupers()-1);
//
//	MPI_Comm colcomm = grid->cscp.comm, rowcomm = grid->rscp.comm;
//	int pkrow = PROW( ksup, grid ), pkcol = PCOL(ksup, grid);
//
//	if( myrow == pkrow ){
//		// Each processor in the row of pkrow receive U(ksup,ksup) block
//		// from (pkrow,pkcol)
//		LBlock Diagbuf; 
//
//		if( mycol == pkcol ){
//			int kb = LBj(ksup, grid);
//			Diagbuf = PMloc.L(kb)[0]; // Diagonal block saved in L
//		}
//		iC(BcastLBlock(Diagbuf, mycol, pkcol, rowcomm));
//
//		// Triangular solve
//		int nrows, ncols;
//
//		int kb = LBi(ksup,grid); // LLIN: IMPORTANT, different from L
//		for(int jb = 0; jb < PMloc.nblku(kb); jb++){
//			UBlock& UB = PMloc.U(kb)[jb];
//			int jsup = UB.blkind();
//			if( jsup > ksup ){
//				nrows = UB.nrows();
//				ncols = UB.ncols();
//		
//#ifdef _USE_COMPLEX_
//				Complex Z_MONE = -Z_ONE;
//				ztrsm_(&FROMLEFT, &UPPER, &NOTRAN, &NOUDIAG, &nrows,&ncols,
//						&Z_MONE, Diagbuf.nzval().data(), &nrows, UB.nzval().data(), 
//						&nrows);
//#else
//				Real   D_MONE = -D_ONE;
//				dtrsm_(&FROMLEFT, &UPPER, &NOTRAN, &NOUDIAG, &nrows,&ncols,
//						&D_MONE, Diagbuf.nzval().data(), &nrows, UB.nzval().data(), 
//						&nrows);
//#endif
//			}
//		}
//	} // if( myrow == pkrow)
//
//	iC( MPI_Barrier(grid->comm) );
//
//	return 0;
//}
//
//
///**********************************************************************
// * Broadcast an LBlock to a row or column sub communicator
// *
// **********************************************************************/
//int BcastLBlock(LBlock& LB, int mykey, int srckey, MPI_Comm comm){ // new version
//	vector<char> str;
//	int sz = 0;
//	stringstream sstm;
//
//	if( mykey == srckey ){
//		vector<int> all(LBlock_Number,1);
//		serialize(LB, sstm, all);
//		const string& sstr = sstm.str();
//		sz = sstr.length();
//		MPI_Bcast(&sz, 1, MPI_INT, srckey, comm);
//		MPI_Bcast((void*)sstr.c_str(), sz, MPI_BYTE, srckey, comm);
//	}
//	else{
//		MPI_Bcast(&sz, 1, MPI_INT, srckey, comm);
//		str.resize(sz);
//		MPI_Bcast(&str[0], sz, MPI_CHAR, srckey, comm);
//
//		// Deseralize
//		vector<int> all(LBlock_Number,1);
//		sstm.write(&str[0], sz);
//		deserialize(LB, sstm, all);
//	}
//
//
//#if ( __DEBUGlevel >= 2 )
//	cerr << LB._nzval << endl;
//#endif
//
//	return 0;
//} 
//
//
///**********************************************************************
// * Broadcast an UBlock to a row or column sub communicator
// **********************************************************************/
//int BcastUBlock(UBlock& UB, int mykey, int srckey, MPI_Comm comm){ // new version
//	vector<char> str;
//	int sz = 0;
//	stringstream sstm;
//
//	if( mykey == srckey ){
//		vector<int> all(UBlock_Number,1);
//		serialize(UB, sstm, all);
//		const string& sstr = sstm.str();
//		sz = sstr.length();
//		MPI_Bcast(&sz, 1, MPI_INT, srckey, comm);
//		MPI_Bcast((void*)sstr.c_str(), sz, MPI_BYTE, srckey, comm);
//	}
//	else{
//		MPI_Bcast(&sz, 1, MPI_INT, srckey, comm);
//		str.resize(sz);
//		MPI_Bcast((void*)&str[0], sz, MPI_BYTE, srckey, comm);
//
//		// Deseralize
//		vector<int> all(UBlock_Number,1);
//		sstm.write(&str[0], sz);
//		deserialize(UB, sstm, all);
//#if ( __DEBUGlevel >= 2 )
//		cerr << sz << endl;
//		cerr << UB._nzval << endl;
//#endif
//	}
//
//	return 0;
//} 
//
//
///**********************************************************************
// * PLT(i) = S^-1(i,j) PL(j)
// *
// * Build PLT in each processor from
// * (ksup+1:nsupers,ksup+1:nsupers)
// *
// * Build PL(i) = -L(i,k) L^-1(k,k) in each processor from
// * (ksup+1:nsupers,ksup+1:nsupers)
// *
// * Local matrix-matrix multiplication
// *
// * Reduce the PLT to processors (ksup+1:nsupers, ksup)
// **********************************************************************/
//int SinvPL(PMatrix& PMloc, vector<vector<int> >& localEtree, int ksup, gridinfo_t* grid){
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	int nprow = grid->nprow,   npcol = grid->npcol;
//	int nsupers = PMloc.nsupers();
//	const int NULLID = -1;
//	const int NDATA = 1;
//	enum {SUPID};
//
//	double t0, t1;
//
//	map<int, LBlock> PL;
//	map<int, LBlock> PLT;
//
//	PL.clear();
//	PLT.clear();
//
//	iA(ksup >= 0 && ksup < nsupers-1);
//	int pkrow = PROW(ksup, grid), pkcol = PCOL(ksup, grid);
//	int kb = LBj(ksup, grid);
//
//	MPI_Comm comm = grid->comm;
//	MPI_Comm colcomm = grid->cscp.comm, rowcomm = grid->rscp.comm;
//
//	// Use localEtree for communication
//	if(1)
//	{
//		vector<int>& nodesksup = localEtree[ksup];
//		int sizeksup = nodesksup.size();
//
//		// Only processors with nonvanishing setksup participates SinvPL.
//		if(sizeksup == 0) return 0;
//
//
//		// Each processor knows the row and column processor id arrays that
//		// is restricted SinvPL.
//		vector<int> prowRestricted;  
//		vector<int> pcolRestricted;  
//		{
//			set<int> prowtmp;
//			set<int> pcoltmp;
//
//			int isup, pirow, picol;
//			for(int i = 0; i < sizeksup; i++){
//				isup  = nodesksup[i];
//				pirow = PROW(isup, grid);
//				picol = PCOL(isup, grid);
//				prowtmp.insert(pirow);
//				pcoltmp.insert(picol);
//			}
//
//			prowRestricted.clear();
//			prowRestricted.insert(prowRestricted.end(), prowtmp.begin(),
//					prowtmp.end());
//			pcolRestricted.clear();
//			pcolRestricted.insert(pcolRestricted.end(), pcoltmp.begin(),
//					pcoltmp.end());
//		}
//		int nprowRestricted = prowRestricted.size();
//		int npcolRestricted = pcolRestricted.size();
//#if (__DEBUGlevel >= 2 )
//		if(mpirank == 0){
//			cerr << "ksup = " << ksup << endl;
//			cerr << "prowRestricted = " << 
//				IntNumVec(nprowRestricted, false, &prowRestricted[0]);
//			cerr << "pcolRestricted = " << 
//				IntNumVec(npcolRestricted, false, &pcolRestricted[0]);
//		}
//#endif
//
//#if (__PRNTlevel >= 3 )
//		if(grid->iam == 0){
//			cout << "SinvPL: initiation passed" << endl;
//		}
//#endif
//
//		// Setup PLT, in preparation for PLT = PM * PL
//		{
//			MPI_Request*    reqs   = new MPI_Request[npcolRestricted];   iA( reqs != NULL ); 
//			MPI_Status *    stats  = new MPI_Status [npcolRestricted];   iA( stats != NULL );
//			vector<int> mask(LBlock_Number,1); // Communication for all
//
//			for(int inode = 0; inode < sizeksup; inode++){
//				int isup  = nodesksup[inode];
//				int pirow = PROW(isup, grid);
//				if( isup > ksup && myrow == pirow ){
//
//					// Sender/Receiver only has at most one copy of the data
//					stringstream   sstm;
//					int            sizeMsg;
//					vector<char>   strchr;
//
//					for(int ip = 0; ip < npcolRestricted; ip++){
//						reqs[ip]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					// Sender
//					if( mycol == pkcol ){
//						LBlock *ptrLB = NULL;
//						for(int ib = 0; ib < PMloc.nblkl(kb); ib++){
//							if( PMloc.L(kb)[ib].blkind() == isup ){
//								ptrLB = &PMloc.L(kb)[ib];
//							}
//						}
//						iA( ptrLB != NULL );
//						PLT[isup] = *ptrLB; // setup local PLT
//						serialize(*ptrLB, sstm, mask);
//						sizeMsg = (sstm.str()).length();
//
//					} 
//
//					// Size information
//					int sendID = PNUM(pirow, pkcol, grid);
//					for(int jp = 0; jp < npcolRestricted; jp++){
//						int pjcol  = pcolRestricted[jp];
//						int recvID = PNUM(pirow, pjcol, grid);
//						if( sendID != recvID ){
//							if( mpirank == sendID ){
//								MPI_Isend(&sizeMsg, 1, MPI_INT, pjcol, pjcol, rowcomm, &(reqs[jp]));
//							}
//							if( mpirank == recvID ){
//								MPI_Irecv(&sizeMsg, 1, MPI_INT, pkcol, pjcol, rowcomm, &(reqs[jp])); 
//							}
//						}
//					}
//					MPI_Waitall(npcolRestricted, &(reqs[0]), &(stats[0]));
//
//
//					for(int ip = 0; ip < npcolRestricted; ip++){
//						reqs[ip]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					for(int jp = 0; jp < npcolRestricted; jp++){
//						int pjcol  = pcolRestricted[jp];
//						int recvID = PNUM(pirow, pjcol, grid);
//						if( sendID != recvID ){
//							if( mpirank == sendID ){
//								const string& sstr = sstm.str();
//								MPI_Isend((void*)sstr.c_str(), sizeMsg, MPI_BYTE, 
//										pjcol, pjcol, rowcomm, &(reqs[jp]));
//							}
//							if( mpirank == recvID ){
//								strchr.clear();
//								strchr.resize(sizeMsg);
//								MPI_Irecv((void*)(&strchr[0]), sizeMsg, MPI_BYTE, 
//										pkcol, pjcol, rowcomm, &(reqs[jp])); 
//							}
//						}
//					}
//					MPI_Waitall(npcolRestricted, &(reqs[0]), &(stats[0]));
//
//					// Deserialize into PLT
//					if( mpirank != sendID ){
//						sstm.write(&strchr[0], sizeMsg);
//						iA(PLT.count(isup) == 0);  // PLT should be empty
//						deserialize(PLT[isup], sstm, mask);
//#if (__DEBUGlevel >= 2 )
//						LBlock& LB = PLT[isup];
//						fprintf(stderr, "mpirank = %5d\n", mpirank);
//						fprintf(stderr, "LB(%5d,%5d) \n",  LB.blkind(), ksup);
//						fprintf(stderr, "nrows = %5d, ncols = %5d\n",
//								LB.nrows(), LB.ncols());
//#endif
//					}
//				} 
//
//			}// for (inode)
//			delete[] reqs;    reqs = NULL;
//			delete[] stats;   stats = NULL;
//		}
//
//#if (__DEBUGlevel >= 2 )
//		if(mpirank == 0){
//			for(map<int,LBlock>::iterator mi = PLT.begin();
//					mi != PLT.end(); mi++){
//				cerr << (*mi).second.blkind() << endl;
//			}
//		}
//#endif
//
//#if (__PRNTlevel >= 3 )
//		if(grid->iam == 0){
//			cout << "SinvPL: PLT distribution passed" << endl;
//		}
//#endif
//
//
//		// For all processors in the column processor group of ksup, if
//		// LBlock(isup, ksup) is nonempty, setup PL locally. Then send the
//		// information of PL (all information) to all processors owning
//		// (jsup, isup), if (jsup, ksup) is in the localEtree. (Do not send
//		// to itself)
//
//		// For all processors owning (isup, jsup) with isup > ksup, jsup >
//		// ksup, receive from (jsup, ksup) if (jsup, ksup) is in the
//		// localEtree, and setup the local PL. Can overlap communication
//		// with computation. (Do not receive from itself)
//		{
//			MPI_Request*    reqs   = new MPI_Request[nprowRestricted];   iA( reqs  != NULL ); 
//			MPI_Status *    stats  = new MPI_Status [nprowRestricted];   iA( stats != NULL );
//			vector<int> mask(LBlock_Number,1);  // all data is needed for PL.
//
//			for(int inode = 0; inode < sizeksup; inode++){
//				int isup  = nodesksup[inode];
//				int pirow = PROW(isup, grid);
//				int pjcol = PCOL(isup, grid);
//				// Sender
//				// The diagonal processor sends PL to column processors
//				int sendID = PNUM(pirow, pjcol, grid); 
//				if( isup > ksup && mpirank == sendID){
//
//					stringstream   sstm;
//					int            sizeMsg;
//
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						reqs[ip]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					PL[isup] = PLT[isup]; // Should have received from PM(isup,ksup)
//					serialize(PL[isup], sstm, mask);
//					sizeMsg = (sstm.str()).length();
//
//					// Size information
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						int pirowRecv  = prowRestricted[ip];
//						int recvID = PNUM(pirowRecv, pjcol, grid);
//						if( sendID != recvID ){
//							MPI_Isend(&sizeMsg, 1, MPI_INT, pirowRecv, pirowRecv, colcomm, &(reqs[ip]));
//						}
//					}
//					MPI_Waitall(nprowRestricted, &(reqs[0]), &(stats[0]));
//
//					// Actual data
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						reqs[ip]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						int pirowRecv  = prowRestricted[ip];
//						int recvID = PNUM(pirowRecv, pjcol, grid);
//						if( sendID != recvID ){
//							const string& sstr = sstm.str();
//							MPI_Isend((void*)sstr.c_str(), sizeMsg, MPI_BYTE, 
//									pirowRecv, pirowRecv, colcomm, &(reqs[ip]));
//						}
//					}
//					MPI_Waitall(nprowRestricted, &(reqs[0]), &(stats[0]));
//				}  // if( mpirank == sendID )
//
//				// Receiver
//				if( isup > ksup && mycol == pjcol ){
//					stringstream   sstm;
//					int            sizeMsg;
//					vector<char>   strchr;
//
//					// Size information
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						reqs[ip]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						int pirowRecv  = prowRestricted[ip];
//						int recvID = PNUM(pirowRecv, pjcol, grid);
//						if( sendID != recvID && mpirank == recvID ){
//							MPI_Irecv(&sizeMsg, 1, MPI_INT, pirow, pirowRecv, colcomm, &(reqs[ip]));
//						}
//					}
//					MPI_Waitall(nprowRestricted, &(reqs[0]), &(stats[0]));
//
//					// Actual data
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						reqs[ip]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						int pirowRecv  = prowRestricted[ip];
//						int recvID     = PNUM(pirowRecv, pjcol, grid);
//						if( sendID != recvID && mpirank == recvID ){
//							strchr.clear();
//							strchr.resize(sizeMsg);
//							MPI_Irecv((void*)(&strchr[0]), sizeMsg, MPI_BYTE,
//									pirow, pirowRecv, colcomm, &(reqs[ip]));
//						}
//					}
//					MPI_Waitall(nprowRestricted, &(reqs[0]), &(stats[0]));
//
//					// Deserialize into PL
//					if( mpirank != sendID ){
//						sstm.write(&strchr[0], sizeMsg);
//						iA(PL.count(isup) == 0);  // PL should be empty
//						deserialize(PL[isup], sstm, mask);
//#if (__DEBUGlevel >= 2 )
//						LBlock& LB = PL[isup];
//						fprintf(stderr, "mpirank = %5d\n", mpirank);
//						fprintf(stderr, "LB(%5d,%5d) \n",  LB.blkind(), ksup);
//						fprintf(stderr, "nrows = %5d, ncols = %5d\n",
//								LB.nrows(), LB.ncols());
//#endif
//					} 
//
//				}
//			}// for (inode)
//			delete[] reqs;    reqs = NULL;
//			delete[] stats;   stats = NULL;
//		}
//
//#if (__DEBUGlevel >= 2 )
//		if(mpirank == 0){
//			cerr << "PL distribution passed" << endl;
//			for(map<int,LBlock>::iterator mi = PL.begin();
//					mi != PL.end(); mi++){
//				cerr << (*mi).second.blkind() << endl;
//			}
//		}
//#endif
//
//#if (__PRNTlevel >= 3 )
//		if(grid->iam == 0){
//			cout << "SinvPL: PL distribution passed" << endl;
//		}
//#endif
//
//
//		// Perform matrix vector multiplication 
//		{
//			// For each processor, clean the PLT for safety. 
//			for(map<int,LBlock>::iterator mi = PLT.begin();
//					mi != PLT.end(); mi++){
//				setvalue((*mi).second._nzval,SCALAR_ZERO); 
//			}
//
//			// GEMM for PLT(i) += Sinv(i,j) * PL(j), L part
//			for( map<int,LBlock>::iterator mj = PL.begin(); mj != PL.end(); mj++ ){
//				int jsup     = (*mj).first;
//				LBlock& PLj  = (*mj).second;
//				int pjcol    = PCOL(jsup, grid);
//
//				if( mycol != pjcol || jsup <= ksup ) continue;  // LLIN: IMPORTANT
//
//				int jb = LBj(jsup, grid);
//
//				for(int ib = 0; ib < PMloc.nblkl(jb); ib++){
//					LBlock& LB = PMloc.L(jb)[ib];
//					int isup = LB.blkind();
//					map<int,LBlock>::iterator mi= PLT.find(isup);
//					if( mi == PLT.end() ) continue; // LLIN: IMPORTANT
//
//					LBlock& PLTi = (*mi).second;
//					// PMloc(isup, jsup) * PL(jsup)
//
//					// Find the subset of row and column indices
//					IntNumVec& LBrows = LB.rows();
//					IntNumVec rowidx(PLTi.nrows()), colidx(PLj.nrows());
//					for(int i = 0; i < PLTi.nrows(); i++){
//						for(int j = 0; j < LB.nrows(); j++){
//							if(PLTi.rows(i) == LBrows(j)){
//								rowidx[i] = j;
//								break;
//							}
//						}
//					}
//
//
//					int colsta = PMloc.xsup(jsup);
//					for(int i = 0; i < PLj.nrows(); i++){
//						colidx[i] = PLj.rows(i) - colsta;
//					}
//
//
//					// Setup the buffer matrix first
//					NumMat<Scalar> BufMat(PLTi.nrows(),PLj.nrows());
//					setvalue(BufMat, SCALAR_ZERO);
//					NumMat<Scalar>& LBMat = LB.nzval();
//
//					for(int j = 0; j < BufMat.n(); j++){
//						for(int i = 0; i < BufMat.m(); i++){
//							BufMat(i,j) = LBMat(rowidx[i],colidx[j]);
//						}
//					}
//
//					// GEMM 
//					{
//						int m = BufMat.m(), k = BufMat.n(), n = PLj.ncols();
//						// FIXME Use a better interface for blas...
//#ifdef _USE_COMPLEX_
//						Complex ZONE = Z_ONE; char notrans = 'N';
//						zgemm_(&notrans, &notrans, 
//								&m, &n, &k, (doublecomplex*)&ZONE, (doublecomplex*)BufMat.data(),
//								&m, (doublecomplex*)PLj.nzval().data(), &k, (doublecomplex*)&ZONE, 
//								(doublecomplex*)PLTi.nzval().data(), &m);  
//#else
//						double DONE = D_ONE; char notrans = 'N';
//						dgemm_(&notrans, &notrans, 
//								&m, &n, &k, &DONE, BufMat.data(),
//								&m, PLj.nzval().data(), &k, &DONE, 
//								PLTi.nzval().data(), &m);  
//#endif
//					}
//
//
//#if (__DEBUGlevel >= 2 )
//					cerr << "Rows: " << PMloc.xsup(isup) + 1 << " -- " <<
//						PMloc.xsup(isup+1) << endl; 
//					cerr << "Cols: " << PMloc.xsup(ksup) + 1 << " -- " <<
//						PMloc.xsup(ksup+1) << endl;
//					cerr << PLTi << endl;
//#endif
//
//				} // for (ib)
//			} // jsup
//
//			// GEMM for PLT(i) += Sinv(i,j) * PL(j), U part
//			for( map<int,LBlock>::iterator mi = PLT.begin(); mi != PLT.end(); mi++){
//				int isup  = (*mi).first;
//				int pirow = PROW(isup, grid);
//				if( myrow != pirow || isup <= ksup ) continue;
//				int ib = LBi(isup, grid);
//				LBlock& PLTi = (*mi).second;
//
//				for(int jb = 0; jb < PMloc.nblku(ib); jb++){
//					UBlock& UB = PMloc.U(ib)[jb];
//					int jsup = UB.blkind();
//					if( PL.count(jsup) == 0 ) continue;
//					map<int,LBlock>::iterator mj = PL.find(jsup);
//
//					LBlock& PLj  = (*mj).second;
//					// Mult(PMloc(isup,jsup),PL(isup),rowindices, column
//					// indices) to PLT; 
//
//					//	  cerr << "mpirank = " << mpirank << endl;
//					//	  cerr << (*mi).first << endl << (*mi).second << endl;
//
//					// Find the subset of row indices
//					IntNumVec& UBcols = UB.cols();
//					IntNumVec rowidx(PLTi.nrows()), colidx(PLj.nrows());
//					int rowsta = PMloc.xsup(isup);
//					for(int i = 0; i < PLTi.nrows(); i++){
//						rowidx[i] = PLTi.rows(i) - rowsta;
//					}
//
//					for(int i = 0; i < PLj.nrows(); i++){
//						for(int j = 0; j < UB.ncols(); j++){
//							if(PLj.rows(i) == UBcols(j)){
//								colidx[i] = j;
//								break;
//							}
//						}
//					}
//
//					// Setup the buffer matrix 
//					NumMat<Scalar> BufMat(PLTi.nrows(),PLj.nrows());
//					setvalue(BufMat, SCALAR_ZERO);
//					NumMat<Scalar>& UBMat = UB.nzval();
//					for(int j = 0; j < BufMat.n(); j++){
//						for(int i = 0; i < BufMat.m(); i++){
//							BufMat(i,j) = UBMat(rowidx[i],colidx[j]);
//						}
//					}
//
//					// GEMM 
//					{
//						int m = BufMat.m(), k = BufMat.n(), n = PLj.ncols();
//						// FIXME Use a better interface for blas...
//						
//#ifdef _USE_COMPLEX_
//						Complex ZONE = Z_ONE; char notrans = 'N';
//						zgemm_(&notrans, &notrans, 
//								&m, &n, &k, (doublecomplex*)&ZONE, (doublecomplex*)BufMat.data(),
//								&m, (doublecomplex*)PLj.nzval().data(), &k, (doublecomplex*)&ZONE, 
//								(doublecomplex*)PLTi.nzval().data(), &m);  
//#else
//						double DONE = D_ONE; char notrans = 'N';
//						dgemm_(&notrans, &notrans, 
//								&m, &n, &k, &DONE, BufMat.data(),
//								&m, PLj.nzval().data(), &k, &DONE, 
//								PLTi.nzval().data(), &m);  
//#endif
//					}
//#if (__DEBUGlevel >= 2 )
//					cerr << "Rows: " << PMloc.xsup(isup) + 1 << " -- " <<
//						PMloc.xsup(isup+1) << endl; 
//					cerr << "Cols: " << PMloc.xsup(ksup) + 1 << " -- " <<
//						PMloc.xsup(ksup+1) << endl;
//					cerr << PLTi << endl;
//#endif
//
//				} // for (ib)
//			} // jsup
//
//		}
//
//#if (__DEBUGlevel >= 2 )
//		if(mpirank == 0){
//			cerr << "Matvec passed" << endl;
//		}
//#endif
//
//#if (__PRNTlevel >= 3 )
//		if(grid->iam == 0){
//			cout << "SinvPL: matvec passed" << endl;
//		}
//#endif
//
//
//
//		// After the computation of each block, send it to
//		// (isup, ksup).
//		// (isup, ksup) receive from all processors, and perform local
//		// reduce. Can overlap communication with computation.
//		{
//			vector<int> mask(LBlock_Number,0);  
//			mask[LBlock_nzval] = 1; // Only need to communicate nzval
//
//			MPI_Request*    reqs   = new MPI_Request[npcolRestricted];   iA( reqs != NULL ); 
//			MPI_Status *    stats  = new MPI_Status [npcolRestricted];   iA( stats != NULL );
//
//			for(int inode = 0; inode < sizeksup; inode++){
//				int isup = nodesksup[inode];
//				int pirow = PROW(isup, grid);
//
//				if( isup > ksup && myrow == pirow ){
//					// Sender/Receiver only has at most one copy of the data
//					// Receiver overlaps communication with computation
//					stringstream   sstm;
//					int            sizeMsg;
//					vector<char>   strchr;
//
//					// Everyone does serialization, including the receiver in
//					// order to know the size of message. The size of messages
//					// should be the same for everyone
//					serialize(PLT[isup], sstm, mask);
//					sizeMsg = (sstm.str()).length();
//
//					int recvID = PNUM(pirow, pkcol, grid);
//
//					// Receiver ID put contribution back to PMloc
//					LBlock* ptrLB = NULL;
//					if( mpirank == recvID ){
//						for(int ib = 0; ib < PMloc.nblkl(kb); ib++){
//							if(PMloc.L(kb)[ib].blkind() == isup) {
//								ptrLB = &PMloc.L(kb)[ib];
//								// PMloc(isup,ksup) is set to be the local PLT(isup) and
//								// to be updated below 
//								{
//									Scalar *ptr0 = (*ptrLB)._nzval.data();
//									Scalar *ptr1 = PLT[isup]._nzval.data();
//									int sz = ptrLB->nrows() * ptrLB->ncols();
//									for(int i = 0; i < sz; i++){
//										*(ptr0++) = *(ptr1++);
//									}
//								}
//								break;
//							}
//						}
//						iA( ptrLB != NULL );
//					} // if( mpirank == recvID )
//
//
//#if (__DEBUGlevel >= 2 )
//					if(mpirank == 0){
//						cerr << "isup = " << isup << endl;
//						cerr << PLT[isup]._rows << endl;
//						cerr << PLT[isup]._nzval << endl;
//					}
//#endif
//
//
//
//
//
//					// Use synchronous Send/Recv to perform reduce operation to
//					// avoid a large buffer.  Rowcomm is necessary here to improve
//					// efficiency
//					for(int jp = 0; jp < npcolRestricted; jp++){
//						int pjcol = pcolRestricted[jp];
//						int sendID = PNUM(pirow, pjcol, grid);
//						if( sendID != recvID ){
//							if( mpirank == sendID ){
//								const string& sstr = sstm.str();
//								MPI_Send((void*)sstr.c_str(), sizeMsg, MPI_BYTE,
//										pkcol, pjcol, rowcomm);
//							}
//							if( mpirank == recvID ){
//								strchr.clear();
//								strchr.resize(sizeMsg);
//								MPI_Recv((void*)(&strchr[0]), sizeMsg, MPI_BYTE,
//										pjcol, pjcol, rowcomm, &(stats[jp]));
//
//								// Clear its own stringstream. IMPORTANT
//								sstm.str("");
//								sstm.write(&strchr[0], sizeMsg);
//								LBlock LBtmp;
//								deserialize(LBtmp, sstm, mask);
//								// Add the contribution to PMloc from other processors
//								{
//									Scalar *ptr0 = ptrLB->_nzval.data();
//									Scalar *ptr1 = LBtmp._nzval.data();
//									int sz = ptrLB->nrows() * ptrLB->ncols();
//									for(int i = 0; i < sz; i++){
//										*(ptr0++) += *(ptr1++);
//									}
//								}
//							}
//						} // if( sendID != recvID )
//					} // for( jp )
//				}
//			} // for( inode )
//			delete[] reqs;    reqs = NULL;
//			delete[] stats;   stats = NULL;
//		}
//
//
//#if (__DEBUGlevel >= 2 )
//		if(mpirank == 0){
//			cerr << "Reduce distribution passed" << endl;
//		}
//#endif
//
//#if (__PRNTlevel >= 3 )
//		if(grid->iam == 0){
//			cout << "SinvPL: Reduce distribution passed" << endl;
//		}
//#endif
//
//	}
//	return 0;
//}
//
//
///**********************************************************************
// * Update the diagonal block by
// *   U_kk^-1 U_ik S^-1_ij L_jk L^-1_kk
// *
// * Blocks owning L(i,k) sending L(i,k) to processors owning block (k,i)
// *
// * Processors owning block (k,i) perform inner product 
// *
// * All processors in the row group of k reduce the sum of the inner
// * products to the processor owning the diagonal block (k,k)
// **********************************************************************/
//int DiagInnerProd(PMatrix& PMloc, vector<vector<int> >& localEtree, int ksup, gridinfo_t* grid){
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int nprow = grid->nprow,   npcol = grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	int nsupers = PMloc.nsupers();
//	int NULLID = -1; // EMPTY means the block is not owned by any processor
//	const int NDATA = 4; 
//	enum {SUPID, SNDID, RCVID, SIZE};
//
//	iA(ksup >= 0 && ksup < nsupers-1);
//
//	MPI_Comm comm = grid->comm, colcomm = grid->cscp.comm, rowcomm = grid->rscp.comm;
//
//	int pkrow = PROW(ksup, grid), pkcol = PCOL(ksup, grid);
//
//	if(1){
//		map<int, LBlock> PL;
//
//		vector<int>& nodesksup = localEtree[ksup];
//		int sizeksup = nodesksup.size();
//
//		// Only processors with nonvanishing setksup participates SinvPL.
//		if(sizeksup == 0) return 0;
//		if( myrow != pkrow && mycol != pkcol ) return 0;
//
//		// Each processor knows the row and column processor id arrays
//		// related to supernode ksup
//		vector<int> prowRestricted;  
//		vector<int> pcolRestricted;  
//		{
//			set<int> prowtmp;
//			set<int> pcoltmp;
//
//			int isup, pirow, picol;
//			for(int i = 0; i < sizeksup; i++){
//				isup  = nodesksup[i];
//				pirow = PROW(isup, grid);
//				picol = PCOL(isup, grid);
//				prowtmp.insert(pirow);
//				pcoltmp.insert(picol);
//			}
//
//			prowRestricted.clear();
//			prowRestricted.insert(prowRestricted.end(), prowtmp.begin(),
//					prowtmp.end());
//			pcolRestricted.clear();
//			pcolRestricted.insert(pcolRestricted.end(), pcoltmp.begin(),
//					pcoltmp.end());
//		}
//		int nprowRestricted = prowRestricted.size();
//		int npcolRestricted = pcolRestricted.size();
//
//
//#if (__PRNTlevel >= 3)
//		if(grid->iam == 1){
//			cerr << "DiagInnerProd: initiated" << endl;
//		}
//#endif
//
//
//		vector<int> mask(LBlock_Number,1);  
//
//
//		// Communication based on nodesksup. If inode does not need to be
//		// communicated, everything is set to be NULL
//		{
//			MPI_Request*    reqs   = new MPI_Request[sizeksup];   iA( reqs != NULL ); 
//			MPI_Status *    stats  = new MPI_Status [sizeksup];   iA( stats != NULL );
//			vector<stringstream*>  sstm(sizeksup);
//			vector<int>            sizeMsg(sizeksup);
//			vector<vector<char> >  strchr(sizeksup);          
//
//			for(int inode = 0; inode < sizeksup; inode++){
//				reqs[inode]  = MPI_REQUEST_NULL;
//				sstm[inode] = NULL;
//				sizeMsg[inode] = 0;
//				strchr[inode].clear();
//			}
//
//			// Sender
//			if( mycol == pkcol ){
//				int kb = LBj(ksup, grid);
//
//				for(int inode = 0; inode < sizeksup; inode++){
//					int isup = nodesksup[inode];
//					int pirow = PROW(isup, grid);
//					int picol = PCOL(isup, grid);
//					if( isup > ksup ){
//
//						// Pack the data to be sent, and send size information
//						int sendID = PNUM(pirow, pkcol, grid);
//						int recvID = PNUM(pkrow, picol, grid);
//
//						// Do not send to itself
//						if( mpirank == sendID ){
//							LBlock *ptrLB = NULL;
//							for(int ib = 0; ib < PMloc.nblkl(kb); ib++){
//								if( PMloc.L(kb)[ib].blkind() == isup ){
//									ptrLB = &PMloc.L(kb)[ib];
//								}
//							}
//							iA( ptrLB != NULL );
//
//							if( sendID == recvID ){
//								PL[isup] = *ptrLB;
//							}
//							else{
//								sstm[inode] = new stringstream();
//								serialize(*ptrLB, *(sstm[inode]), mask);
//								sizeMsg[inode] = (sstm[inode]->str()).length();
//								MPI_Isend(&sizeMsg[inode], 1, MPI_INT, recvID, inode, comm, &(reqs[inode]));
//							}
//						}
//
//					}
//				} // for (inode)
//				MPI_Waitall(sizeksup, &(reqs[0]), &(stats[0]));
//
//				for(int inode = 0; inode < sizeksup; inode++){
//					reqs[inode]  = MPI_REQUEST_NULL;
//				}
//
//				for(int inode = 0; inode < sizeksup; inode++){
//					int isup = nodesksup[inode];
//					int pirow = PROW(isup, grid);
//					int picol = PCOL(isup, grid);
//					if( isup > ksup ){
//						int sendID = PNUM(pirow, pkcol, grid);
//						int recvID = PNUM(pkrow, picol, grid);
//						if( sendID != recvID && mpirank == sendID){
//							const string& sstr = sstm[inode]->str();
//							MPI_Isend((void*)sstr.c_str(), sizeMsg[inode], MPI_BYTE,
//									recvID, inode, comm, &(reqs[inode]));
//						}
//					}
//				} // for(inode)
//				MPI_Waitall(sizeksup, &(reqs[0]), &(stats[0]));
//			} // if (mycol == pkcol)
//
//			//Receiver
//			if( myrow == pkrow){
//				for(int inode = 0; inode < sizeksup; inode++){
//					int isup = nodesksup[inode];
//					int pirow = PROW(isup, grid);
//					int picol = PCOL(isup, grid);
//					if( isup > ksup ){
//						int sendID = PNUM(pirow, pkcol, grid);
//						int recvID = PNUM(pkrow, picol, grid);
//						if( sendID != recvID && mpirank == recvID ){
//							MPI_Irecv(&sizeMsg[inode], 1, MPI_INT, sendID, 
//									inode, comm, &(reqs[inode]));
//						}
//					}
//				} // for(inode)
//				MPI_Waitall(sizeksup, &(reqs[0]), &(stats[0]));
//
//				for(int inode = 0; inode < sizeksup; inode++){
//					reqs[inode]  = MPI_REQUEST_NULL;
//				}
//
//				for(int inode = 0; inode < sizeksup; inode++){
//					int isup = nodesksup[inode];
//					int pirow = PROW(isup, grid);
//					int picol = PCOL(isup, grid);
//					if( isup > ksup ){
//						int sendID = PNUM(pirow, pkcol, grid);
//						int recvID = PNUM(pkrow, picol, grid);
//						if( sendID != recvID && mpirank == recvID ){
//							strchr[inode].clear();
//							strchr[inode].resize(sizeMsg[inode]);
//							MPI_Irecv((void*)(&(strchr[inode][0])), sizeMsg[inode], MPI_BYTE, 
//									sendID, inode, comm, &(reqs[inode]));
//						}
//					}
//				} // for(inode)
//				MPI_Waitall(sizeksup, &(reqs[0]), &(stats[0]));
//
//				// Deserialize into PL
//				for(int inode = 0; inode < sizeksup; inode++){
//					int isup = nodesksup[inode];
//					if(strchr[inode].size() > 0 ){
//						sstm[inode] = new stringstream();
//						sstm[inode]->write(&(strchr[inode][0]), sizeMsg[inode]);
//						iA(PL.count(isup) == 0); // PL should be empty
//						deserialize(PL[isup], *(sstm[inode]), mask);
//					}
//				}
//			}
//
//			for(int inode = 0; inode < sizeksup; inode++){
//				strchr[inode].clear();
//				if( sstm[inode] != NULL )   delete sstm[inode];
//				delete[] reqs;    reqs = NULL;
//				delete[] stats;   stats = NULL;
//			}	
//		}
//
//#if (__PRNTlevel >= 3)
//		if(grid->iam == 1){
//			cerr << "DiagInnerProd: Communication of PL passed" << endl;
//		}
//#endif
//
//		// Inner product: NOTE: The rows and columns of the diagonal blocks 
//		// are ordered
//
//		if( myrow == pkrow ){
//			NumMat<Scalar> tDiagMat;
//			tDiagMat.resize(PMloc.supsize(ksup), PMloc.supsize(ksup)); 
//			setvalue(tDiagMat, SCALAR_ZERO);
//
//			int kb = LBi(ksup, grid);
//
//
//			for(map<int, LBlock>::iterator mi = PL.begin(); 
//					mi != PL.end(); mi++){
//				LBlock& LB = (*mi).second;
//				int jsup = LB.blkind();
//				for(int jb = 0; jb < PMloc.nblku(kb); jb++){
//					UBlock& UB = PMloc.U(kb)[jb];
//					int tjsup = UB.blkind();
//					if( jsup == tjsup ){
//						// GEMM
//						int tm = UB.nrows(), tn = LB.ncols(), tk = LB.nrows();
//						iA(tm == PMloc.supsize(ksup) && 
//								tn == PMloc.supsize(ksup) );
//						// Match the indices of U
//						// Set the column of UBMat to be 0 if no columns is matching
//						// the rows of LB
//						NumMat<Scalar> UBMat(tm, tk);  
//						setvalue(UBMat, SCALAR_ZERO);
//						for(int ti = 0; ti < LB.nrows(); ti++){
//							for(int tj = 0; tj < UB.ncols(); tj++){
//								if( LB.rows(ti) == UB.cols(tj) ){
//									int ONE = 1;
//#ifdef _USE_COMPLEX_
//									zlacpy_("A", &tm, &ONE, UB.nzval().clmdata(tj), &tm, 
//											UBMat.clmdata(ti), &tm);
//#else
//									dlacpy_("A", &tm, &ONE, UB.nzval().clmdata(tj), &tm, 
//											UBMat.clmdata(ti), &tm);
//#endif
//									break;
//								}
//							}
//						}
//
//						// GEMM
//#ifdef _USE_COMPLEX_
//						char notrans = 'N'; Complex ZONE = Z_ONE;
//						zgemm_(&notrans, &notrans, &tm, &tn, &tk, 
//								(doublecomplex*)&ZONE, (doublecomplex*)UBMat.data(), &tm,
//								(doublecomplex*)LB.nzval().data(), &tk, 
//								(doublecomplex*)&ZONE, (doublecomplex*)tDiagMat.data(),
//								&tm);  
//#else
//						char notrans = 'N'; double DONE = D_ONE;
//						dgemm_(&notrans, &notrans, &tm, &tn, &tk, &DONE, UBMat.data(), &tm,
//								LB.nzval().data(), &tk, &DONE, tDiagMat.data(),
//								&tm);  
//#endif
//					}
//
//				} // for(jb)
//			} // for (mi)
//
//
//
//			// Reduce tDiagMat to the diagonal block
//			{
//				// Update the diagonal block
//				LBlock *ptrLB;
//				MPI_Status stat;
//
//				if(myrow == pkrow && mycol == pkcol){
//					int kb = LBj(ksup, grid);
//					ptrLB = &(PMloc.L(kb)[0]); // Diagonal block
//					// IMPORTANT: DO NOT clear this block to 0
//					{
//						Scalar *ptr0 = ptrLB->_nzval.data();
//						Scalar *ptr1 = tDiagMat.data();
//						int sz = ptrLB->nrows() * ptrLB->ncols();
//						for(int i = 0; i < sz; i++){
//							*(ptr0++) += *(ptr1++);
//						}
//					}
//				}
//
//				for(int jp = 0; jp < npcolRestricted; jp++){
//					int pjcol = pcolRestricted[jp];
//					int sendID = PNUM(pkrow, pjcol, grid);
//					int recvID = PNUM(pkrow, pkcol, grid);
//					if( sendID != recvID ){
//						if(mpirank == sendID ){
//							int sz = tDiagMat.m() * tDiagMat.n();
//							MPI_Send((void*)tDiagMat.data(), 
//									sz*sizeof(Scalar), MPI_BYTE,
//									recvID, 1, comm);
////							MPI_Send(tDiagMat.data(), sz, MPI_DOUBLE,
////									recvID, 1, comm);
//						}
//						if(mpirank == recvID){
//							// Set tDiagMat to be 0 to receive from other processors.
//							setvalue(tDiagMat, SCALAR_ZERO);
//							int sz = tDiagMat.m() * tDiagMat.n();
//							MPI_Recv((void*)tDiagMat.data(), 
//										sz*sizeof(Scalar), MPI_BYTE,
//										sendID, 1, comm, &stat);
////							MPI_Recv(tDiagMat.data(), sz, MPI_DOUBLE,
////									sendID, 1, comm, &stat);
//							{
//								Scalar *ptr0 = ptrLB->_nzval.data();
//								Scalar *ptr1 = tDiagMat.data();
//								int sz = ptrLB->nrows() * ptrLB->ncols();
//								for(int i = 0; i < sz; i++){
//									*(ptr0++) += *(ptr1++);
//								}
//							}
//						}
//					} // if (sendID != recvID)
//				} // for(jp) 
//			}
//
//
//#if (__PRNTlevel >= 3)
//			if(grid->iam == 1){
//				cerr << "DiagInnerProd: Diagonal updated" << endl;
//			}
//#endif
//
//
//		} // if( myrow == pkrow )
//
//	}
//	return 0;
//}
//
//
///**********************************************************************
// * PUT(i) = PU(i) S^-1(i,j)
// *
// * Build PUT in each processor from
// * (ksup+1:nsupers,ksup+1:nsupers)
// *
// * Build PU(i) = -U^-1(k,k) U(i,k) in each processor from
// * (ksup+1:nsupers,ksup+1:nsupers)
// *
// * Local matrix-matrix multiplication
// *
// * Reduce the PUT to processors (ksup, ksup+1:nsupers)
// **********************************************************************/
//int SinvPU(PMatrix& PMloc, vector<vector<int> >& localEtree, int ksup, gridinfo_t* grid){
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	int nprow = grid->nprow,   npcol = grid->npcol;
//	int nsupers = PMloc.nsupers();
//	const int NULLID = -1;
//	const int NDATA = 1;
//	enum {SUPID};
//
//	double t0, t1;
//
//	map<int, UBlock> PU;
//	map<int, UBlock> PUT;
//
//	PU.clear();
//	PUT.clear();
//
//	iA(ksup >= 0 && ksup < nsupers-1);
//	int pkrow = PROW(ksup, grid), pkcol = PCOL(ksup, grid);
//	int kb = LBi(ksup, grid);
//
//	MPI_Comm comm = grid->comm;
//	MPI_Comm colcomm = grid->cscp.comm, rowcomm = grid->rscp.comm;
//
//	// Use localEtree for communication
//	if(1)
//	{
//		vector<int>& nodesksup = localEtree[ksup];
//		int sizeksup = nodesksup.size();
//
//		// Only processors with nonvanishing setksup participates SinvPU.
//		if(sizeksup == 0) return 0;
//
//
//		// Each processor knows the row and column processor id arrays that
//		// is restricted SinvPL.
//		vector<int> prowRestricted;  
//		vector<int> pcolRestricted;  
//		{
//			set<int> prowtmp;
//			set<int> pcoltmp;
//
//			int isup, pirow, picol;
//			for(int i = 0; i < sizeksup; i++){
//				isup  = nodesksup[i];
//				pirow = PROW(isup, grid);
//				picol = PCOL(isup, grid);
//				prowtmp.insert(pirow);
//				pcoltmp.insert(picol);
//			}
//
//			prowRestricted.clear();
//			prowRestricted.insert(prowRestricted.end(), prowtmp.begin(),
//					prowtmp.end());
//			pcolRestricted.clear();
//			pcolRestricted.insert(pcolRestricted.end(), pcoltmp.begin(),
//					pcoltmp.end());
//		}
//		int nprowRestricted = prowRestricted.size();
//		int npcolRestricted = pcolRestricted.size();
//
//#if (__PRNTlevel >= 3 )
//		if(grid->iam == 0){
//			cout << "SinvPU: initiation passed" << endl;
//		}
//#endif
//
//		// Setup PLT, in preparation for PLT = PM * PL
//		{
//			MPI_Request*    reqs   = new MPI_Request[nprowRestricted];   iA( reqs != NULL ); 
//			MPI_Status *    stats  = new MPI_Status [nprowRestricted];   iA( stats != NULL );
//			vector<int> mask(UBlock_Number,1); // Communication for all
//
//			for(int jnode = 0; jnode < sizeksup; jnode++){
//				int jsup  = nodesksup[jnode];
//				int pjcol = PCOL(jsup, grid);
//				if( jsup > ksup && mycol == pjcol ){
//
//					// Sender/Receiver only has at most one copy of the data
//					stringstream   sstm;
//					int            sizeMsg;
//					vector<char>   strchr;
//
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						reqs[ip]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					// Sender
//					if( myrow == pkrow ){
//						UBlock *ptrUB = NULL;
//						for(int jb = 0; jb < PMloc.nblku(kb); jb++){
//							if( PMloc.U(kb)[jb].blkind() == jsup ){
//								ptrUB = &PMloc.U(kb)[jb];
//							}
//						}
//						iA( ptrUB != NULL );
//						PUT[jsup] = *ptrUB; // setup local PUT
//						serialize(*ptrUB, sstm, mask);
//						sizeMsg = (sstm.str()).length();
//
//					} 
//
//					// Size information
//					int sendID = PNUM(pkrow, pjcol, grid);
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						int pirow  = prowRestricted[ip];
//						int recvID = PNUM(pirow, pjcol, grid);
//						if( sendID != recvID ){
//							if( mpirank == sendID ){
//								MPI_Isend(&sizeMsg, 1, MPI_INT, pirow, pirow, colcomm, &(reqs[ip]));
//							}
//							if( mpirank == recvID ){
//								MPI_Irecv(&sizeMsg, 1, MPI_INT, pkrow, pirow, colcomm, &(reqs[ip])); 
//							}
//						}
//					}
//					MPI_Waitall(nprowRestricted, &(reqs[0]), &(stats[0]));
//
//
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						reqs[ip]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						int pirow  = prowRestricted[ip];
//						int recvID = PNUM(pirow, pjcol, grid);
//						if( sendID != recvID ){
//							if( mpirank == sendID ){
//								const string& sstr = sstm.str();
//								MPI_Isend((void*)sstr.c_str(), sizeMsg, MPI_BYTE, 
//										pirow, pirow, colcomm, &(reqs[ip]));
//							}
//							if( mpirank == recvID ){
//								strchr.clear();
//								strchr.resize(sizeMsg);
//								MPI_Irecv((void*)(&strchr[0]), sizeMsg, MPI_BYTE, 
//										pkrow, pirow, colcomm, &(reqs[ip])); 
//							}
//						}
//					}
//					MPI_Waitall(nprowRestricted, &(reqs[0]), &(stats[0]));
//
//					// Deserialize into PUT
//					if( mpirank != sendID ){
//						sstm.write(&strchr[0], sizeMsg);
//						iA(PUT.count(jsup) == 0);  // PLT should be empty
//						deserialize(PUT[jsup], sstm, mask);
//					}
//				} 
//
//			}// for (jnode)
//			delete[] reqs;    reqs = NULL;
//			delete[] stats;   stats = NULL;
//		}
//
//#if (__PRNTlevel >= 3 )
//		if(grid->iam == 0){
//			cout << "SinvPU: PUT distribution passed" << endl;
//		}
//#endif
//
//		// Setup PU
//		{
//			MPI_Request*    reqs   = new MPI_Request[npcolRestricted];   iA( reqs  != NULL ); 
//			MPI_Status *    stats  = new MPI_Status [npcolRestricted];   iA( stats != NULL );
//			vector<int> mask(UBlock_Number,1);  // all data is needed for PL.
//
//			for(int jnode = 0; jnode < sizeksup; jnode++){
//				int jsup  = nodesksup[jnode];
//				int pirow = PROW(jsup, grid);
//				int pjcol = PCOL(jsup, grid);
//				// Sender
//				// The diagonal processor sends PU to column processors
//				int sendID = PNUM(pirow, pjcol, grid); 
//				if( jsup > ksup && mpirank == sendID){
//
//					stringstream   sstm;
//					int            sizeMsg;
//
//					for(int ip = 0; ip < npcolRestricted; ip++){
//						reqs[ip]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					PU[jsup] = PUT[jsup]; // Should have received from PM(ksup,jsup)
//					serialize(PU[jsup], sstm, mask);
//					sizeMsg = (sstm.str()).length();
//
//					// Size information
//					for(int jp = 0; jp < npcolRestricted; jp++){
//						int pjcolRecv  = pcolRestricted[jp];
//						int recvID = PNUM(pirow, pjcolRecv, grid);
//						if( sendID != recvID ){
//							MPI_Isend(&sizeMsg, 1, MPI_INT, pjcolRecv, pjcolRecv, rowcomm, &(reqs[jp]));
//						}
//					}
//					MPI_Waitall(npcolRestricted, &(reqs[0]), &(stats[0]));
//
//					// Actual data
//					for(int jp = 0; jp < npcolRestricted; jp++){
//						reqs[jp]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					for(int jp = 0; jp < npcolRestricted; jp++){
//						int pjcolRecv  = pcolRestricted[jp];
//						int recvID = PNUM(pirow, pjcolRecv, grid);
//						if( sendID != recvID ){
//							const string& sstr = sstm.str();
//							MPI_Isend((void*)sstr.c_str(), sizeMsg, MPI_BYTE, 
//									pjcolRecv, pjcolRecv, rowcomm, &(reqs[jp]));
//						}
//					}
//					MPI_Waitall(npcolRestricted, &(reqs[0]), &(stats[0]));
//				}  // if( mpirank == sendID )
//
//				// Receiver
//				if( jsup > ksup && myrow == pirow ){
//					stringstream   sstm;
//					int            sizeMsg;
//					vector<char>   strchr;
//
//					// Size information
//					for(int jp = 0; jp < npcolRestricted; jp++){
//						reqs[jp]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					for(int jp = 0; jp < npcolRestricted; jp++){
//						int pjcolRecv  = pcolRestricted[jp];
//						int recvID = PNUM(pirow, pjcolRecv, grid);
//						if( sendID != recvID && mpirank == recvID ){
//							MPI_Irecv(&sizeMsg, 1, MPI_INT, pjcol, pjcolRecv, rowcomm, &(reqs[jp]));
//						}
//					}
//					MPI_Waitall(npcolRestricted, &(reqs[0]), &(stats[0]));
//
//					// Actual data
//					for(int jp = 0; jp < npcolRestricted; jp++){
//						reqs[jp]   = MPI_REQUEST_NULL;  // LLIN: IMPORTANT
//					}
//
//					for(int jp = 0; jp < npcolRestricted; jp++){
//						int pjcolRecv  = pcolRestricted[jp];
//						int recvID     = PNUM(pirow, pjcolRecv, grid);
//						if( sendID != recvID && mpirank == recvID ){
//							strchr.clear();
//							strchr.resize(sizeMsg);
//							MPI_Irecv((void*)(&strchr[0]), sizeMsg, MPI_BYTE,
//									pjcol, pjcolRecv, rowcomm, &(reqs[jp]));
//						}
//					}
//					MPI_Waitall(npcolRestricted, &(reqs[0]), &(stats[0]));
//
//					// Deserialize into PL
//					if( mpirank != sendID ){
//						sstm.write(&strchr[0], sizeMsg);
//						iA(PU.count(jsup) == 0);  // PU should be empty
//						deserialize(PU[jsup], sstm, mask);
//					} 
//
//				}
//			}// for (jnode)
//			delete[] reqs;    reqs = NULL;
//			delete[] stats;   stats = NULL;
//		}
//
//#if (__PRNTlevel >= 3 )
//		if(grid->iam == 0){
//			cout << "SinvPU: PU distribution passed" << endl;
//		}
//#endif
//
//		{
//			// For each processor, clean the PUT.  LLIN: IMPORTANT
//			for(map<int,UBlock>::iterator mi = PUT.begin();
//					mi != PUT.end(); mi++){
//				setvalue((*mi).second._nzval, SCALAR_ZERO); 
//			}
//
//			// GEMM for PUT(j) += PU(i) * Sinv(i,j), L part
//			for( map<int,UBlock>::iterator mj = PUT.begin(); mj != PUT.end();
//					mj++ ){
//				int jsup     = (*mj).first;
//				UBlock& PUTj = (*mj).second;
//				int pjcol    = PCOL(jsup, grid);
//
//
//				if( mycol != pjcol || jsup <= ksup ) continue;  // LLIN: IMPORTANT
//				int jb = LBj(jsup, grid);
//
//				for(int ib = 0; ib < PMloc.nblkl(jb); ib++){
//					LBlock& LB = PMloc.L(jb)[ib];
//					int isup = LB.blkind();
//					if( PU.count(isup) == 0 ) continue; // LLIN: IMPORTANT
//					map<int,UBlock>::iterator mi = PU.find(isup);
//
//					UBlock& PUi  = (*mi).second;
//					// PU(isup) * PMloc(isup, jsup)
//
//					// Find the subset of row and column indices
//					IntNumVec& LBrows = LB.rows();
//					IntNumVec rowidx(PUi.ncols()), colidx(PUTj.ncols());
//					setvalue(rowidx, 0);  setvalue(colidx, 0);
//					for(int i = 0; i < PUi.ncols(); i++){
//						for(int j = 0; j < LB.nrows(); j++){
//							if(PUi.cols(i) == LBrows(j)){
//								rowidx[i] = j;
//								break;
//							}
//						}
//					}
//
//					int colsta = PMloc.xsup(jsup);
//					for(int i = 0; i < PUTj.ncols(); i++){
//						colidx[i] = PUTj.cols(i) - colsta;
//					}
//
//
//					// Setup the buffer matrix first
//					NumMat<Scalar> BufMat(PUi.ncols(),PUTj.ncols());
//					setvalue(BufMat, SCALAR_ZERO);
//					NumMat<Scalar>& LBMat = LB.nzval();
//					for(int j = 0; j < BufMat.n(); j++){
//						for(int i = 0; i < BufMat.m(); i++){
//							BufMat(i,j) = LBMat(rowidx[i],colidx[j]);
//						}
//					}
//
//					// GEMM 
//					{
//						int m = PUi.nrows(), n = PUTj.ncols(), k = PUi.ncols();
//#ifdef _USE_COMPLEX_
//						Complex ZONE = Z_ONE; char notrans='N';
//						zgemm_(&notrans, &notrans, &m, &n, &k, (doublecomplex*)&ZONE, 
//								(doublecomplex*)PUi.nzval().data(), &m, 
//								(doublecomplex*)BufMat.data(), &k, 
//								(doublecomplex*)&ZONE, (doublecomplex*)PUTj.nzval().data(), &m);
//#else
//						double DONE = D_ONE; char notrans='N';
//						dgemm_(&notrans, &notrans, &m, &n, &k, &DONE, 
//								PUi.nzval().data(), &m, 
//								BufMat.data(), &k, 
//								&DONE, PUTj.nzval().data(), &m);
//#endif
//					}
//#if (__DEBUGlevel >= 2 )
//					cerr << "Rows: " << PMloc.xsup(ksup) + 1 << " -- " <<
//						PMloc.xsup(ksup+1) << endl; 
//					cerr << "Cols: " << PMloc.xsup(jsup) + 1 << " -- " <<
//						PMloc.xsup(jsup+1) << endl;
//					cerr << PUTj << endl;
//#endif
//
//				} // for (ib)
//			} // jsup
//
//			// GEMM for PUT(j) += PU(i) * Sinv(i,j), U part
//			for( map<int,UBlock>::iterator mi = PU.begin(); mi != PU.end();
//					mi++ ){
//				int isup     = (*mi).first;
//				UBlock& PUi  = (*mi).second;
//				int pirow    = PROW(isup, grid);
//
//				if( myrow != pirow || isup <= ksup ) continue;
//				int ib = LBi(isup, grid);
//
//				for(int jb = 0; jb < PMloc.nblku(ib); jb++){
//					UBlock& UB = PMloc.U(ib)[jb];
//					int jsup = UB.blkind();
//					if( PUT.count(jsup) == 0 ) continue; // LLIN: IMPORTANT
//					map<int,UBlock>::iterator mj = PUT.find(jsup);
//
//					UBlock& PUTj = (*mj).second;
//					// PU(isup) * PMloc(isup, jsup)
//
//					// Find the subset of row indices
//					IntNumVec& UBcols = UB.cols();
//					IntNumVec rowidx(PUi.ncols()), colidx(PUTj.ncols());
//					setvalue(rowidx, 0);  setvalue(colidx, 0);
//
//					for(int i = 0; i < PUTj.ncols(); i++){
//						for(int j = 0; j < UB.ncols(); j++){
//							if(PUTj.cols(i) == UBcols(j)){
//								colidx[i] = j;
//								break;
//							}
//						}
//					}
//
//					int rowsta = PMloc.xsup(isup);
//					for(int i = 0; i < PUi.ncols(); i++){
//						rowidx[i] = PUi.cols(i) - rowsta;
//					}
//
//					// Setup the buffer matrix 
//					NumMat<Scalar> BufMat(PUi.ncols(),PUTj.ncols());
//					setvalue(BufMat, SCALAR_ZERO);
//					NumMat<Scalar>& UBMat = UB.nzval();
//					for(int j = 0; j < BufMat.n(); j++){
//						for(int i = 0; i < BufMat.m(); i++){
//							BufMat(i,j) = UBMat(rowidx[i],colidx[j]);
//						}
//					}
//
//					// GEMM 
//					{
//						int m = PUi.nrows(), n = PUTj.ncols(), k = PUi.ncols();
//#ifdef _USE_COMPLEX_
//						Complex ZONE = Z_ONE; char notrans = 'N';
//						zgemm_(&notrans, &notrans, &m, &n, &k, (doublecomplex*)&ZONE, 
//								(doublecomplex*)PUi.nzval().data(), &m, 
//								(doublecomplex*)BufMat.data(), &k, 
//								(doublecomplex*)&ZONE, (doublecomplex*)PUTj.nzval().data(), &m);
//#else
//						double DONE = D_ONE; char notrans = 'N';
//						dgemm_(&notrans, &notrans, &m, &n, &k, &DONE, 
//								PUi.nzval().data(), &m, 
//								BufMat.data(), &k, 
//								&DONE, PUTj.nzval().data(), &m);
//#endif
//					}
//
//#if (__DEBUGlevel >= 2 )
//					cerr << "Rows: " << PMloc.xsup(ksup) + 1 << " -- " <<
//						PMloc.xsup(ksup+1) << endl; 
//					cerr << "Cols: " << PMloc.xsup(jsup) + 1 << " -- " <<
//						PMloc.xsup(jsup+1) << endl;
//					cerr << PUTj << endl;
//#endif
//
//				} // for (jb)
//			} // 
//		}
//
//#if (__PRNTlevel >= 3 )
//		if(grid->iam == 0){
//			cout << "SinvPU: matvec passed" << endl;
//		}
//#endif
//
//
//
//		// After the computation of each block, send it to
//		// (isup, ksup).
//		// (isup, ksup) receive from all processors, and perform local
//		// reduce. Can overlap communication with computation.
//		{
//			vector<int> mask(UBlock_Number,0);  
//			mask[UBlock_nzval] = 1; // Only need to communicate nzval
//
//			MPI_Request*    reqs   = new MPI_Request[nprowRestricted];   iA( reqs != NULL ); 
//			MPI_Status *    stats  = new MPI_Status [nprowRestricted];   iA( stats != NULL );
//
//			for(int jnode = 0; jnode < sizeksup; jnode++){
//				int jsup = nodesksup[jnode];
//				int pjcol = PCOL(jsup, grid);
//
//				if( jsup > ksup && mycol == pjcol ){
//					// Sender/Receiver only has at most one copy of the data
//					// Receiver overlaps communication with computation
//					stringstream   sstm;
//					int            sizeMsg;
//					vector<char>   strchr;
//
//					// Everyone does serialization, including the receiver in
//					// order to know the size of message. The size of messages
//					// should be the same for everyone
//					serialize(PUT[jsup], sstm, mask);
//					sizeMsg = (sstm.str()).length();
//
//					int recvID = PNUM(pkrow, pjcol, grid);
//
//					// Receiver ID put contribution back to PMloc
//					UBlock* ptrUB = NULL;
//					if( mpirank == recvID ){
//						for(int jb = 0; jb < PMloc.nblku(kb); jb++){
//							if(PMloc.U(kb)[jb].blkind() == jsup) {
//								ptrUB = &PMloc.U(kb)[jb];
//								// PMloc(ksup,jsup) is set to be the local PUT(jsup) and
//								// to be updated below 
//								{
//									Scalar *ptr0 = (*ptrUB)._nzval.data();
//									Scalar *ptr1 = PUT[jsup]._nzval.data();
//									int sz = ptrUB->nrows() * ptrUB->ncols();
//									for(int i = 0; i < sz; i++){
//										*(ptr0++) = *(ptr1++);
//									}
//								}
//								break;
//							}
//						}
//						iA( ptrUB != NULL );
//					} // if( mpirank == recvID )
//
//					// Use synchronous Send/Recv to perform reduce operation to
//					// avoid a large buffer.  Rowcomm is necessary here to improve
//					// efficiency
//					for(int ip = 0; ip < nprowRestricted; ip++){
//						int pirow = prowRestricted[ip];
//						int sendID = PNUM(pirow, pjcol, grid);
//						if( sendID != recvID ){
//							if( mpirank == sendID ){
//								const string& sstr = sstm.str();
//								MPI_Send((void*)sstr.c_str(), sizeMsg, MPI_BYTE,
//										pkrow, pirow, colcomm);
//							}
//							if( mpirank == recvID ){
//								strchr.clear();
//								strchr.resize(sizeMsg);
//								MPI_Recv((void*)(&strchr[0]), sizeMsg, MPI_BYTE,
//										pirow, pirow, colcomm, &(stats[ip]));
//
//								// Clear its own stringstream. IMPORTANT
//								sstm.str("");
//								sstm.write(&strchr[0], sizeMsg);
//								UBlock UBtmp;
//								deserialize(UBtmp, sstm, mask);
//								// Add the contribution to PMloc from other processors
//								{
//									Scalar *ptr0 = ptrUB->_nzval.data();
//									Scalar *ptr1 = UBtmp._nzval.data();
//									int sz = ptrUB->nrows() * ptrUB->ncols();
//									for(int i = 0; i < sz; i++){
//										*(ptr0++) += *(ptr1++);
//									}
//								}
//							}
//						} // if( sendID != recvID )
//					} // for( ip )
//				}
//			} // for( jnode )
//			delete[] reqs;    reqs = NULL;
//			delete[] stats;   stats = NULL;
//		}
//
//#if (__PRNTlevel >= 3 )
//		if(grid->iam == 0){
//			cout << "SinvPU: Reduce distribution passed" << endl;
//		}
//#endif
//
//	}
//	return 0;
//}
//
//
//
///**********************************************************************
// * Output the localEtree
// **********************************************************************/
//int DumpLocalEtree(vector<vector<int> >& localEtree, gridinfo_t* grid){
//	int mpirank;  mpirank = grid->iam; 
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int nsupers = localEtree.size();
//
//	// FIXME
//	if(mpirank == 0){
//		for(int ksup = 0; ksup < nsupers; ksup++){
//			for(int i = 0; i < localEtree[ksup].size(); i++){
//				fprintf(stderr, "(%5d,%5d)  ", localEtree[ksup][i], ksup);
//			}
//			fprintf(stderr, "\n");
//		} 
//	}
//	return 0;
//}
//
//
///**********************************************************************
// * Estimate the condition number for the diagonal blocks
// **********************************************************************/
////int PMatrix::CondDiagBlock(gridinfo_t* grid){
////	MPI_Comm comm = grid->comm;
////
////	int mpirank;  mpirank = grid->iam;
////	int mpisize;  mpisize = grid->nprow * grid->npcol;
////	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
////	iC( MPI_Barrier(comm) );
////
////	int globalIndexMaxCond, localIndexMaxCond;
////	double globalValueMaxCond, localValueMaxCond;
////
////	localIndexMaxCond = -1;
////	localValueMaxCond = -1.0;
////
////	for(int ib = 0; ib < this->nbc(); ib++){
////		int bnum = ( (ib) * grid->npcol ) + mycol;
////		if( bnum >= this->nsupers() ) continue;
////		if( this->nblkl(ib) ){ // Not an empty column
////			for(int iblk = 0; iblk < this->nblkl(ib); iblk++){
////				LBlock LB = this->L(ib)[iblk]; // IMPORTANT: Make a copy here rather than reference
////				if(bnum == LB.blkind()){
////					// Diagonal block computes the condition number
////					char norm = '1';
////					int N = LB.nrows();
////					int info;
////					double anorm = 1.0; // FIXME
////					double rcond;
////					DblNumVec work(4*N);
////					IntNumVec iwork(N);
////					iA( LB.nrows() == LB.ncols() );
////					dgecon_(&norm, &N, LB.nzval().data(), &N, &anorm, &rcond, 
////							work.data(), iwork.data(), &info);
////					double CondLB = 1.0 / rcond;
////					if( localValueMaxCond < 1.0 / rcond ){
////						localValueMaxCond = 1.0 / rcond;
////						localIndexMaxCond = bnum;
////					}
////				} 
////			} // for(iblk) 
////		}  
////	} 
////	// for(ib)
////	MPI_Barrier(comm);
////	MPI_Reduce(&localValueMaxCond, &globalValueMaxCond, 1, MPI_DOUBLE,
////			MPI_MAX, 0, comm);
////	std::cerr.precision(15);
////	if( mpirank == 0 ){
////		std::cerr << "Maximum condition number " << 
////			globalValueMaxCond << std::endl;
////	}
////	return 0;
////}
//
///**********************************************************************
// * Dump the diagonal elements in MATLAB format.  
// *
// * NOTE: The output is 1-based to be the same as the MATLAB format 
// **********************************************************************/
//int PMatrix::DumpDiagVec(NumVec<Scalar>& globalDiagVec, 
//		string filename, gridinfo_t *grid){
//	MPI_Comm comm = grid->comm;
//
//	int mpirank;  mpirank = grid->iam;
//	int mpisize;  mpisize = grid->nprow * grid->npcol;
//	int myrow = MYROW( grid->iam, grid ), mycol = MYCOL( grid->iam, grid );
//	iC( MPI_Barrier(comm) );
//
//
//	NumVec<Scalar> localDiagVec(this->ndof());
//	if( globalDiagVec.m() != this->ndof() )
//		globalDiagVec.resize(this->ndof());
//	setvalue(localDiagVec, SCALAR_ZERO);
//	setvalue(globalDiagVec, SCALAR_ZERO);
//
//	for(int ib = 0; ib < this->nbc(); ib++){
//		int bnum = ( (ib) * grid->npcol ) + mycol;
//		if( bnum >= this->nsupers() ) continue;
//		if( this->nblkl(ib) ){ // Not an empty column
//			for(int iblk = 0; iblk < this->nblkl(ib); iblk++){
//				LBlock& LB = this->L(ib)[iblk];
//				if(bnum == LB.blkind()){
//					// Diagonal block 
//					for(int i = 0; i < LB.nrows(); i++){
//						localDiagVec[LB.rows(i)] = LB.nzval(i,i);
//					} // for (i)
//				} 
//			} // for(iblk) 
//		}  
//	} 
//	// for(ib)
//
//#ifdef _USE_COMPLEX_
//	MPI_Reduce((double*)localDiagVec.data(), (double*)globalDiagVec.data(), 
//			this->ndof()*2, MPI_DOUBLE, MPI_SUM, 0, comm);
//#else
//	MPI_Reduce(localDiagVec.data(), globalDiagVec.data(), 
//			this->ndof(), MPI_DOUBLE, MPI_SUM, 0, comm);
//#endif
////	MPI_Reduce((void*)localDiagVec.data(), 
////			(void*)globalDiagVec.data(), 
////			this->ndof()*sizeof(Scalar), MPI_BYTE, MPI_SUM, 0, comm);
////	MPI_Reduce(localDiagVec.data(), globalDiagVec.data(), 
////			this->ndof(), MPI_DOUBLE, MPI_SUM, 0, comm);
//
//	if( mpirank == 0 ){
//		ofstream fid;
//		fid.open(filename.c_str(), ios::out | ios::trunc);
//		for(int i = 0; i < this->ndof(); i++){
//			fid << setw(25) << scientific << globalDiagVec[i] << std::endl;
//		}
//		fid.close();
//	}
//	MPI_Barrier(comm);
//	return 0;
//}
//
//// *********************************************************************
//// Convert PMatrix structure to a CSC matrix in order to perform trace
//// operations. 
//// *********************************************************************
//void PMatrixToCSCMatrix(int n, gridinfo_t *grid, PMatrix& PMloc){
//	// Each processor constructs the data structure for saving the local
//	// CSC format 
//	
//	// Each processor goes through the the L and U structure in the
//	// PMatrix it owns 
//	//
//	//
//	
//	// All column processors
//
//}
