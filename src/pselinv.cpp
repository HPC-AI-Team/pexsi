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

#if ( _DEBUGlevel_ >= 1 )
	statusOFS << std::endl << "PMatrix is constructed. The grid information: " << std::endl;
	statusOFS << "mpirank = " << MYPROC(grid_) << std::endl;
	statusOFS << "myrow   = " << MYROW(grid_) << std::endl; 
	statusOFS << "mycol   = " << MYCOL(grid_) << std::endl; 
#endif

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
		for( Int jsup = ksup + 1; jsup < numSuper; jsup++ ){
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
		for( Int isup = ksup + 1; isup < numSuper; isup++ ){
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
  
	for( Int ksup = 0; ksup < numSuper - 1; ksup++ ){
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
	PushCallStack("PMatrix::SelInv");
#endif
  Int numSuper = this->NumSuper(); 

	// Preparation for the selected inversion
	this->PreSelInv();

	// Main loop
	for( Int ksup = numSuper - 2; ksup >= 0; ksup-- ){
#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "ksup = " << ksup << std::endl << std::endl; 
#endif
#ifndef _RELEASE_
		PushCallStack("PMatrix::SelInv::UpdateL");
#endif
#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Communication to the Schur complement." << std::endl << std::endl; 
#endif

		// Communication for the U part.
#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Communication for the U part." << std::endl << std::endl; 
#endif
		std::vector<MPI_Request> mpireqsSendToDown( 2 * grid_->numProcRow, MPI_REQUEST_NULL ); 

		// Senders
		if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
			// Pack the data in U
			std::stringstream sstm;
			std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
			std::vector<UBlock>&  Urow = this->U( LBi(ksup, grid_) );
			// All blocks are to be sent down.
			serialize( (Int)Urow.size(), sstm, NO_MASK );
			for( Int jb = 0; jb < Urow.size(); jb++ ){
				serialize( Urow[jb], sstm, mask );
			}
			for( Int iProcRow = 0; iProcRow < grid_->numProcRow; iProcRow++ ){
				if( MYROW( grid_ ) != iProcRow &&
						isSendToDown_( ksup, iProcRow ) == true ){
					// Use Isend to send to multiple targets
					const std::string& sstr = sstm.str();
					Int sizeStm = sstr.length();
					MPI_Isend( &sizeStm, 1, MPI_INT,  iProcRow, grid_->numProcRow + iProcRow, 
							grid_->colComm, &mpireqsSendToDown[grid_->numProcRow + iProcRow] );
					MPI_Isend( (void*)sstr.c_str(), sizeStm, MPI_BYTE, iProcRow, iProcRow, 
							grid_->colComm, &mpireqsSendToDown[iProcRow] );
				} // Send 
			} // for (iProcRow)
		} // if I am the sender

		// Communication for the L part.
#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Communication for the L part." << std::endl << std::endl; 
#endif

		std::vector<MPI_Request> mpireqsSendToRight( 2 * grid_->numProcCol, MPI_REQUEST_NULL ); 

		// Senders
		if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
			// Pack the data in L 
			std::stringstream sstm;
			std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
			mask[LBlockMask::NZVAL] = 0; // nzval is excluded 

			std::vector<LBlock>&  Lcol = this->L( LBj(ksup, grid_) );
			// All blocks except for the diagonal block are to be sent right
			if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
				serialize( (Int)Lcol.size() - 1, sstm, NO_MASK );
			else
				serialize( (Int)Lcol.size(), sstm, NO_MASK );

			for( Int ib = 0; ib < Lcol.size(); ib++ ){
				if( Lcol[ib].blockIdx > ksup ){
					serialize( Lcol[ib], sstm, mask );
				}
			}
			
			for( Int iProcCol = 0; iProcCol < grid_->numProcCol ; iProcCol++ ){
				if( MYCOL( grid_ ) != iProcCol &&
						isSendToRight_( ksup, iProcCol ) == true ){
					// Use Isend to send to multiple targets
					const std::string& sstr = sstm.str();
					Int sizeStm = sstr.length();
					MPI_Isend( &sizeStm, 1, MPI_INT,  iProcCol, grid_->numProcCol + iProcCol, 
							grid_->rowComm, &mpireqsSendToRight[grid_->numProcCol + iProcCol] );
					MPI_Isend( (void*)sstr.c_str(), sizeStm, MPI_BYTE, iProcCol, iProcCol, 
							grid_->rowComm, &mpireqsSendToRight[iProcCol] );
				} // Send 
			} // for (iProcCol)
		} // if I am the sender

#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Receive from both up and from left." << std::endl << std::endl; 
#endif
    // Receive
		std::vector<MPI_Request> mpireqsRecvFromUp( 2, MPI_REQUEST_NULL ); 
		std::vector<MPI_Request> mpireqsRecvFromLeft( 2, MPI_REQUEST_NULL ); 
		Int sizeStmFromUp, sizeStmFromLeft;

		std::stringstream     sstmLcolRecv;
		std::stringstream     sstmUrowRecv;
		std::vector<char>     sstrLcolRecv;
		std::vector<char>     sstrUrowRecv;

		// Receive the size first

		if( isRecvFromUp_( ksup ) && 
				MYROW( grid_ ) != PROW( ksup, grid_ ) ){
			MPI_Irecv( &sizeStmFromUp, 1, MPI_INT, PROW( ksup, grid_ ), grid_->numProcRow + MYROW( grid_ ),
					grid_->colComm, &mpireqsRecvFromUp[0] );
		} // if I need to receive from up


		if( isRecvFromLeft_( ksup ) &&
				MYCOL( grid_ ) != PCOL( ksup, grid_ ) ){
			MPI_Irecv( &sizeStmFromLeft, 1, MPI_INT, PCOL( ksup, grid_ ), grid_->numProcCol + MYCOL( grid_ ),
					grid_->rowComm, &mpireqsRecvFromLeft[0] );
		} // if I need to receive from left

		// Wait to obtain size information
		mpi::Wait( mpireqsRecvFromUp[0] );
		mpi::Wait( mpireqsRecvFromLeft[0] );


		if( isRecvFromUp_( ksup ) && 
				MYROW( grid_ ) != PROW( ksup, grid_ ) ){
			sstrUrowRecv.resize( sizeStmFromUp );
			MPI_Irecv( &sstrUrowRecv[0], sizeStmFromUp, MPI_BYTE, PROW( ksup, grid_ ), 
					MYROW( grid_ ), grid_->colComm, 
					&mpireqsRecvFromUp[1] );
		} // if I need to receive from up

		if( isRecvFromLeft_( ksup ) &&
				MYCOL( grid_ ) != PCOL( ksup, grid_ ) ){
			sstrLcolRecv.resize( sizeStmFromLeft );
			MPI_Irecv( &sstrLcolRecv[0], sizeStmFromLeft, MPI_BYTE, PCOL( ksup, grid_ ), 
					MYCOL( grid_ ), grid_->rowComm,
					&mpireqsRecvFromLeft[1] );
		} // if I need to receive from left

		// Wait for all communication to finish
#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Wait for all communication to be done." << std::endl << std::endl; 
#endif
		// Wait to obtain packed information in a string and then write into stringstream
		mpi::Wait( mpireqsRecvFromUp[1] );
		mpi::Wait( mpireqsRecvFromLeft[1] );

		if( isRecvFromUp_( ksup ) && 
				MYROW( grid_ ) != PROW( ksup, grid_ ) ){
			sstmUrowRecv.write( &sstrUrowRecv[0], sizeStmFromUp );
			sstrUrowRecv.clear();
		} // if I need to receive from up 


		if( isRecvFromLeft_( ksup ) &&
				MYCOL( grid_ ) != PCOL( ksup, grid_ ) ){
			sstmLcolRecv.write( &sstrLcolRecv[0], sizeStmFromLeft );
			sstrLcolRecv.clear();
		} // if I need to receive from left

		// Overlap the communication with computation.  All processors move
		// to Gemm phase when ready 

#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Unpack the received data for processors participate in Gemm. " << std::endl << std::endl; 
#endif
		std::vector<LBlock>   LcolRecv;
		std::vector<UBlock>   UrowRecv;
		if( isRecvFromUp_( ksup ) && isRecvFromLeft_( ksup ) ){
			// U part
			if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
				std::vector<Int> mask( UBlockMask::TOTAL_NUMBER, 1 );
				Int numUBlock;
				deserialize( numUBlock, sstmUrowRecv, NO_MASK );
				UrowRecv.resize( numUBlock );
				for( Int jb = 0; jb < numUBlock; jb++ ){
					deserialize( UrowRecv[jb], sstmUrowRecv, mask );
				} 
			} // sender is not the same as receiver
			else{
				// U is obtained locally, just make a copy. Include everything
				// (there is no diagonal block)
				UrowRecv = this->U( LBi( ksup, grid_ ) );
			} // sender is the same as receiver

			if( MYCOL( grid_ ) != PCOL( ksup, grid_ ) ){
				std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
				mask[LBlockMask::NZVAL] = 0; // nzval is excluded
				Int numLBlock;
				deserialize( numLBlock, sstmLcolRecv, NO_MASK );
				LcolRecv.resize( numLBlock );
				for( Int ib = 0; ib < numLBlock; ib++ ){
					deserialize( LcolRecv[ib], sstmLcolRecv, mask );
				}
			} // sender is not the same as receiver
			else{
				// L is obtained locally, just make a copy. 
				// Do not include the diagonal block
				std::vector<LBlock>& Lcol =  this->L( LBj( ksup, grid_ ) );
				if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
					LcolRecv.resize( Lcol.size() );
					for( Int ib = 0; ib < Lcol.size(); ib++ ){
						LcolRecv[ib] = Lcol[ib];
					}
				}
				else{
					LcolRecv.resize( Lcol.size() - 1 );
					for( Int ib = 0; ib < Lcol.size() - 1; ib++ ){
						LcolRecv[ib] = Lcol[ib+1];
					}
				}
			} // sender is the same as receiver
		} // if I am a receiver


#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Main work: Gemm" << std::endl << std::endl; 
#endif

		// Save all the data to be updated for { L( isup, ksup ) | isup > ksup }.
		// The size will be updated in the Gemm phase and the reduce phase
		NumMat<Scalar>  LUpdateBuf;
    
		// Only the processors received information participate in the Gemm 
		if( isRecvFromUp_( ksup ) && isRecvFromLeft_( ksup ) ){
			// rowPtr[ib] gives the row index in LUpdateBuf for the first
			// nonzero row in LcolRecv[ib]. The total number of rows in
			// LUpdateBuf is given by rowPtr[end]-1
			std::vector<Int> rowPtr(LcolRecv.size() + 1);
			// colPtr[jb] gives the column index in UBuf for the first
			// nonzero column in UrowRecv[jb]. The total number of rows in
			// UBuf is given by colPtr[end]-1
			std::vector<Int> colPtr(UrowRecv.size() + 1);

			rowPtr[0] = 0;
			for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
				rowPtr[ib+1] = rowPtr[ib] + LcolRecv[ib].numRow;
			}
			colPtr[0] = 0;
			for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
				colPtr[jb+1] = colPtr[jb] + UrowRecv[jb].numCol;
			}

			Int numRowAinvBuf = *rowPtr.rbegin();
			Int numColAinvBuf = *colPtr.rbegin();

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << "AinvBuf ~ " << numRowAinvBuf << " x " << numColAinvBuf << std::endl;
			statusOFS << "rowPtr:" << std::endl << rowPtr << std::endl;
			statusOFS << "colPtr:" << std::endl << colPtr << std::endl;
#endif
			// Allocate for the computational storage
			NumMat<Scalar> AinvBuf( numRowAinvBuf, numColAinvBuf );
			LUpdateBuf.Resize( numRowAinvBuf, SuperSize( ksup, super_ ) );
			NumMat<Scalar> UBuf( SuperSize( ksup, super_ ), numColAinvBuf );
			SetValue( AinvBuf, SCALAR_ZERO );
			SetValue( LUpdateBuf, SCALAR_ZERO );
			SetValue( UBuf, SCALAR_ZERO );

			// Fill UBuf first.  Make the transpose later in the Gemm phase.
			for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
				UBlock& UB = UrowRecv[jb];
				if( UB.numRow != SuperSize(ksup, super_) ){
					throw std::logic_error( "The size of UB is not right.  Something is seriously wrong." );
				}
				lapack::Lacpy( 'A', UB.numRow, UB.numCol, UB.nzval.Data(),
					UB.numRow, UBuf.VecData( colPtr[jb] ), UBuf.m() );	
			}


			// Calculate the relative indices for (isup, jsup)
			// Fill AinvBuf with the information in L or U block.
			for( Int jb = 0; jb < UrowRecv.size(); jb++ ){
				for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
					LBlock& LB = LcolRecv[ib];
					UBlock& UB = UrowRecv[jb];
					Int isup = LB.blockIdx;
					Int jsup = UB.blockIdx;
					Scalar* nzvalAinv = &AinvBuf( rowPtr[ib], colPtr[jb] );
					Int     ldAinv    = AinvBuf.m();
					
					// Pin down the corresponding block in the part of Sinv.
					if( isup >= jsup ){
						std::vector<LBlock>&  LcolSinv = this->L( LBj(jsup, grid_ ) );
						bool isBlockFound = false;
						for( Int ibSinv = 0; ibSinv < LcolSinv.size(); ibSinv++ ){
							// Found the (isup, jsup) block in Sinv
							if( LcolSinv[ibSinv].blockIdx == isup ){
								LBlock& SinvB = LcolSinv[ibSinv];
								
								// Row relative indices
								std::vector<Int> relRows( LB.numRow );
								Int* rowsLBPtr    = LB.rows.Data();
								Int* rowsSinvBPtr = SinvB.rows.Data();
								for( Int i = 0; i < LB.numRow; i++ ){
									bool isRowFound = false;
									for( Int i1 = 0; i1 < SinvB.numRow; i1++ ){
										if( rowsLBPtr[i] == rowsSinvBPtr[i1] ){
											isRowFound = true;
											relRows[i] = i1;
											break;
										}
									}
									if( isRowFound == false ){
										std::ostringstream msg;
										msg << "Row " << rowsLBPtr[i] << 
											" in LB cannot find the corresponding row in SinvB" << std::endl
										 << "LB.rows    = " << LB.rows << std::endl
									   << "SinvB.rows = " << SinvB.rows << std::endl;
										throw std::runtime_error( msg.str().c_str() );
									}
								}

								// Column relative indicies
								std::vector<Int> relCols( UB.numCol );
                Int SinvColsSta = FirstBlockCol( jsup, super_ );
								for( Int j = 0; j < UB.numCol; j++ ){
									relCols[j] = UB.cols[j] - SinvColsSta;
								}

								statusOFS << "relRows = " << relRows << std::endl;
								statusOFS << "relCols = " << relCols << std::endl;

								// Transfer the values from Sinv to AinvBlock
								Scalar* nzvalSinv = SinvB.nzval.Data();
								Int     ldSinv    = SinvB.numRow;
								for( Int j = 0; j < UB.numCol; j++ ){
									for( Int i = 0; i < LB.numRow; i++ ){
										nzvalAinv[i+j*ldAinv] =
										 nzvalSinv[relRows[i] + relCols[j] * ldSinv];
									}
								}
								
								isBlockFound = true;
								break;
							}	
						} // for (ibSinv )
						if( isBlockFound == false ){
							std::ostringstream msg;
							msg << "Block(" << isup << ", " << jsup 
								<< ") did not find a matching block in Sinv." << std::endl;
							throw std::runtime_error( msg.str().c_str() );
						}
					} // if (isup, jsup) is in L
					else{
						std::vector<UBlock>&   UrowSinv = this->U( LBi( isup, grid_ ) );
						bool isBlockFound = false;
						for( Int jbSinv = 0; jbSinv < UrowSinv.size(); jbSinv++ ){
							// Found the (isup, jsup) block in Sinv
							if( UrowSinv[jbSinv].blockIdx == jsup ){
								UBlock& SinvB = UrowSinv[jbSinv];

								// Row relative indices
								std::vector<Int> relRows( LB.numRow );
								Int SinvRowsSta = FirstBlockCol( isup, super_ );
								for( Int i = 0; i < LB.numRow; i++ ){
									relRows[i] = LB.rows[i] - SinvRowsSta;
								}

								// Column relative indices
								std::vector<Int> relCols( UB.numCol );
								Int* colsUBPtr    = UB.cols.Data();
								Int* colsSinvBPtr = SinvB.cols.Data();
								for( Int j = 0; j < UB.numCol; j++ ){
									bool isColFound = false;
									for( Int j1 = 0; j1 < SinvB.numCol; j1++ ){
										if( colsUBPtr[j] == colsSinvBPtr[j1] ){
											isColFound = true;
											relCols[j] = j1;
											break;
										}
									}
									if( isColFound == false ){
										std::ostringstream msg;
										msg << "Col " << colsUBPtr[j] << 
											" in UB cannot find the corresponding row in SinvB" << std::endl
										 << "UB.cols    = " << UB.cols << std::endl
									   << "UinvB.cols = " << SinvB.cols << std::endl;
										throw std::runtime_error( msg.str().c_str() );
									}
								}

								// Trasnfer the values from Sinv to AinvBlock
								Scalar* nzvalSinv = SinvB.nzval.Data();
								Int     ldSinv    = SinvB.numRow;
								for( Int j = 0; j < UB.numCol; j++ ){
									for( Int i = 0; i < LB.numRow; i++ ){
										nzvalAinv[i+j*ldAinv] =
											nzvalSinv[relRows[i] + relCols[j] * ldSinv];
									}
								}

								isBlockFound = true;
								break;
							}
						} // for (jbSinv)
						if( isBlockFound == false ){
							std::ostringstream msg;
							msg << "Block(" << isup << ", " << jsup 
								<< ") did not find a matching block in Sinv." << std::endl;
							throw std::runtime_error( msg.str().c_str() );
						}
					} // if (isup, jsup) is in U

				} // for( ib )
			} // for ( jb )

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "AinvBuf: " << AinvBuf << std::endl;
      statusOFS << std::endl << "UBuf   : " << UBuf << std::endl;
#endif

			// Gemm for LUpdateBuf = AinvBuf * UBuf^T
			blas::Gemm( 'N', 'T', AinvBuf.m(), UBuf.m(), AinvBuf.n(), SCALAR_MINUS_ONE, 
					AinvBuf.Data(), AinvBuf.m(), 
					UBuf.Data(), UBuf.m(), SCALAR_ZERO,
					LUpdateBuf.Data(), LUpdateBuf.m() ); 

#if ( _DEBUGlevel_ >= 1 )
      statusOFS << std::endl << "LUpdateBuf: " << LUpdateBuf << std::endl;
#endif

		} // if Gemm is to be done locally

		// Now all the Isend / Irecv should have finished.
		mpi::Waitall( mpireqsSendToRight );
		mpi::Waitall( mpireqsSendToDown );

		// Reduce LUpdateBuf across all the processors in the same processor row.
		
#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Reduce LUpdateBuf to PCOL( ksup, grid_ )" << std::endl << std::endl; 
#endif

		NumMat<Scalar> LUpdateBufReduced;

		// Processor column of ksup collects the symbolic data for LUpdateBuf.
		std::vector<Int>  rowLocalPtr;
		Int numRowLUpdateBuf;
		if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
			std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
			if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
				rowLocalPtr.resize( Lcol.size() + 1 );
				rowLocalPtr[0] = 0;
				for( Int ib = 0; ib < Lcol.size(); ib++ ){
					rowLocalPtr[ib+1] = rowLocalPtr[ib] + Lcol[ib].numRow;
				}
			} // I do not own the diaogonal block
			else{
				rowLocalPtr.resize( Lcol.size() );
				rowLocalPtr[0] = 0;
				for( Int ib = 1; ib < Lcol.size(); ib++ ){
					rowLocalPtr[ib] = rowLocalPtr[ib-1] + Lcol[ib].numRow;
				}
			} // I owns the diagonal block, skip the diagonal block
			numRowLUpdateBuf = *rowLocalPtr.rbegin();
			if( numRowLUpdateBuf > 0 ){
				LUpdateBufReduced.Resize( numRowLUpdateBuf, SuperSize( ksup, super_ ) );
				SetValue( LUpdateBufReduced, SCALAR_ZERO );
			}
		} 

		// Processor column sends the total row dimension to all processors
		// in the same row to prepare for reduce
		MPI_Bcast( &numRowLUpdateBuf, 1, MPI_INT, PCOL( ksup, grid_ ), grid_->rowComm );

		// If LUpdatebuf has not been constructed, resize and fill with zero
		if( numRowLUpdateBuf > 0 ){
			if( LUpdateBuf.m() == 0 && LUpdateBuf.n() == 0 ){
				LUpdateBuf.Resize( numRowLUpdateBuf, SuperSize( ksup, super_ ) );
				// Fill zero is important
				SetValue( LUpdateBuf, SCALAR_ZERO );
			}

			mpi::Reduce( LUpdateBuf.Data(), LUpdateBufReduced.Data(),
					numRowLUpdateBuf * SuperSize( ksup, super_ ), MPI_SUM, 
					PCOL( ksup, grid_ ), grid_->rowComm );

		} // Perform reduce for nonzero block rows in the column of ksup

#if ( _DEBUGlevel_ >= 1 )
		if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) && numRowLUpdateBuf > 0 )
			statusOFS << std::endl << "LUpdateBufReduced: " << LUpdateBufReduced << std::endl << std::endl; 
#endif
		

#ifndef _RELEASE_
		PopCallStack();
#endif

#ifndef _RELEASE_
		PushCallStack("PMatrix::SelInv::UpdateD");
#endif
#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Update the diagonal block" << std::endl << std::endl; 
#endif
		if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
			NumMat<Scalar> DiagBuf( SuperSize( ksup, super_ ), SuperSize( ksup, super_ ) );
			SetValue( DiagBuf, SCALAR_ZERO );
			std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
			if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
				for( Int ib = 0; ib < Lcol.size(); ib++ ){
					blas::Gemm( 'T', 'N', SuperSize( ksup, super_ ), SuperSize( ksup, super_ ), Lcol[ib].numRow, 
							SCALAR_MINUS_ONE, &LUpdateBufReduced( rowLocalPtr[ib], 0 ), LUpdateBufReduced.m(),
							Lcol[ib].nzval.Data(), Lcol[ib].nzval.m(), 
							SCALAR_ONE, DiagBuf.Data(), DiagBuf.m() );
				}
			} // I do not own the diaogonal block
			else{
				for( Int ib = 1; ib < Lcol.size(); ib++ ){
					statusOFS << "LB = " << Lcol[ib].nzval << std::endl;
					blas::Gemm( 'T', 'N', SuperSize( ksup, super_ ), SuperSize( ksup, super_ ), Lcol[ib].numRow, 
							SCALAR_MINUS_ONE, &LUpdateBufReduced( rowLocalPtr[ib-1], 0 ), LUpdateBufReduced.m(),	
							Lcol[ib].nzval.Data(), Lcol[ib].nzval.m(), 
							SCALAR_ONE, DiagBuf.Data(), DiagBuf.m() );
				}
			} // I owns the diagonal block, skip the diagonal block

			NumMat<Scalar> DiagBufReduced( SuperSize( ksup, super_ ), SuperSize( ksup, super_ ) );

			if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
				SetValue( DiagBufReduced, SCALAR_ZERO );

			mpi::Reduce( DiagBuf.Data(), DiagBufReduced.Data(), 
					SuperSize( ksup, super_ ) * SuperSize( ksup, super_ ),
					MPI_SUM, PROW( ksup, grid_ ), grid_->colComm );
			
#if ( _DEBUGlevel_ >= 1 )
			statusOFS << std::endl << "DiagBuf: " << DiagBuf << std::endl << std::endl; 
#endif
			// Add DiagBufReduced to diagonal block.
			if( MYROW( grid_ ) == PROW( ksup, grid_ ) ){
#if ( _DEBUGlevel_ >= 1 )
				statusOFS << std::endl << "DiagBufReduced: " << DiagBufReduced << std::endl << std::endl; 
#endif
        LBlock&  LB = this->L( LBj( ksup, grid_ ) )[0];
				blas::Axpy( LB.numRow * LB.numCol, 1.0, DiagBufReduced.Data(),
						1, LB.nzval.Data(), 1 );
#if ( _DEBUGlevel_ >= 1 )
				statusOFS << std::endl << "Diag of Ainv: " << LB.nzval << std::endl << std::endl; 
#endif
			}
		} // Update the diagonal in the processor column of ksup. All processors participate



#ifndef _RELEASE_
		PopCallStack();
#endif

#ifndef _RELEASE_
		PushCallStack("PMatrix::SelInv::UpdateU");
#endif
#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Update the upper triangular block" << std::endl << std::endl; 
#endif
		// Send LUpdateBufReduced to the cross diagonal blocks. 
		// NOTE: This assumes square processor grid
	  if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) && isSendToCrossDiagonal_( ksup ) ){
			Int dest = PNUM( PROW( ksup, grid_ ), MYROW( grid_ ), grid_ );
			if( MYPROC( grid_ ) != dest	){
				std::stringstream sstm;
				serialize( LUpdateBufReduced, sstm, NO_MASK );
				mpi::Send( sstm, dest, ksup + numSuper, ksup, grid_->comm );
			}
		} // sender

		if( MYROW( grid_ ) == PROW( ksup, grid_ ) && isRecvFromCrossDiagonal_( ksup ) ){
			Int src = PNUM( MYCOL( grid_ ), PCOL( ksup, grid_ ), grid_ );
			NumMat<Scalar> UUpdateBuf;
			if( MYPROC( grid_ ) != src ){
				std::stringstream sstm;
				mpi::Recv( sstm, src, ksup + numSuper, ksup, grid_->comm );
			  deserialize( UUpdateBuf, sstm, NO_MASK );	
			} // sender is not the same as receiver
			else{
				UUpdateBuf = LUpdateBufReduced;
			} // sender is the same as receiver
#if ( _DEBUGlevel_ >= 1 )
			statusOFS << std::endl << "UUpdateBuf:" << UUpdateBuf << std::endl << std::endl; 
#endif

			// Update U
			std::vector<UBlock>&  Urow = this->U( LBi( ksup, grid_ ) );
			Int cntRow = 0;
			for( Int jb = 0; jb < Urow.size(); jb++ ){
				UBlock& UB = Urow[jb];
				NumMat<Scalar> Ltmp( UB.numCol, UB.numRow );
				lapack::Lacpy( 'A', Ltmp.m(), Ltmp.n(), &UUpdateBuf( cntRow, 0 ),
						UUpdateBuf.m(), Ltmp.Data(), Ltmp.m() );

				Transpose( Ltmp, UB.nzval );
				cntRow += UB.numCol;
			} // for (jb)
			if( cntRow != UUpdateBuf.m() ){
				std::ostringstream msg;
				msg << "The number of rows received from L is " << UUpdateBuf.m()
					<< ", which does not match the total number of columns in U which is "
					<< cntRow <<  std::endl;
				throw std::runtime_error( msg.str().c_str() );
			}
		} // receiver

#ifndef _RELEASE_
		PopCallStack();
#endif

#ifndef _RELEASE_
		PushCallStack("PMatrix::SelInv::UpdateLFinal");
#endif
#if ( _DEBUGlevel_ >= 1 )
		statusOFS << std::endl << "Finish updating the L part by filling LUpdateBufReduced back to L" << std::endl << std::endl; 
#endif
#ifndef _RELEASE_
		PopCallStack();
#endif
		if( MYCOL( grid_ ) == PCOL( ksup, grid_ ) && numRowLUpdateBuf > 0 ){
			std::vector<LBlock>&  Lcol = this->L( LBj( ksup, grid_ ) );
			if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
        for( Int ib = 0; ib < Lcol.size(); ib++ ){
					LBlock& LB = Lcol[ib];
					lapack::Lacpy( 'A', LB.numRow, LB.numCol, &LUpdateBufReduced( rowLocalPtr[ib], 0 ),
							LUpdateBufReduced.m(), LB.nzval.Data(), LB.numRow );
				}
			} // I do not own the diagonal block
			else{
				for( Int ib = 1; ib < Lcol.size(); ib++ ){
					LBlock& LB = Lcol[ib];
					lapack::Lacpy( 'A', LB.numRow, LB.numCol, &LUpdateBufReduced( rowLocalPtr[ib-1], 0 ),
							LUpdateBufReduced.m(), LB.nzval.Data(), LB.numRow );
				}
			} // I owns the diagonal block
		} // Finish updating L	

    MPI_Barrier( grid_-> comm );
	} // for (ksup) : Main loop

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PMatrix::SelInv  ----- 


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

		// Sender
		if( isSendToCrossDiagonal_[ksup]  &&
			  MYPROC( grid_ ) !=  PNUM( ksupProcRow, MYROW( grid_ ), grid_ ) ){
      // Pack L data
			std::stringstream sstm;
			std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
			std::vector<LBlock>&  Lcol = this->L( LBj(ksup, grid_) );
			// All blocks except for the diagonal block are to be sent right
			if( MYROW( grid_ ) == PROW( ksup, grid_ ) )
				serialize( (Int)Lcol.size() - 1, sstm, NO_MASK );
			else
				serialize( (Int)Lcol.size(), sstm, NO_MASK );

			for( Int ib = 0; ib < Lcol.size(); ib++ ){
				if( Lcol[ib].blockIdx > ksup ) 
					serialize( Lcol[ib], sstm, mask );
			}
			// Send/Recv is possible here due to the one to one correspondence
			// in the case of square processor grid
			mpi::Send( sstm, PNUM( ksupProcRow, MYROW( grid_ ), grid_ ), 
					ksup + numSuper, ksup, grid_->comm );
		} // if I am a sender

		// Receiver
		if( isRecvFromCrossDiagonal_[ksup] ){

			std::vector<LBlock> LcolRecv;
			if( PNUM( MYCOL( grid_ ), ksupProcCol, grid_ ) != MYPROC( grid_ ) ){
				std::stringstream sstm;
				mpi::Recv( sstm, PNUM( MYCOL( grid_ ), ksupProcCol, grid_ ), 
						ksup + numSuper, ksup, grid_->comm );
				
				// Unpack L data.  
				Int numLBlock;
				std::vector<Int> mask( LBlockMask::TOTAL_NUMBER, 1 );
				deserialize( numLBlock, sstm, NO_MASK );
				LcolRecv.resize(numLBlock);
				for( Int ib = 0; ib < numLBlock; ib++ ){
					deserialize( LcolRecv[ib], sstm, mask );
				}
			} // sender is not the same as receiver
			else{
				// L is obtained locally, just make a copy. Do not include the diagonal block
				std::vector<LBlock>& Lcol = this->L( LBj( ksup, grid_ ) );
				if( MYROW( grid_ ) != PROW( ksup, grid_ ) ){
					LcolRecv.resize( Lcol.size() );
					for( Int ib = 0; ib < Lcol.size(); ib++ ){
						LcolRecv[ib] = Lcol[ib];
					}
				}
				else{
					LcolRecv.resize( Lcol.size() - 1 );
					for( Int ib = 0; ib < Lcol.size() - 1; ib++ ){
						LcolRecv[ib] = Lcol[ib+1];
					}
				}
			} // sender is the same as receiver
				
			// Update U
			// Make sure that the size of L and the corresponding U blocks match.
			std::vector<UBlock>& Urow = this->U( LBi( ksup, grid_ ) );
			for( Int ib = 0; ib < LcolRecv.size(); ib++ ){
				LBlock& LB = LcolRecv[ib];
				if( LB.blockIdx <= ksup ){
					throw std::logic_error( "LcolRecv contains the wrong blocks." );
				}
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
			} // for (ib)
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

	for( Int ksup = 0; ksup < numSuper; ksup++ ){
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
			statusOFS << "Factorized A (" << ksup << ", " << ksup << "): " << LB.nzval << std::endl;
#endif
			lapack::Getri( SuperSize( ksup, super_ ), LB.nzval.Data(), 
					SuperSize( ksup, super_ ), ipiv.Data() );
#if ( _DEBUGlevel_ >= 1 )
			// Check the correctness of the matrix inversion for the first local column
			statusOFS << "Inversed   A (" << ksup << ", " << ksup << "): " << LB.nzval << std::endl;
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
PMatrix::Diagonal	( NumVec<Scalar>& diagNaturalOrder )
{
#ifndef _RELEASE_
	PushCallStack("PMatrix::Diagonal");
#endif
  Int numSuper = this->NumSuper(); 

  Int numCol = this->NumCol();

	NumVec<Scalar> diagLocal( numCol );
	SetValue( diagLocal, SCALAR_ZERO );
	NumVec<Scalar> diag( numCol );
	SetValue( diag, SCALAR_ZERO );

	diagNaturalOrder.Resize( numCol );
	SetValue( diagNaturalOrder, SCALAR_ZERO );


	for( Int ksup = 0; ksup < numSuper; ksup++ ){
	  // I own the diagonal block	
		if( MYROW( grid_ ) == PROW( ksup, grid_ ) &&
				MYCOL( grid_ ) == PCOL( ksup, grid_ ) ){
			LBlock& LB = this->L( LBj( ksup, grid_ ) )[0];
			statusOFS << "LB(" << ksup << ", " << ksup << ").rows = " << LB.rows << std::endl;
      for( Int i = 0; i < LB.numRow; i++ ){
        diagLocal( LB.rows(i) ) = LB.nzval( i, i );
			}
		}
	}

	// All processors own diag
	mpi::Allreduce( diagLocal.Data(), diag.Data(), numCol, MPI_SUM, grid_->comm );

	// Permute diag back to natural order
	for( Int i = 0; i < numCol; i++ ){
		diagNaturalOrder(i) = diag( super_->perm(i) );
	}

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PMatrix::Diagonal  ----- 

} // namespace PEXSI
