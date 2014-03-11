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
/// @file new_ppexsi.cpp
/// @brief New implementation of the parallel %PEXSI.
/// @date 2014-03-09  Revise for the new interface.

// FIXME
#define _DEBUGlevel_ 0

#include "new_ppexsi.hpp"

namespace PEXSI{

PPEXSINewData::PPEXSINewData	( 
    MPI_Comm   comm,
    Int        numProcRow, 
    Int        numProcCol, 
    Int        outputFileIndex ){
#ifndef _RELEASE_
	PushCallStack("PPEXSINewData::PPEXSINewData");
#endif

  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );

  Int npPerPole = numProcRow * numProcCol;
  if( mpisize % npPerPole != 0 ){
    std::ostringstream msg;
    msg 
      << "mpisize    = " << mpisize << std::endl
      << "npPerPole = " << npPerPole << std::endl
      << "mpisize is not divisible by npPerPole!" << std::endl;
    throw std::runtime_error( msg.str().c_str() );
  }

  gridPole_     = new GridType( comm, mpisize / npPerPole, npPerPole );
	gridSuperLU_  = new SuperLUGrid( gridPole_->rowComm, 
      numProcRow, numProcCol );
	gridSelInv_   = new GridType( gridPole_->rowComm, 
      numProcRow, numProcCol );

  // Start the log file. Append to previous log files
  std::stringstream ss;
	ss << "logPEXSI" << outputFileIndex;
	statusOFS.open( ss.str().c_str(), std::ios_base::app );
  
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSINewData::PPEXSINewData  ----- 


PPEXSINewData::~PPEXSINewData	(  )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSINewData::~PPEXSINewData");
#endif
  if( gridPole_    != NULL ){
    delete gridPole_;
  }

	if( gridSuperLU_ != NULL ){
		delete gridSuperLU_;
	}
	
	if( gridSelInv_ != NULL ){
		delete gridSelInv_;
	}

  // Close the log file
  statusOFS.close();

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSINewData::~PPEXSINewData  ----- 


void
PPEXSINewData::LoadRealSymmetricMatrix	(
    Int           nrows,                        
    Int           nnz,                          
    Int           nnzLocal,                     
    Int           numColLocal,                  
    Int*          colptrLocal,                  
    Int*          rowindLocal,                  
    Real*         HnzvalLocal,                  
    Int           isSIdentity,                  
    Real*         SnzvalLocal )
{
#ifndef _RELEASE_
  PushCallStack("PPEXSINewData::LoadRealSymmetricMatrix");
#endif
  std::vector<char> sstr;
  Int sizeStm;
  if( MYROW( gridPole_ ) == 0 ){
    std::stringstream sstm;

    HRealMat_.size        = nrows;
    HRealMat_.nnz         = nnz;
    HRealMat_.nnzLocal    = nnzLocal;
    // The first row processor does not need extra copies of the index /
    // value of the matrix. 
    HRealMat_.colptrLocal = IntNumVec( numColLocal+1, false, colptrLocal );
    HRealMat_.rowindLocal = IntNumVec( nnzLocal,      false, rowindLocal );
    // H value
    HRealMat_.nzvalLocal  = DblNumVec( nnzLocal,      false, HnzvalLocal );
    HRealMat_.comm        = gridPole_->rowComm;

    // Serialization will copy the values regardless of the ownership
    serialize( HRealMat_, sstm, NO_MASK );

    // S value
    if( isSIdentity ){
      SRealMat_.size = 0;
      SRealMat_.nnz  = 0;
      SRealMat_.nnzLocal = 0;
      SRealMat_.comm = HRealMat_.comm; 
    }
    else{
      CopyPattern( HRealMat_, SRealMat_ );
      SRealMat_.comm = gridPole_->rowComm;
      SRealMat_.nzvalLocal  = DblNumVec( nnzLocal,      false, SnzvalLocal );
      serialize( SRealMat_.nzvalLocal, sstm, NO_MASK );
    }

    sstr.resize( Size( sstm ) );
    sstm.read( &sstr[0], sstr.size() ); 	
    sizeStm = sstr.size();
  }

  MPI_Bcast( &sizeStm, 1, MPI_INT, 0, gridPole_->colComm );

#if ( _DEBUGlevel_ >= 0 )
  statusOFS << "sizeStm = " << sizeStm << std::endl;
#endif

  if( MYROW( gridPole_ ) != 0 ) sstr.resize( sizeStm );

  MPI_Bcast( (void*)&sstr[0], sizeStm, MPI_BYTE, 0, gridPole_->colComm );

  if( MYROW( gridPole_ ) != 0 ){
    std::stringstream sstm;
    sstm.write( &sstr[0], sizeStm );
    deserialize( HRealMat_, sstm, NO_MASK );
    // Communicator
    HRealMat_.comm = gridPole_->rowComm;
    if( isSIdentity ){
      SRealMat_.size = 0;    // Means S is an identity matrix
      SRealMat_.nnz  = 0;
      SRealMat_.nnzLocal = 0;
      SRealMat_.comm = HRealMat_.comm;
    }
    else{
      CopyPattern( HRealMat_, SRealMat_ );
      SRealMat_.comm = gridPole_->rowComm;
      deserialize( SRealMat_.nzvalLocal, sstm, NO_MASK );
    }
  }
  sstr.clear();


#if ( _DEBUGlevel_ >= 0 )
  statusOFS << "H.size     = " << HRealMat_.size     << std::endl;
  statusOFS << "H.nnzLocal = " << HRealMat_.nnzLocal << std::endl;
  statusOFS << "S.size     = " << SRealMat_.size     << std::endl;
  statusOFS << "S.nnzLocal = " << SRealMat_.nnzLocal << std::endl;
#endif


#ifndef _RELEASE_
  PopCallStack();
#endif

  return ;
}    	// -----  end of method PPEXSINewData::LoadRealSymmetricMatrix  ----- 


void
PPEXSINewData::DFTDriver (
    Real       numElectronExact,
    Real       temperature,
    Real       gap,
    Real       deltaE,
    Int        numPole, 
    Int        isInertiaCount,
    Int        maxPEXSIIter,
    Real       muMin0,
    Real       muMax0,
    Real       muInertiaTolerance,
    Real       muPEXSISafeGuard,
    Real       numElectronPEXSITolerance,
    Int        matrixType,
    Int        ordering,
    Int        numProcSymbFact,
    Int        verbosity,
    Real&      muPEXSI,                   
    Real&      numElectronPEXSI,         
    Real&      muMinInertia,              
    Real&      muMaxInertia,             
    Int&       numTotalInertiaIter,   
    Int&       numTotalPEXSIIter )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSINewData::DFTDriver");
#endif

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSINewData::DFTDriver  ----- 

} //  namespace PEXSI
