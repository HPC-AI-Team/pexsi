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
    Real       muInertiaExpansion,
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
  numTotalInertiaIter = 0;
  numTotalPEXSIIter   = 0;

	Real timeSta, timeEnd;

  // Initial setup
  muMinInertia = muMin0;
  muMaxInertia = muMax0;

  std::string colPerm;
  switch (ordering){
    case 0:
      colPerm = "PARMETIS";
      break;
    case 1:
      colPerm = "METIS_AT_PLUS_A";
      break;
    case 2:
      colPerm = "MMD_AT_PLUS_A";
      break;
    default:
      throw std::logic_error("Unsupported ordering strategy.");
  }


  // Main loop: inertia count + PEXSI
   
  while(1){

    // Perform inertia count
    if( isInertiaCount == 1 ){
      bool     isBadBound = false;

      // Number of shifts is exactly determined by the number of
      // independent groups to minimize the cost
      // However, the minimum number of shifts is 10 to accelerate
      // convergence.
      Int numShift = std::max( gridPole_->numProcRow, 10 );
      std::vector<Real>  shiftVec( numShift );
      std::vector<Real>  inertiaVec( numShift );   // Zero temperature
      std::vector<Real>  inertiaFTVec( numShift ); // Finite temperature
      Int iter;
      Int maxInertiaIter = std::max( 1, iround( 
            std::log( (muMax0 - muMin0) / muInertiaTolerance ) /
            std::log( numShift ) ) ); 

      if( verbosity == 1 ){
        statusOFS << std::endl
          << "From" << std::endl
          << "(muMin0, muMax0)   = " << "(" << muMin0 << ", " << muMax0 
          << ")" << std::endl
          << "muInertiaTolerance = " << muInertiaTolerance << std::endl
          << "numShift           = " << numShift << std::endl
          << "we have " << std::endl
          << "maxInertiaIter     = " << maxInertiaIter << std::endl << std::endl;
      }

      for( iter = 1; iter <= maxInertiaIter; iter++ ){
        for( Int l = 0; l < numShift; l++ ){
          shiftVec[l] = muMinInertia + l * (muMaxInertia - muMinInertia) / (numShift-1);
        }


        GetTime( timeSta );
        CalculateNegativeInertiaReal(
            shiftVec,
            inertiaVec,
            HRealMat_,
            SRealMat_,
            colPerm,
            numProcSymbFact );

        GetTime( timeEnd );

        // Inertia is multiplied by 2.0 to reflect the doubly occupied
        // orbitals.
        for( Int l = 0; l < numShift; l++ ){
          inertiaVec[l] *= 2.0;
        }

        // Using linear interpolation procedure to compute the finite
        // temperature cumulative DOS
        {
          Int numInp = 1000;   // Number of interpolation points
          Real shiftExpand = 20 * temperature;   // Expand the interval for interpolation
          std::vector<Real>  shiftInpVec( numInp );
          std::vector<Real>  inertiaInpVec( numInp ); 
          std::vector<Real>  fdDrvVec( numInp );
          for( Int l = 0; l < numShift; l++ ){
            Real shiftInp0 = shiftVec[l] - shiftExpand;
            Real shiftInp1 = shiftVec[l] + shiftExpand; 
            Real h = (shiftInp1 - shiftInp0) / (numInp-1);
            for( Int i = 0; i < numInp; i++ ){
              shiftInpVec[i] = shiftInp0 + h * i;
              // fdDrvMu(x) = beta * exp(beta*x)/(1+exp(beta*z))^2
              // Note: compared to fdDrvMu used in pole.cpp, the factor 2.0 is not
              // present, because it is included in inertiaVec. 
              fdDrvVec[i]    = 1.0 / ( 2.0 * temperature * (
                    1.0 + std::cosh( ( shiftInpVec[i] - shiftVec[l] ) / temperature ) ) );
            }
            LinearInterpolation( shiftVec, inertiaVec, 
                shiftInpVec, inertiaInpVec );

            inertiaFTVec[l] = 0.0;
            for( Int i = 0; i < numInp; i++ ){
              inertiaFTVec[l] += fdDrvVec[i] * inertiaInpVec[i] * h;
            }

            if( verbosity == 1 ){
              statusOFS << std::setiosflags(std::ios::left) 
                << std::setw(LENGTH_VAR_NAME) << "Shift      = "
                << std::setw(LENGTH_VAR_DATA) << shiftVec[l]
                << std::setw(LENGTH_VAR_NAME) << "Inertia    = "
                << std::setw(LENGTH_VAR_DATA) << inertiaVec[l]
                << std::setw(LENGTH_VAR_NAME) << "InertiaFT  = "
                << std::setw(LENGTH_VAR_DATA) << inertiaFTVec[l]
                << std::endl;
            }
          } // for (l)
        }

        if( verbosity == 1 ){
          statusOFS << std::endl << "Time for iteration " 
            << iter << " of the inertia count is " 
            << timeEnd - timeSta << std::endl;
        }

        // If the chemical potential does not fall into the
        // (muMinInertia, muMaxInertia) bracket, expand the interval.

        if( inertiaFTVec[0] > numElectronExact - 0.1 ){
          isBadBound = true;
          muMaxInertia = muMinInertia;
          muMinInertia = muMinInertia - muInertiaExpansion;
        }

        if( inertiaFTVec[numShift-1] < numElectronExact + 0.1 ){
          isBadBound = true;
//          muMinInertia = 
        }

        if( inertiaFTVec[0] >= numElectronExact ||
            inertiaFTVec[numShift-1] <= numElectronExact ){
          std::ostringstream msg;
          msg 
            << "The solution is not in the prescribed (muMin, muMax) interval." << std::endl
            << "(muMin, muMax) ~ (" << shiftVec[0] << " , " << shiftVec[numShift-1] << " ) " << std::endl
            << "(Ne(muMin), Ne(muMax)) ~ (" << inertiaFTVec[0] << " , " << inertiaFTVec[numShift-1] 
            << " ) " << std::endl
            << "NeExact = " << numElectronExact << std::endl
            << "Try to increase numElectronTolerance or initial search interval for mu." << std::endl;
          throw std::runtime_error( msg.str().c_str() );
        }
//
//
//        // Update muMin, muMax
//        // NOTE: divide by 4 is just heuristics so that the interval does
//        // not get to be too small for the PEXSI iteration
//
//        // First find the smallest interval
//        const Real EPS = 1e-1;
//        Int idxMin = 0, idxMax = numShift-1;
//        for( Int i = 0; i < numShift; i++ ){
//          if( ( inertiaFTVec[i] < numElectronExact ) &&
//              ( inertiaVec[i]   < numElectronExact - EPS ) ){
//            idxMin = ( idxMin < i ) ? i : idxMin;
//          }
//          if( ( inertiaFTVec[i] > numElectronExact ) &&
//              ( inertiaVec[i]   > numElectronExact + EPS ) ){
//            idxMax = ( idxMax > i ) ? i : idxMax;
//          }
//        }
//
//#if ( _DEBUGlevel_ >= 0 )
//        statusOFS << "idxMin = " << idxMin << "inertiaVec = " << inertiaVec[idxMin] << std::endl;
//        statusOFS << "idxMax = " << idxMax << "inertiaVec = " << inertiaVec[idxMax] << std::endl;
//#endif
//
//
//
//        if( inertiaFTVec[idxMax] - inertiaFTVec[idxMin] >= 
//            numElectronTolerance ){
//          // The smallest interval still contain too many eigenvalues
//          muMin = shiftVec[idxMin];
//          muMax = shiftVec[idxMax];
//          Real NeMin, NeMax;
//          NeMin = std::min( inertiaFTVec[idxMin], 
//              numElectronExact - numElectronTolerance / 2 + EPS );
//          NeMin = std::max( inertiaFTVec[0] + EPS, NeMin );
//          NeMax = std::max( inertiaFTVec[idxMax],
//              numElectronExact + numElectronTolerance / 2 - EPS );
//          NeMax = std::min( inertiaFTVec[numShift-1] - EPS, NeMax );
//          muMin = MonotoneRootFinding( shiftVec, inertiaFTVec, NeMin );
//          muMax = MonotoneRootFinding( shiftVec, inertiaFTVec, NeMax );
//        }
//        else{
//          // Root finding using linear interpolation
//          Real NeMin, NeMax;
//          NeMin = numElectronExact - numElectronTolerance / 2 + EPS;
//          NeMin = std::max( inertiaFTVec[0] + EPS, NeMin );
//          NeMax = numElectronExact + numElectronTolerance / 2 - EPS;
//          NeMax = std::min( inertiaFTVec[numShift-1] - EPS, NeMax );
//          muMin = MonotoneRootFinding( shiftVec, inertiaFTVec, NeMin );
//          muMax = MonotoneRootFinding( shiftVec, inertiaFTVec, NeMax );
//        }
//
//
//        double mu0 = (muMin + muMax)/2.0;
//
//
//        Real NeLower = numElectronExact - EPS;
//        Real NeUpper = numElectronExact + EPS;
//
//
//
//        muLow = MonotoneRootFinding( shiftVec, inertiaFTVec, NeLower );
//        muUpp = MonotoneRootFinding( shiftVec, inertiaFTVec, NeUpper );
//
//        // Convert the internal variables to output parameters
//        // This early output is mainly for debugging purpose.
//        {
//          *muMinInertia = muMin;
//          *muMaxInertia = muMax;
//          *muLowerEdge  = muLow;
//          *muUpperEdge  = muUpp;
//          *numIter      = iter;
//        }
//
//        // Output some information after EACH STEP for dynamical adjustament.
//        {
//#if ( _DEBUGlevel_ >= 0 )
//          statusOFS << std::endl << "After inertia count iteration " << iter
//            << std::endl;
//          Print( statusOFS, "muLowerEdge   = ", muLow );
//          Print( statusOFS, "muUppperEdge  = ", muUpp );
//          Print( statusOFS, "muMin         = ", muMin );
//          Print( statusOFS, "muMax         = ", muMax );
//          statusOFS << std::endl;
//#endif
//        }
//
//
//        if( inertiaFTVec[numShift-1] - inertiaFTVec[0] < numElectronTolerance ){
//          isConverged = true;
//          break;
//        }
      } // for (iter)
    }
    // FIXME
    break;
  }

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSINewData::DFTDriver  ----- 

} //  namespace PEXSI
