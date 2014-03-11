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
/// @file new_interface.cpp
/// @brief New interface subroutines of PPEXSI that can be called by both C and FORTRAN.
///
/// This file will eventually merge with interface.cpp.
/// @date 2014-03-07
#include "c_pexsi_interface.h"
#include "c_pexsi_new_interface.h"
#include "new_ppexsi.hpp"
#include "blas.hpp"

// FIXME
// Error handling used in the C interface that is different from the
// throw/catch system.
#define iC(fun)  { int ierr=fun; if(ierr!=0) exit(1); }
#define iA(expr) { if((expr)==0) { std::cerr<<"wrong "<<__LINE__<<" in " <<__FILE__<<std::endl; std::cerr.flush(); exit(1); } }

using namespace PEXSI;

extern "C"
void PPEXSISetDefaultOptions(
    PPEXSIOptions*   options ){
  options->temperature           = 0.0019;   // 300K 
  options->gap                   = 0.0;      // no gap 
  options->deltaE                = 10.0; 
  options->numPole               = 40;
  options->isInertiaCount        = 1;
  options->maxPEXSIIter          = 3;
  options->muMin0                = -10.0; 
  options->muMax0                = +10.0; 
  options->muInertiaTolerance    = 0.05;
  options->muPEXSISafeGuard      = 0.05;
  options->numElectronPEXSITolerance = 0.01;
  options->matrixType            = 0;
  options->ordering              = 0;
  options->npSymbFact            = 1;
  options->verbosity             = 1;
}   // -----  end of function PPEXSISetDefaultOptions  ----- 


extern "C"
PPEXSIPlan PPEXSIPlanInitialize(
    MPI_Comm      comm,
    int           numProcRow,
    int           numProcCol,
    int           outputFileIndex, 
    int*          info ){

  Int mpirank, mpisize;
  MPI_Comm_rank( comm, &mpirank );
  MPI_Comm_size( comm, &mpisize );

  *info = 0;
  PPEXSINewData *ptrData;

  try{
    ptrData = new PPEXSINewData( comm, numProcRow, numProcCol, mpirank );
  }
	catch( std::exception& e )
	{
		statusOFS << std::endl << "ERROR!!! Proc " << mpirank << " caught exception with message: "
			<< std::endl << e.what() << std::endl;
		*info = 1;
#ifndef _RELEASE_
		DumpCallStack();
#endif
	}

  return reinterpret_cast<PPEXSIPlan>(ptrData);
}   // -----  end of function PPEXSIPlanInitialize  ----- 


extern "C"
void PPEXSILoadRealSymmetricHSMatrix(
    PPEXSIPlan    plan,
    int           nrows,                        
    int           nnz,                          
    int           nnzLocal,                     
    int           numColLocal,                  
    int*          colptrLocal,                  
    int*          rowindLocal,                  
    double*       HnzvalLocal,                  
    int           isSIdentity,                  
    double*       SnzvalLocal, 
    int*          info ){

  const GridType* gridPole = 
    reinterpret_cast<PPEXSINewData*>(plan)->GridPole();

  *info = 0;

  try{
    reinterpret_cast<PPEXSINewData*>(plan)->
      LoadRealSymmetricMatrix(
          nrows,                        
          nnz,                          
          nnzLocal,                     
          numColLocal,                  
          colptrLocal,                  
          rowindLocal,                  
          HnzvalLocal,                  
          isSIdentity,                  
          SnzvalLocal );
  }
	catch( std::exception& e )
	{
		statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank 
      << " caught exception with message: "
			<< std::endl << e.what() << std::endl;
		*info = 1;
#ifndef _RELEASE_
		DumpCallStack();
#endif
	}

  return;
}   // -----  end of function PPEXSILoadRealSymmetricHSMatrix  ----- 


extern "C"
void PPEXSIDFTDriver(
    /* Input parameters */
    PPEXSIPlan        plan,
    double            numElectronExact,
    PPEXSIOptions     options,
    /* Output parameters */
		double*           muPEXSI,                   
		double*           numElectronPEXSI,         
    double*           muMinInertia,              
		double*           muMaxInertia,             
		int*              numTotalInertiaIter,   
		int*              numTotalPEXSIIter,   
    int*              info ){
  *info = 0;
  const GridType* gridPole = 
    reinterpret_cast<PPEXSINewData*>(plan)->GridPole();
  
  try{
    reinterpret_cast<PPEXSINewData*>(plan)->DFTDriver(
        numElectronExact,
        options.temperature,
        options.gap,
        options.deltaE,
        options.numPole,
        options.isInertiaCount,
        options.maxPEXSIIter,
        options.muMin0,
        options.muMax0,
        options.muInertiaTolerance,
        options.muPEXSISafeGuard,
        options.numElectronPEXSITolerance,
        options.matrixType,
        options.ordering,
        options.npSymbFact,
        options.verbosity,
        *muPEXSI,
        *numElectronPEXSI,
        *muMinInertia,
        *muMaxInertia,
        *numTotalInertiaIter,
        *numTotalPEXSIIter );
  }
	catch( std::exception& e )
	{
		statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank 
      << " caught exception with message: "
			<< std::endl << e.what() << std::endl;
		*info = 1;
#ifndef _RELEASE_
		DumpCallStack();
#endif
	}
  return;

  return;
}   // -----  end of function PPEXSIDFTDriver  ----- 


//extern "C"
//void PPEXSIRetrieveRealSymmetricDFTMatrix(
//    PPEXSIPlan        plan,
//		double*      DMnzvalLocal,
//		double*     EDMnzvalLocal,
//		double*     FDMnzvalLocal,
//    int*              info ){
//  *info = 0;
//
//	// Synchronize the info among all processors. 
//	// If any processor gets error message, info = 1
//	Int infoAll = 0;
//	mpi::Allreduce( info, &infoAll, 1, MPI_MAX, comm  );
//	*info = infoAll;
//
//  return;
//}


extern "C"
void PPEXSIPlanFinalize( 
    PPEXSIPlan    plan,
    int*          info ){

  *info = 0;
  const GridType* gridPole = 
    reinterpret_cast<PPEXSINewData*>(plan)->GridPole();
  
  try{
    delete reinterpret_cast<PPEXSINewData*>(plan);
  }
	catch( std::exception& e )
	{
		statusOFS << std::endl << "ERROR!!! Proc " << gridPole->mpirank 
      << " caught exception with message: "
			<< std::endl << e.what() << std::endl;
		*info = 1;
#ifndef _RELEASE_
		DumpCallStack();
#endif
	}
  return;
}   // -----  end of function PPEXSIPlanFinalize  ----- 
