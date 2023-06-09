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
/// @file run_fermi_complex.cpp
/// @brief Calculate the Fermi operator given complex Hermitian H and S matrices.
///
/// @date 2016-09-05  Original version.
/// @date 2018-06-14  Compatible with the interface at version 1.0

#include "ppexsi.hpp"
#include "c_pexsi_interface.h"

using namespace PEXSI;
using namespace std;


void Usage(){
  std::cout 
    << "run_fermi" << std::endl;
}

int main(int argc, char **argv) 
{
  MPI_Init(&argc, &argv);
  int mpirank, mpisize;
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

  // PEXSI parameters
  int           isSIdentity;                  

  double        numElectronExact;

  int           numColLocal;

  char*         HfileR;
  char*         HfileI;
  char*         Hfile;
  char*         Sfile;
  int           isFormatted;

  PPEXSIOptions options;

  int           nprow, npcol;
  MPI_Comm      readComm;
  int           isProcRead;
  int           outputFileIndex;

  Real          numElectron;
  Real          numElectronDrvMu;

  std::string colPerm, rowPerm;

  try{
    // *********************************************************************
    // Input parameter
    // *********************************************************************

    if(1){
      numElectronExact    = 12.0;
      nprow               = 1;
      npcol               = 1;
      Hfile              = "lap2dc.matrix";
      //HfileI              = "lap2dc_imag.matrix";
      Sfile               = "";
      isFormatted         = 1;
      isSIdentity         = 1;
      colPerm             = "PARMETIS";
      rowPerm             = "NOROWPERM";
    }
    if(0){
      // Si64k matrix from SIESTA
      numElectronExact    = 256.0;
      nprow               = 1;
      npcol               = 1;
      Hfile               = "/data/matrices/Si64k/Hk.00002.matrix";
      Sfile               = "/data/matrices/Si64k/Sk.00002.matrix";
      isFormatted         = 1;
      isSIdentity         = 0;
      colPerm             = "PARMETIS";
      rowPerm             = "NOROWPERM";
    }


    /* Split the processors to read matrix */
    if( mpirank < nprow * npcol )
      isProcRead = 1;
    else
      isProcRead = 0;

    MPI_Comm_split( MPI_COMM_WORLD, isProcRead, mpirank, &readComm );

    DistSparseMatrix<Complex> HMat;
    DistSparseMatrix<Complex> SMat;

    if( isProcRead == 1 ){
      printf("Proc %5d is reading file...\n", mpirank );
      /* Read the matrix head for allocating memory */
      ReadDistSparseMatrixFormatted( Hfile, 
          HMat, readComm ); 

      numColLocal = HMat.colptrLocal.m() - 1;

      if( mpirank == 0 ){
        printf("On processor 0...\n");
        printf("nrows       = %d\n", HMat.size );
        printf("nnz         = %d\n", HMat.nnz );
        printf("nnzLocal    = %d\n", HMat.nnzLocal );
        printf("numColLocal = %d\n", numColLocal );
      }

      if( isSIdentity == 0 ){
        if( isFormatted == 1 ){
          ReadDistSparseMatrixFormatted( Sfile, SMat, readComm ); 
        }
//        else{
//          ReadDistSparseMatrix( Sfile, SMat, readComm ); 
//        }
      }

      if( mpirank == 0 ){ 
        printf("Finish reading the matrix.\n");
      }
    } // Read the matrix


    // *********************************************************************
    // Check the input parameters
    // *********************************************************************

    // Initialize

    /* Set the outputFileIndex to be the pole index */
    /* The first processor for each pole outputs information */

    if( mpirank % (nprow * npcol) == 0 ){
      outputFileIndex = mpirank / (nprow * npcol);
    }
    else{
      outputFileIndex = -1;
    }

    Int npPerPole = nprow * npcol;
    PPEXSIData pexsi( MPI_COMM_WORLD, nprow, npcol, outputFileIndex );


    /* Step 1. Initialize PEXSI */

    PPEXSISetDefaultOptions( &options );
    options.muMin0 = -1.0;
    options.muMax0 =  1.0;
    options.mu0    = +0.300;
    options.npSymbFact = 1;
    options.ordering = 0;
    options.solver = 0;
    options.isInertiaCount = 1;
    options.verbosity = 1;
    options.deltaE   = 10.0;
    options.numElectronPEXSITolerance = 0.00001;
    options.isSymbolicFactorize = 1;
    options.nPoints = 2;
    options.spin = 2.0;

    options.method = 1;
    options.numPole  = 40;
    options.temperature  = 0.00095; // 300K

    /*
    FILE * fp;
    fp = fopen("input.txt", "r");
    int temp;
    double temp1;
    int method;
    rewind(fp);
    fscanf(fp, "%d %d %lf", &temp, &method, &temp1);
    if(mpirank == 0) 
      printf(" PoleNum: %d method: %d temperature: %lf\n", temp, method, temp1);
    fflush(stdout);
    options.numPole  = temp;
    options.temperature  = temp1;
    options.method = method;
  
    fclose(fp);
    */

    pexsi.LoadComplexMatrix(
        HMat.size,                        
        HMat.nnz,                          
        HMat.nnzLocal,                     
        numColLocal,                  
        HMat.colptrLocal.Data(),
        HMat.rowindLocal.Data(),
        HMat.nzvalLocal.Data(),
        isSIdentity,
        SMat.nzvalLocal.Data(),
        options.solver,
        options.verbosity );

    if( mpirank == 0 ){
      printf("After loading the matrix.\n");
    }

    if( options.matrixType == 0 ){

      // No permutation
      pexsi.SymbolicFactorizeComplexSymmetricMatrix( 
          options.solver,
          options.symmetricStorage,
          colPerm, 
          options.npSymbFact,
          options.verbosity );

      // Can have permutation
      pexsi.SymbolicFactorizeComplexUnsymmetricMatrix( 
          options.solver,
          colPerm, 
          rowPerm,
          options.npSymbFact,
          options.transpose,
          reinterpret_cast<double*>(HMat.nzvalLocal.Data()),
          options.verbosity );
    }

    // Inertia counting
    Int numShift = options.numPole;
    std::vector<double> shiftVec( numShift );
    for( Int i = 0; i < numShift; i++ ){
      shiftVec[i] = options.muMin0 + 
        ( options.muMax0 - options.muMin0 ) / numShift * double(i);
    }
    std::vector<double> inertiaVec( numShift );
    pexsi.CalculateNegativeInertiaComplex(
        shiftVec,
        inertiaVec,
        options.solver,
        options.verbosity );

    if( mpirank == 0 ){
      for( Int i = 0; i < numShift; i++ )
        printf( "Shift = %25.15f  inertia = %25.1f\n", 
            shiftVec[i], inertiaVec[i] );
    }


    // Compute Fermi operator
    pexsi.CalculateFermiOperatorComplex(
        options.numPole,
        options.temperature,
        options.gap,
        options.deltaE,
        options.mu0,
        numElectronExact, 
        options.numElectronPEXSITolerance,
        options.solver,
        options.verbosity,
        numElectron,
        numElectronDrvMu,
        options.method,
        options.nPoints,
        options.spin );

    if( mpirank == 0 ){
      printf("numElectron       = %25.15f\n", numElectron);
      printf("numElectronDrvMu  = %25.15f\n", numElectronDrvMu);
    }

  }
  catch( std::exception& e )
  {
    statusOFS << std::endl << " ERROR!!! Proc " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
    statusOFS.close();
    statusOFS << std::endl << " ERROR!!! Proc " << mpirank << " caught exception with message: "
      << e.what() << std::endl;
  }

  MPI_Finalize();

  return 0;
}
