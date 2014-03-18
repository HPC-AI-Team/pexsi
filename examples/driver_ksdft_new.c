/*
   Copyright (c) 2012 The Regents of the University of California,
   through Lawrence Berkeley National Laboratory.  

   Authors: Lin Lin

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
/**
 * @file driver_ksdft_new.c
 * @brief Example for using the new driver interface for performing KSDFT
 * calculations.
 *
 * This file is eventually going to be merged with the driver_ksdft.c
 *
 *
 * @date 2014-03-07
 */
#include  <stdio.h>
#include  <stdlib.h>
#include  "c_pexsi_new_interface.h"

int main(int argc, char **argv) 
{
  int mpirank, mpisize;
  int           nrows;
  int           nnz;                          
  int           nnzLocal;                     
  int           numColLocal;                  
  int*          colptrLocal;                  
  int*          rowindLocal;                  
  double*       HnzvalLocal;                  
  int           isSIdentity;                  
  double*       SnzvalLocal;                  
  double*       AinvnzvalLocal;

  double*       DMnzvalLocal;
  double*       EDMnzvalLocal;
  double*       FDMnzvalLocal;
  double*       localDOSnzvalLocal;  

  double        Energy;
  double        eta;

  double        numElectronExact;

  double        muPEXSI;
  double        numElectronPEXSI;
  double        muMinInertia;
  double        muMaxInertia;
  int           numTotalInertiaIter;
  int           numTotalPEXSIIter;

  double        totalEnergyH;
  double        totalEnergyS;
  double        totalFreeEnergy;
  
  char*         Hfile;
  char*         Sfile;
  int           isFormatted;


  int           i, j;
  int           nprow, npcol;
  MPI_Comm      readComm;
  int           isProcRead;
  int           info;


  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

  /* Below is the data used for the toy g20 matrix */

  numElectronExact    = 12.0;
  nprow               = 1;
  npcol               = 1;
  Hfile               = "lap2dr.matrix";
  Sfile               = "";
  isFormatted         = 1;
  isSIdentity         = 1;
  Energy              = 1.0;
  eta                 = 0.001;

  /* Split the processors to read matrix */
  if( mpirank < nprow * npcol )
    isProcRead = 1;
  else
    isProcRead = 0;

  MPI_Comm_split( MPI_COMM_WORLD, isProcRead, mpirank, &readComm );

  if( isProcRead == 1 ){
    printf("Proc %5d is reading file...\n", mpirank );
    /* Read the matrix head for allocating memory */
    if( isFormatted == 1 ){
      ReadDistSparseMatrixFormattedHeadInterface(
          Hfile,
          &nrows,
          &nnz,
          &nnzLocal,
          &numColLocal,
          readComm );
    }
    else{
      ReadDistSparseMatrixHeadInterface(
          Hfile,
          &nrows,
          &nnz,
          &nnzLocal,
          &numColLocal,
          readComm );
    }
    
    if( mpirank == 0 ){
      printf("On processor 0...\n");
      printf("nrows       = %d\n", nrows );
      printf("nnz         = %d\n", nnz );
      printf("nnzLocal    = %d\n", nnzLocal );
      printf("numColLocal = %d\n", numColLocal );
    }


    /* Allocate memory visible to processors in the group of readComm */
    colptrLocal             = (int*) malloc( sizeof(int) * (numColLocal+1) );
    rowindLocal             = (int*) malloc( sizeof(int) * nnzLocal );
    HnzvalLocal             = (double*) malloc( sizeof(double) * nnzLocal );
    SnzvalLocal             = (double*) malloc( sizeof(double) * nnzLocal );
    DMnzvalLocal            = (double*) malloc( sizeof(double) * nnzLocal );
    EDMnzvalLocal           = (double*) malloc( sizeof(double) * nnzLocal );
    FDMnzvalLocal           = (double*) malloc( sizeof(double) * nnzLocal );
    localDOSnzvalLocal      = (double*) malloc( sizeof(double) * nnzLocal );

    /* Actually read the matrix */
    if( isFormatted == 1 ){
      ReadDistSparseMatrixFormattedInterface(
          Hfile,
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          HnzvalLocal,
          readComm );
    }
    else{
      ParaReadDistSparseMatrixInterface(
          Hfile,
          nrows,
          nnz,
          nnzLocal,
          numColLocal,
          colptrLocal,
          rowindLocal,
          HnzvalLocal,
          readComm );
    }

    if( isSIdentity == 0 ){
      if( isFormatted == 1 ){
        ReadDistSparseMatrixFormattedInterface(
            Hfile,
            nrows,
            nnz,
            nnzLocal,
            numColLocal,
            colptrLocal,
            rowindLocal,
            SnzvalLocal,
            readComm );
      }
      else{
        ParaReadDistSparseMatrixInterface(
            Hfile,
            nrows,
            nnz,
            nnzLocal,
            numColLocal,
            colptrLocal,
            rowindLocal,
            SnzvalLocal,
            readComm );
      }
    }
    
    if( mpirank == 0 ){ 
      printf("Finish reading the matrix.\n");
    }
  }

  /* Call the PEXSI interface */

  /* Step 1. Initialize PEXSI */

  PPEXSIOptions  options;
  PPEXSISetDefaultOptions( &options );
  options.muMin0 = 0.0;
  options.muMax0 = 0.5;
  options.npSymbFact = 1;
  options.ordering = 0;
  options.isInertiaCount = 1;
  options.verbosity = 1;
  options.deltaE   = 20.0;
  options.numPole  = 60;
  options.temperature  = 0.019; // 3000K
  options.muPEXSISafeGuard  = 0.2; 
  options.numElectronPEXSITolerance = 0.0001;

  PPEXSIPlan   plan;

  plan = PPEXSIPlanInitialize( 
      MPI_COMM_WORLD, 
      nprow,
      npcol,
      mpirank, 
      &info );

  PPEXSILoadRealSymmetricHSMatrix( 
      plan, 
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      isSIdentity,
      SnzvalLocal,
      &info );

  /* Step 2. PEXSI Solve */

  PPEXSIDFTDriver(
      plan,
      numElectronExact,
      options,
      &muPEXSI,                   
      &numElectronPEXSI,         
      &muMinInertia,              
      &muMaxInertia,             
      &numTotalInertiaIter,   
      &numTotalPEXSIIter,   
      &info );

  if( isProcRead == 1 ){
    PPEXSIRetrieveRealSymmetricDFTMatrix(
        plan,
        DMnzvalLocal,
        EDMnzvalLocal,
        FDMnzvalLocal,
        &totalEnergyH,
        &totalEnergyS,
        &totalFreeEnergy,
        &info );

    if( mpirank == 0 ){
      printf("Output from the main program\n");
      printf("Total energy (H*DM)         = %15.5f\n", totalEnergyH);
      printf("Total energy (S*EDM)        = %15.5f\n", totalEnergyS);
      printf("Total free energy           = %15.5f\n", totalFreeEnergy);
    }
  }

//  if( info != 0 ){
//    if( mpirank == 0 ){
//      printf("PEXSI solve routine gives info = %d. Exit now.\n", info );
//    }
//    MPI_Finalize();
//    return info;
//  }
//
//  if( mpirank == 0 ){ 
//    printf("PEXSI Solve finished. \n");
//  }


  /* Step 3. Post processing */
  
  /* Compute the density of states (DOS) via inertia counting (without
   * including finite temperature effects) */

  /* Step 4. Clean up */

  PPEXSIPlanFinalize( 
      plan,
      &info );
  
  if( mpirank == 0 ){ 
    printf("\nAll calculation is finished. Exit the program.\n");
  }


  if( isProcRead == 1 ){
    free( colptrLocal );
    free( rowindLocal );
    free( HnzvalLocal );
    free( SnzvalLocal );
    free( DMnzvalLocal );
    free( EDMnzvalLocal );
    free( FDMnzvalLocal );
    free( localDOSnzvalLocal );
  }

  
  MPI_Comm_free( &readComm );
  MPI_Finalize();

  return 0;
}
