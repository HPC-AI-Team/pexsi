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
 * @file driver_ksdft.c
 * @brief Example for using the driver interface for performing KSDFT
 * calculations.
 *
 * See "Tutorial" in the documentation for explanation of the
 * parameters.
 *
 * @see f_ppexsi.f90
 * @date 2013-12-14
 */
#include  <stdio.h>
#include  <stdlib.h>
#include  "c_pexsi_interface.h"

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
  int           ordering;                
  int           npSymbFact;                   
  double*       AinvnzvalLocal;

  int*          inertiaListInt;
  double*       DMnzvalLocal;
  double*       EDMnzvalLocal;
  double*       FDMnzvalLocal;
  double*       muList;
  double*       numElectronList;
  double*       numElectronDrvList;
  double*       shiftList;
  double*       inertiaList;
  double*       localDOSnzvalLocal;  

  double        Energy;
  double        eta;

  int           numPole;
  int           npPerPole;

  double        temperature;
  double        numElectronExact;
  double        numElectron;
  double        gap;
  double        deltaE;
  double        muMin0;
  double        muMax0;
  double        muInertia;
  double        muMinInertia;
  double        muMaxInertia;
  double        muLowerEdge;
  double        muUpperEdge;
  double        muPEXSI;
  double        muMinPEXSI;
  double        muMaxPEXSI;

  int           inertiaMaxIter;
  int           inertiaIter;
  int           muMaxIter;
  int           muIter;
  int           isInertiaCount;
  
  char*         Hfile;
  char*         Sfile;
  int           isFormatted;

  double        PEXSINumElectronTolerance;
  double        inertiaNumElectronTolerance;


  int           i, j, irow, jcol;
  int           numColLocalFirst, firstCol;
  MPI_Comm      readComm;
  int           isProcRead;
  int           info;


  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

  /* Below is the data used for the toy g20 matrix */

  temperature         = 0.0019;
  numElectronExact    = 12.0;
  numPole             = 40;
  gap                 = 0.0;
  deltaE              = 30.0;
  muMin0              = 0.0;
  muMax0              = 10.0;
  inertiaMaxIter      = 5;
  muMaxIter           = 5;
  inertiaNumElectronTolerance = 10.0;
  PEXSINumElectronTolerance   = 0.01;
  npPerPole           = 1;
  npSymbFact          = 1;
  Hfile               = "lap2dr.matrix";
  Sfile               = "";
  isFormatted         = 1;
  isSIdentity         = 1;
  ordering            = 1;
  Energy              = 1.0;
  eta                 = 0.001;

  /* Split the processors to read matrix */
  if( mpirank < npPerPole )
    isProcRead = 1;
  else
    isProcRead = 0;

  MPI_Comm_split( MPI_COMM_WORLD, isProcRead, mpirank, &readComm );

  /* Allocate memory visible to all processors */
  muList                  = (double*) malloc( sizeof(double) * muMaxIter );
  numElectronList         = (double*) malloc( sizeof(double) * muMaxIter );
  numElectronDrvList      = (double*) malloc( sizeof(double) * muMaxIter );
  shiftList               = (double*) malloc( sizeof(double) * numPole   );
  inertiaList             = (double*) malloc( sizeof(double) * numPole   );
  inertiaListInt          = (int*) malloc( sizeof(int) * numPole   );


  if( isProcRead == 1 ){
    printf("Proc %5d is reading file...", mpirank );
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

  /* Step 1. Estimate the range of chemical potential */
  
  PPEXSIInertiaCountInterface(
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      isSIdentity,
      SnzvalLocal,
      temperature,
      numElectronExact,
      muMin0,
      muMax0,
      numPole,
      inertiaMaxIter,
      inertiaNumElectronTolerance,
      ordering,
      npPerPole,
      npSymbFact,
      MPI_COMM_WORLD,
      &muMinInertia,
      &muMaxInertia,
      &muLowerEdge,
      &muUpperEdge,
      &inertiaIter,
      shiftList,
      inertiaList,
      &info);

  if( info != 0 ){
    if( mpirank == 0 ){
      printf("Inertia count routine gives info = %d. Exit now.\n", info );
    }
    MPI_Finalize();
    return info;
  }

  muInertia = (muLowerEdge + muUpperEdge)/2.0;

  if( mpirank == 0 ){ 
    printf("The computed finite temperature inertia = \n");
    for( i = 0; i < numPole; i++ )
      printf( "Shift = %25.15f, inertia = %25.15f\n", 
          shiftList[i], inertiaList[i] );
  }
  
  /* Step 2. Solve KSDFT using PEXSI */
  PPEXSISolveInterface(
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      isSIdentity,
      SnzvalLocal,
      temperature,
      numElectronExact,
      muInertia,
      muMinInertia,
      muMaxInertia,
      gap,
      deltaE,
      numPole,
      muMaxIter,
      PEXSINumElectronTolerance,
      ordering,
      npPerPole,
      npSymbFact,
      MPI_COMM_WORLD,
      DMnzvalLocal,
      EDMnzvalLocal,
      FDMnzvalLocal,
      &muPEXSI,
      &numElectron,
      &muMinPEXSI,
      &muMaxPEXSI,
      &muIter,
      muList,
      numElectronList,
      numElectronDrvList,
      &info );

  if( info != 0 ){
    if( mpirank == 0 ){
      printf("PEXSI solve routine gives info = %d. Exit now.\n", info );
    }
    MPI_Finalize();
    return info;
  }

  if( mpirank == 0 ){ 
    printf("PEXSI Solve finished. \n");
  }


  /* Step 3. Post processing */
  
  /* Compute the density of states (DOS) via inertia counting (without
   * including finite temperature effects) */

  PPEXSIRawInertiaCountInterface(
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      isSIdentity,
      SnzvalLocal,
      muMinInertia,
      muMaxInertia,
      numPole,
      ordering,
      npPerPole,
      npSymbFact,
      MPI_COMM_WORLD,
      shiftList,
      inertiaListInt,
      &info);

  if( info != 0 ){
    if( mpirank == 0 ){
      printf("PEXSI raw inertia counting routine gives info = %d. Exit now.\n", info );
    }
    MPI_Finalize();
    return info;
  }


  if( mpirank == 0 ){ 
    printf("PEXSI raw inertia counting routine finished. \n");
    printf("The computed raw inertia (zero temperature) = \n");
    for( i = 0; i < numPole; i++ )
      printf( "Shift = %25.15f, inertia = %15d\n", 
          shiftList[i], inertiaListInt[i] );
  }


  /* Compute the local density of states (LDOS).
   *
   * Only the first pole group participates in the computation of the
   * selected inversion for a single shift. */

  if( isProcRead == 1 ){ 
    PPEXSILocalDOSInterface(
        nrows,
        nnz,
        nnzLocal,
        numColLocal,
        colptrLocal,
        rowindLocal,
        HnzvalLocal,
        isSIdentity,
        SnzvalLocal,
        Energy,
        eta,
        ordering,
        npSymbFact,
        readComm,
        localDOSnzvalLocal,
        &info);

    if( info != 0 ){
      if( mpirank == 0 ){
        printf("PEXSI LDOS calculation gives info = %d. Exit now.\n", info );
      }
      MPI_Finalize();
      return info;
    }
  }

  if( mpirank == 0 ){ 
    printf("Local DOS calculation has finished.\n");
  }

  if( mpirank == 0 ){ 
    printf("\nAll calculation is finished. Exit the program.\n");
  }


  /* Deallocate memory */
  free( muList );
  free( numElectronList );
  free( numElectronDrvList );
  free( shiftList );
  free( inertiaList );
  free( inertiaListInt );

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
