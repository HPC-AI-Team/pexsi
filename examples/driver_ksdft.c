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
  Hfile               = "g20.matrix";
  Sfile               = "";
  isFormatted         = 1;
  isSIdentity         = 1;
  ordering            = 2;

  /* Split the processors to read matrix */
  if( mpirank < npPerPole ) 
    isProcRead = 1;
  else
    isProcRead = 0;

  MPI_Comm_split( MPI_COMM_WORLD, isProcRead, mpirank, &readComm );

  /* Allocate memory */
  muList                  = (double*) malloc( sizeof(double) * muMaxIter );
  numElectronList         = (double*) malloc( sizeof(double) * muMaxIter );
  numElectronDrvList      = (double*) malloc( sizeof(double) * muMaxIter );
  shiftList               = (double*) malloc( sizeof(double) * numPole   );
  inertiaList             = (double*) malloc( sizeof(double) * numPole   );

  if( isProcRead == 1 ){

  }

//  /* Read the matrix */
//  ReadDistSparseMatrixFormattedHeadInterface(
//      Hfile,
//      &nrows,
//      &nnz,
//      &nnzLocal,
//      &numColLocal,
//      MPI_COMM_WORLD );
//      
//  if( mpirank == 0 ){
//    printf("On processor 0...\n");
//    printf("nrows       = %d\n", nrows );
//    printf("nnz         = %d\n", nnz );
//    printf("nnzLocal    = %d\n", nnzLocal );
//    printf("numColLocal = %d\n", numColLocal );
//  }
//
//  /* Allocate memory */
//  colptrLocal = (int*)malloc( (numColLocal+1) * sizeof(int) );
//  rowindLocal = (int*)malloc( nnzLocal * sizeof(int) );
//  HnzvalLocal = (double*)malloc( nnzLocal * sizeof(double) );
//  AinvnzvalLocal = (double*)malloc( 2*nnzLocal * sizeof(double) );
//
//  /* Read the matrix */
//  ReadDistSparseMatrixFormattedInterface(
//      Hfile,
//      nrows,
//      nnz,
//      nnzLocal,
//      numColLocal,
//      colptrLocal,
//      rowindLocal,
//      HnzvalLocal,
//      MPI_COMM_WORLD );
//
//
//  /* Other parameters */
//  isSIdentity = 1;
//  ordering    = 0;
//  zShift[0]   = 0.5;
//  zShift[1]   = 0.0;
//  npSymbFact  = 1;
//
//  PPEXSISelInvInterface(
//      nrows,
//      nnz,
//      nnzLocal,
//      numColLocal,
//      colptrLocal,
//      rowindLocal,
//      HnzvalLocal,
//      isSIdentity,
//      SnzvalLocal,
//      zShift,
//      ordering,
//      npSymbFact,
//      MPI_COMM_WORLD,
//      AinvnzvalLocal,
//      &info );
//
//  /* The first processor output the diagonal elements in natural order
//   */
//  if( mpirank == 0 ){
//    numColLocalFirst = nrows / mpisize;
//    firstCol         = mpirank * numColLocalFirst;
//    for( j = 0; j < numColLocal; j++ ){
//      jcol = firstCol + j + 1;
//      for( i = colptrLocal[j]-1; 
//           i < colptrLocal[j+1]-1; i++ ){
//        irow = rowindLocal[i];
//        if( irow == jcol ){
//          printf("Ainv[%5d,%5d] = %15.10e + %15.10e i\n", 
//              irow, irow,
//              AinvnzvalLocal[2*i],
//              AinvnzvalLocal[2*i+1]);
//        }
//      }
//    } // for (j)
//  }
  

  /* Deallocate memory */
  free( muList );
  free( numElectronList );
  free( numElectronDrvList );
  free( shiftList );
  free( inertiaList );


//  free( colptrLocal );
//  free( rowindLocal );
//  free( HnzvalLocal );
//  free( AinvnzvalLocal );

  
  MPI_Finalize();

  return 0;
}
