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
 * @file driver_pselinv.c
 * @brief Example for using the driver interface for parallel selected
 * inversion of a complex symmetric matrix.
 *
 * @see PPEXSISelInvInterface 
 * @date 2013-11-10
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
  double        zShift[2];                       
  int           ordering;                
  int           npSymbFact;                   
  double*       AinvnzvalLocal;
  int           info;
  char*         Hfile = "g20.matrix";

  int           i, j, irow, jcol;
  int           numColLocalFirst, firstCol;

  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

  /* Read the matrix */
  ReadDistSparseMatrixFormattedHeadInterface(
      Hfile,
      &nrows,
      &nnz,
      &nnzLocal,
      &numColLocal,
      MPI_COMM_WORLD );
      
  if( mpirank == 0 ){
    printf("On processor 0...\n");
    printf("nrows       = %d\n", nrows );
    printf("nnz         = %d\n", nnz );
    printf("nnzLocal    = %d\n", nnzLocal );
    printf("numColLocal = %d\n", numColLocal );
  }

  /* Allocate memory */
  colptrLocal = (int*)malloc( (numColLocal+1) * sizeof(int) );
  rowindLocal = (int*)malloc( nnzLocal * sizeof(int) );
  HnzvalLocal = (double*)malloc( nnzLocal * sizeof(double) );
  AinvnzvalLocal = (double*)malloc( 2*nnzLocal * sizeof(double) );

  /* Read the matrix */
  ReadDistSparseMatrixFormattedInterface(
      Hfile,
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      MPI_COMM_WORLD );


  /* Other parameters */
  isSIdentity = 1;
  ordering    = 0;
  zShift[0]   = 0.5;
  zShift[1]   = 0.0;
  npSymbFact  = 1;

  PPEXSISelInvInterface(
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      HnzvalLocal,
      isSIdentity,
      SnzvalLocal,
      zShift,
      ordering,
      npSymbFact,
      MPI_COMM_WORLD,
      AinvnzvalLocal,
      &info );

  /* The first processor output the diagonal elements in natural order
   */
  if( mpirank == 0 ){
    numColLocalFirst = nrows / mpisize;
    firstCol         = mpirank * numColLocalFirst;
    for( j = 0; j < numColLocal; j++ ){
      jcol = firstCol + j + 1;
      for( i = colptrLocal[j]-1; 
           i < colptrLocal[j+1]-1; i++ ){
        irow = rowindLocal[i];
        if( irow == jcol ){
          printf("Ainv[%5d,%5d] = %15.10e + %15.10e i\n", 
              irow, irow,
              AinvnzvalLocal[2*i],
              AinvnzvalLocal[2*i+1]);
        }
      }
    } // for (j)
  }
  

  /* Deallocate memory */
  free( colptrLocal );
  free( rowindLocal );
  free( HnzvalLocal );
  free( AinvnzvalLocal );

  
  MPI_Finalize();

  return 0;
}
