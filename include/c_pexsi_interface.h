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
/**
 * @file c_pexsi_interface.h
 * @brief Interface subroutines of PEXSI that can be called by C/FORTRAN
 * 
 * This file is to be updated in the end when the interface is fixed,
 * and when C linkage with PEXSI is needed.
 *
 * This file is NOT needed for FORTRAN interface. All subroutines for
 * the FORTRAN interface is contained in src/interface.cpp. 
 *
 * @date 2013-01-31
 */
#ifndef _C_PEXSI_INTERFACE_H_ 
#define _C_PEXSI_INTERFACE_H_
#include "mpi.h"

// FIXME This file is to be updated in the end when the interface is fixed, and when C linkage with PEXSI is needed.

#ifdef __cplusplus
extern "C"{
#endif

/// @brief Read the sizes of a DistSparseMatrix in unformatted form
/// (csc) for allocating memory in C.
void ReadDistSparseMatrixHeadInterface (
                                        char*    filename,
                                        int*     size,
                                        int*     nnz,
                                        int*     nnzLocal,
                                        int*     numColLocal,
                                        MPI_Comm comm );

/// @brief Actual reading the data of a DistSparseMatrix using MPI-IO,
/// assuming that the arrays have been allocated outside this
/// subroutine.
///
/// This routine can be much faster than the sequential reading of a
/// DistSparseMatrix in the presence of a large number of processors.
void ParaReadDistSparseMatrixInterface(
                                       char*     filename,
                                       int       size,
                                       int       nnz,
                                       int       nnzLocal,
                                       int       numColLocal,
                                       int*      colptrLocal,
                                       int*      rowindLocal,
                                       double*   nzvalLocal,
                                       MPI_Comm  comm );


#ifdef __cplusplus
}// extern "C"
#endif

#endif // _C_PEXSI_INTERFACE_H_

