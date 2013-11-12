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
 * @brief Interface subroutines of %PEXSI that can be called by C.
 *
 * The FORTRAN interface are wrappers of the C interface.  The header
 * file is not needed here, and the interface routines are given
 * directly in interface.cpp.
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


/**
 * @brief Read the sizes of a DistSparseMatrix in formatted form (txt)
 * for allocating memory in C.
 *
 * @param[in] filename (global) Filename for the input matrix.
 * @param[out] size (global) Number of rows and columns of the matrix.
 * @param[out] nnz (global) Total number of nonzeros.
 * @param[out] nnzLocal (local) Number of local nonzeros.
 * @param[out] numColLocal (local) Number of local columns.
 * @param[in]  comm (global) MPI communicator.
 */
void ReadDistSparseMatrixFormattedHeadInterface ( 
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    MPI_Comm comm );

/**
 * @brief Reading the data of a formatted DistSparseMatrix. 
 *
 * This routine assumes that the arrays have been allocated outside this
 * subroutine.
 *
 * @param[in] filename (global) Filename for the input matrix.
 * @param[in] size (global) Number of rows and columns of the matrix.
 * @param[in] nnz (global) Total number of nonzeros.
 * @param[in] nnzLocal (local) Number of local nonzeros.
 * @param[in] numColLocal (local) Number of local columns.
 * @param[out] colptrLocal (local) Dimension: numColLocal+1. Local column
 * pointer in CSC format.
 * @param[out] rowindLocal (local) Dimension: nnzLocal. Local row index
 * pointer in CSC format.
 * @param[out] nzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values in CSC format.
 * @param[in]  comm (global) MPI communicator.
 */
void ReadDistSparseMatrixFormattedInterface(
		char*     filename,
		int       size,
		int       nnz,
		int       nnzLocal,
		int       numColLocal,
		int*      colptrLocal,
		int*      rowindLocal,
		double*   nzvalLocal,
		MPI_Comm  comm );


/**
 * @brief Read the sizes of a DistSparseMatrix in unformatted form
 * (csc) for allocating memory in C.
 *
 * @param[in] filename (global) Filename for the input matrix.
 * @param[out] size (global) Number of rows and columns of the matrix.
 * @param[out] nnz (global) Total number of nonzeros.
 * @param[out] nnzLocal (local) Number of local nonzeros.
 * @param[out] numColLocal (local) Number of local columns.
 * @param[in]  comm (global) MPI communicator.
 */
void ReadDistSparseMatrixHeadInterface ( 
    char*    filename,
    int*     size,
    int*     nnz,
    int*     nnzLocal,
    int*     numColLocal,
    MPI_Comm comm );

/**
 * @brief Actual reading the data of a DistSparseMatrix using MPI-IO,
 * assuming that the arrays have been allocated outside this
 * subroutine.
 *
 * This routine can be much faster than reading a DistSparseMatrix
 * sequentially, especially compared to the version using formatted
 * input @ref ReadDistSparseMatrixFormattedInterface.
 *
 * @param[in] filename (global) Filename for the input matrix.
 * @param[in] size (global) Number of rows and columns of the matrix.
 * @param[in] nnz (global) Total number of nonzeros.
 * @param[in] nnzLocal (local) Number of local nonzeros.
 * @param[in] numColLocal (local) Number of local columns.
 * @param[out] colptrLocal (local) Dimension: numColLocal+1. Local column
 * pointer in CSC format.
 * @param[out] rowindLocal (local) Dimension: nnzLocal. Local row index
 * pointer in CSC format.
 * @param[out] nzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values in CSC format.
 * @param[in]  comm (global) MPI communicator.
 */
void ParaReadDistSparseMatrixInterface ( 
    char*     filename,
    int       size,
    int       nnz,
    int       nnzLocal,
    int       numColLocal,
    int*      colptrLocal,
    int*      rowindLocal,
    double*   nzvalLocal,
    MPI_Comm  comm );



/**
 * @brief Driver interface for computing the selected  elements of a matrix.  
 *
 * Evaluate the selected elements of 
 * \f$(H - z S)^{-1}\f$ for a shift \f$z\in\mathbb{C}\f$.
 * 
 * **Note** 
 * - H and S are always real symmetric matrices, are saved in the
 *   compressed sparse column (CSC) format, and have the same sparsity
 *   pattern.
 * - Complex arithmetic operation is used for the factorization
 *   and selected inversion even if z is real.
 * - The input of S is optional, and the shift z can be 0.
 *
 * @param[in] nrows (global) Number of rows and columns of the matrix.
 * @param[in] nnz (global) Total number of nonzeros of H.
 * @param[in] nnzLocal (local) Number of local nonzeros of H.
 * @param[in] numColLocal (local) Number of local columns for H.
 * @param[in] colptrLocal (local) Dimension: numColLocal+1. Local column
 * pointer in CSC format.
 * @param[in] rowindLocal (local) Dimension: nnzLocal. Local row index
 * pointer in CSC format.
 * @param[in] HnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * values of H in CSC format.
 * @param[in] isSIdentity (global) Whether S is an identity matrix. 
 * If so, the variable SnzvalLocal is omitted.
 * @param[in] SnzvalLocal (local) Dimension: nnzLocal. Local nonzero
 * value of S in CSC format.
 * @param[in] zShift  (global) Dimension: 2. Use 2 real numbers to
 * denote a complex shift. 
 * - Real part: zShift[0]
 * - Imag part: zShift[1]    
 * @param[in] ordering (global) Ordering strategy for factorization and selected
 * inversion.  
 * - = 0   : Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS
 *   option in SuperLU_DIST).
 * - = 1   : Sequential ordering using METIS (METIS_AT_PLUS_A
 *   option in SuperLU_DIST).
 * - = 2   : Multiple minimum degree ordering (MMD_AT_PLUS_A
 *   option in SuperLU_DIST).
 * @param[in] npSymbFact (global) Number of processors for PARMETIS/PT-SCOTCH.  Only used if the ordering = 0.
 * @param[in] comm (global) MPI communicator.
 * @param[out] AinvnzvalLocal (local) Dimension: 2*nnzLocal. Local nonzero
 * values of the selected elements of of \f$(H - z S)^{-1}\f$.
 * - Use 2 double for one complex number. This ensures the compatibility with FORTRAN.
 * - Real part: AinvnzvalLocal[2*k]. Imag part: AinvnzvalLocal[2*k+1].
 * @param[out] info 
 * - = 0: successful exit.  
 * - > 0: unsuccessful.
 */
void PPEXSISelInvInterface ( 
    int           nrows,                        
    int           nnz,                          
    int           nnzLocal,                     
    int           numColLocal,                  
    int*          colptrLocal,                  
    int*          rowindLocal,                  
    double*       HnzvalLocal,                  
    int           isSIdentity,                  
    double*       SnzvalLocal,                  
    double*       zShift,                       
    int           ordering,                
    int           npSymbFact,                   
    MPI_Comm	    comm,                         
    double*       AinvnzvalLocal,
    int*          info
    );




#ifdef __cplusplus
}// extern "C"
#endif

#endif // _C_PEXSI_INTERFACE_H_

