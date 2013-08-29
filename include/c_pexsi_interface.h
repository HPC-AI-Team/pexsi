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
 * @author Lin Lin
 * @date 2013-01-31
 */
#ifndef _C_PEXSI_INTERFACE_H_ 
#define _C_PEXSI_INTERFACE_H_
#include <mpi.h>

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

