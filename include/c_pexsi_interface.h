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

/**
 * @brief Dummy interface for test purpose.
 *
 */
void DummyInterface( MPI_Comm comm, int a );

/**
 * @brief Main C-interface with PEXSI.
 *
 */

#ifdef __cplusplus
}// extern "C"
#endif

#endif // _C_PEXSI_INTERFACE_H_

