/**
 * @file c_pexsi_interface.h
 * @brief Interface subroutines of PEXSI that can be called by C/FORTRAN
 * @author Lin Lin
 * @version 
 * @date 2013-01-31
 */
#ifndef _C_PEXSI_INTERFACE_H_ 
#define _C_PEXSI_INTERFACE_H_
#include <mpi.h>

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
void PEXSIInterface( int*          nrows, 
										 int*          nnzLocal,
										 int*          numColLocal,
										 int*          colptrLocal,
										 int*          rowindLocal,
										 double*       HnzvalLocal,
										 double*       SnzvalLocal,
										 double*      DMnzvalLocal,
										 double*     EDMnzvalLocal,
										 double*     FDMnzvalLocal,
										 int*          numPole,
										 double*       temperature,
										 double*       numElectronExact,
										 double*       gap,
										 double*       deltaE,
										 double*       mu0,
										 double*       muMin,
										 double*       muMax,
										 int*          muMaxIter,
										 double*       poleTolerance,
										 double*       numElectronTolerance,
										 int*          isConverged );

#ifdef __cplusplus
}// extern "C"
#endif

#endif // _C_PEXSI_INTERFACE_H_

