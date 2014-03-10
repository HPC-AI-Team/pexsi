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
 * @file c_pexsi_new_interface.h
 * @brief New interface subroutines of %PEXSI that can be called by C.
 * 
 * The existence of this file is temporary.  The new interface will
 * eventually be merged together with c_pexsi_interface.h
 * 
 * @date 2014-03-07
 */
#ifndef _C_PEXSI_NEW_INTERFACE_H_ 
#define _C_PEXSI_NEW_INTERFACE_H_
#include "mpi.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C"{
#endif

// *********************************************************************
// The following routines belong to the second version of the interface
// *********************************************************************

/**
 * @brief A handle for holding the internal %PEXSI data structure.  
 *
 * @note This handle can be also used with FORTRAN, with the `INTEGER*8`
 * data structure, or the `INTEGER(C_INTPTR_T)` structure if
 * ISO_C_BINDING is used.
 */
typedef intptr_t  PPEXSIPlan;

/**
 * @struct PPEXSIDFTOptions
 * @brief Structure for the input parameters in DFT calculations.
 */
typedef struct {
    /** 
     * @brief  Temperature, in the same unit as H 
     */ 
    double        temperature;  
    /** 
     * @brief  An upper bound for the spectral radius of \f$S^{-1} H\f$.
     */ 
    double        deltaE;
    /** 
     * @brief  Number of terms in the pole expansion.
     */ 
    int           numPole;
    /** 
     * @brief  Whether inertia counting is used at the very beginning.
     */ 
    int           isInertiaCount;
    /** 
     * @brief  Maximum number of PEXSI iterations after each inertia
     * counting procedure.
     */ 
    int           maxPEXSIIter;
    /** 
     * @brief  Initial guess of lower bound for mu.
     */ 
    double        muMin0;
    /** 
     * @brief  Initial guess of upper bound for mu.
     */ 
    double        muMax0;
    /** 
     * @brief  Stopping criterion in terms of the chemical potential
     * for the inertia counting procedure.
     */ 
    double        muInertiaTolerance;
    /** 
     * @brief  Safe guard criterion in terms of the chemical potential
     * to reinvoke the inertia counting procedure.
     */ 
    double        muPEXSISafeGuard;
    /** 
     * @brief  Stopping criterion of the %PEXSI iteration in terms of the
     * number of electrons compared to numElectronExact.
     */ 
    double        numElectronPEXSITolerance;
    /** 
     * @brief  Ordering strategy for factorization and selected
     * inversion.  
     * - = 0   : Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS
     *   option in SuperLU_DIST).
     * - = 1   : Sequential ordering using METIS (METIS_AT_PLUS_A
     *   option in SuperLU_DIST).
     * - = 2   : Multiple minimum degree ordering (MMD_AT_PLUS_A
     *   option in SuperLU_DIST).
     */ 
    int           ordering;
    /** 
     * @brief  Number of processors for PARMETIS/PT-SCOTCH.  Only used
     * if the ordering = 0.
     */ 
    int           npSymbFact;
} PPEXSIDFTOptions;


/**
 * @brief Set the default options for DFT driver.
 *
 * All default values assume the input unit (for H) is Rydberg.
 *
 * @param[in] options (global) Pointer to the options containing input
 * parameters for the driver.  
 */
void PPEXSISetDefaultDFTOptions(
    PPEXSIDFTOptions*   options );


/**
 * @brief Initialize the %PEXSI plan.
 *
 * @todo A specicial FORTRAN interface for this routine.  This is
 * because the comm between C and FORTRAN are different.  All subsequent
 * interface routines do not require anymore interface in this routine,
 * but should be taken care of using ISO_C_BINDING.
 *
 * @param[in] comm  (global) Communicator used for the entire %PEXSI procedure.  The
 * size of this communicator should be a multiple of npPerPole =
 * numProcRow * numProcCol, which is the total number of processors used
 * by each individual pole.
 * @param[in] numProcRow (global) Number of processors in the row communication group
 * for each pole.
 * @param[in] numProcCol (global) Number of processors in the column communication group
 * for each pole.
 * @param[in] outputFileIndex (local) The index for the %PEXSI output file.  For
 * instance, if this index is 1, then the corresponding processor will
 * output to the file `logPEXSI1`.  
 *
 * @note Each processor must output to a **different** file.  By
 * default, outputFileIndex can be set as mpirank.
 *
 *
 * @return (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 */
PPEXSIPlan PPEXSIPlanInitialize(
    MPI_Comm      comm,
    int           numProcRow, 
    int           numProcCol,
    int           outputFileIndex );


/**
 * @brief Load the H and S matrix into the %PEXSI internal data
 * structure.
 *
 * All subsequent operations require this routine to be called first
 *
 * @note Only input from the processors associated with the first pole
 * is required. The information will be broadcast to the other
 * processors in the communicator.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
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
 */
void PPEXSILoadMatrix(
    PPEXSIPlan    plan,
    int           nrows,                        
    int           nnz,                          
    int           nnzLocal,                     
    int           numColLocal,                  
    int*          colptrLocal,                  
    int*          rowindLocal,                  
    double*       HnzvalLocal,                  
    int           isSIdentity,                  
    double*       SnzvalLocal );




/**
 * @brief Simplified driver for solving Kohn-Sham DFT.
 *
 * This function contains both the inertia counting step for estimating
 * the chemical potential, and the Newton's iteration for updating the
 * chemical potential.  Heuristics are built into this routine.  Expert
 * users and developers can modify this routine to obtain better
 * heuristics.
 *
 * The input parameter options are controlled through the structure
 * PPEXSIDFTOptions.  The default value can be obtained through
 * PPEXSISetDefaultDFTOptions.
 *
 *
 * **Basic strategy of the heuristics**
 *
 * - If isInertiaCount == 1, then the inertia counting procedure is
 * invoked until the chemical potential interval is refined from size
 * (muMax0-muMin0) to muInertiaTolerance, or the estimated band gap is
 * larger than muInertiaTolerance.
 * - If the change of step in Newton's iteration is larger than
 * muPEXSISafeGuard, the the inertia counting procedure is invoked
 * again, starting from (muMin0, muMax0).  If Newton's iteration fails
 * again, the subroutine returns error message with info = 1.
 * - The number of shifts in the inertia count is always automatically
 * chosen to be (mpisize / npPerPole). This minimizes the wall clock
 * time of the inertia counting procedure.
 *
 *
 * @todo Estimate deltaE automatically from H and S.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 * @param[in] numElectronExact (global) Exact number of electrons, i.e.
 * \f$N_e(\mu_{\mathrm{exact}})\f$.
 * @param[in] options (global) The options containing input parameters
 * for the driver.  
 * @param[out]  DMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero value
 * of density matrix in CSC format.
 * @param[out] EDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of energy density matrix in CSC format.
 * @param[out] FDMnzvalLocal (local)  Dimension: nnzLocal.  Nonzero
 * value of free energy density matrix in CSC format.
 * @param[out] muPEXSI      (global) Chemical potential after the last
 * iteration.
 * - In the case that convergence is reached within maxPEXSIIter steps, the
 * value of muPEXSI is the last mu used to achieve accuracy within
 * numElectronPEXSITolerance.
 * - In the case that convergence is not reached within maxPEXSIIter steps,
 * and the update from Newton's iteration does not exceed
 * muPEXSISafeGuard, the value of muPEXSI is the last mu plus the update
 * from Newton's iteration.
 *
 * @param[out] numElectronPEXSI (global) Number of electrons
 * evaluated at the last step.  
 * @note In the case that convergence is not reached within maxPEXSIIter steps,
 * and numElectron does not correspond to the number of electrons
 * evaluated at muPEXSI.
 * @param[out] muMinInertia (global) Lower bound for mu after the last
 * inertia count procedure.
 * @param[out] muMaxInertia (global) Upper bound for mu after the last
 * inertia count procedure.
 * @param[out] numTotalInertiaIter (global) Number of total inertia
 * counting procedure.
 * @param[out] numTotalPEXSIIter (global) Number of total %PEXSI
 * evaluation procedure.
 */
void PPEXSIDFTDriver(
    /* Input parameters */
    PPEXSIPlan        plan,
    double            numElectronExact,
    PPEXSIDFTOptions  options,
    /* Output parameters */
		double*      DMnzvalLocal,                  
		double*     EDMnzvalLocal,                 
		double*     FDMnzvalLocal,                
		double*       muPEXSI,                   
		double*       numElectronPEXSI,         
    double*       muMinInertia,              
		double*       muMaxInertia,             
		int*          numTotalInertiaIter,   
		int*          numTotalPEXSIIter,   
    int*          info );    


/**
 * @brief Release the memory used by %PEXSI.
 *
 * @param[in] plan (local) The plan holding the internal data structure for the %PEXSI
 * data structure.
 */
void PPEXSIPlanFinalize( PPEXSIPlan plan );


#ifdef __cplusplus
}// extern "C"
#endif

#endif // _C_PEXSI_NEW_INTERFACE_H_
