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
/// @file new_ppexsi.hpp
/// @brief New main class for parallel %PEXSI.
/// @date 2014-03-09  Revise for the new interface.
#ifndef _NEW_PPEXSI_HPP_
#define _NEW_PPEXSI_HPP_
#include  "environment.hpp"
#include  "sparse_matrix_impl.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "pole.hpp"
#include	"mpi_interf.hpp"
#include  "superlu_dist_interf.hpp"
#include	"pselinv.hpp"
#include  "c_pexsi_new_interface.h"

namespace PEXSI{

	/// @class PPEXSINewData
	///
	/// @brief Main class for parallel %PEXSI.
	///
	class PPEXSINewData{
	private:
		// *********************************************************************
		// Computational variables
		// *********************************************************************

		std::vector<Complex>  zshift_;      // Complex shift for the pole expansion
		std::vector<Complex>  zweightRho_;  // Complex weight for the pole expansion for density
		std::vector<Complex>  zweightRhoDrvMu_;  // Complex weight for the pole expansion for derivative of the Fermi-Dirac with respect to the chemical potential
		std::vector<Complex>  zweightRhoDrvT_;   // Complex weight for the pole expansion for derivative of the Fermi-Dirac with respect to the temperature T (1/beta, in au)
		std::vector<Complex>  zweightHelmholtz_;  // Complex shift for the pole expansion for Helmholtz free energy
		std::vector<Complex>  zweightForce_;  // Complex weight for the pole expansion for force

    // Outer layer communicator. Also used for distributing the
    // DistSparseMatrix.  Each DistSparseMatrix is replicated in the row
    // (numPoleGroup) direction of gridPole.
		const GridType*           gridPole_;          
		const GridType*           gridSelInv_;        // Inner layer communicator for SelInv
		const SuperLUGrid*        gridSuperLU_;       // Inner layer communicator for SuperLU factorization

		SuperNodeType             super_;             // Supernode partition

    DistSparseMatrix<Real>     HRealMat_;
    DistSparseMatrix<Real>     SRealMat_;

		DistSparseMatrix<Real>     rhoMat_;                   // Density matrix 
		DistSparseMatrix<Real>     rhoDrvMuMat_;              // Derivative of the Fermi-Dirac with respect to mu
		DistSparseMatrix<Real>     rhoDrvTMat_;               // Derivative of the Fermi-Dirac with respect to T
		DistSparseMatrix<Real>     freeEnergyDensityMat_;     // Helmholtz free energy density matrix
		DistSparseMatrix<Real>     energyDensityMat_;         // Energy density matrix for computing the Pulay force

		// Saves all the indices of diagonal elements in H, so that
		// H.nzvalLocal(diagIdxLocal_[j]) are diagonal elements for all j.
		// This is manly used when S is implicitly given as an identity matrix.
		std::vector<Int>           diagIdxLocal_;    

		// *********************************************************************
		// Member functions
		// *********************************************************************
		/// @brief Calculate the new chemical potential using the Newton's
		/// method with finite difference calculation of the derivative.
		Real CalculateChemicalPotentialNewtonFD(
				const Int iter, 
				const Real numElectronExact, 
				const Real numElectronTolerance,
				const std::vector<Real>& muList,
				const std::vector<Real>& numElectronList );

		/// @brief Calculate the new chemical potential using the Newton's
		/// method with bisection as safeguard
		/// @param[in] numElectronExact  Exact number of electrons.
		/// @param[in] numElectron       Computed number of electrons.
		/// @param[in] numElectronDrvMu  The derivative of the number of
		/// electrons with respect to the chemical potential.
		/// @param[in] mu Chemical potential and its update (output).
		/// @param[in] muMin Lower bound of the chemical potential.
		/// @param[in] muMax Upper bound of the chemical potential.
		/// @return Update of the chemical potential.
		Real CalculateChemicalPotentialNewtonBisection(
				const Real numElectronExact, 
				const Real numElectron,
				const Real numElectronDrvMu,
				const Real mu,
				const Real muMin,
				const Real muMax );

	public:
    PPEXSINewData(
        MPI_Comm   comm,
        Int        numProcRow, 
        Int        numProcCol, 
        Int        outputFileIndex );

		~PPEXSINewData();

    void LoadRealSymmetricMatrix(
        Int           nrows,                        
        Int           nnz,                          
        Int           nnzLocal,                     
        Int           numColLocal,                  
        Int*          colptrLocal,                  
        Int*          rowindLocal,                  
        Real*         HnzvalLocal,                  
        Int           isSIdentity,                  
        Real*         SnzvalLocal );


    void DFTDriver(
        Real       numElectronExact,
        Real       temperature,
        Real       gap,
        Real       deltaE,
        Int        numPole, 
        Int        isInertiaCount,
        Int        maxPEXSIIter,
        Real       muMin0,
        Real       muMax0,
        Real       muInertiaTolerance,
        Real       muInertiaExpansion,
        Real       muPEXSISafeGuard,
        Real       numElectronPEXSITolerance,
        Int        matrixType,
        Int        ordering,
        Int        numProcSymbFact,
        Int        verbosity,
        Real&      muPEXSI,                   
        Real&      numElectronPEXSI,         
        Real&      muMinInertia,              
        Real&      muMaxInertia,             
        Int&       numTotalInertiaIter,   
        Int&       numTotalPEXSIIter );


		/// @brief Compute the negative inertia (the number of eigenvalues
		/// below a shift) using real arithemetic factorization routine.
		///
		/// This subroutine computes the negative inertia of the matrix
		///
		/// I = H - shift * S
		///
		/// where I is the same as the number of eigenvalues lambda for
		///
		/// H x = lambda S x
		///
		/// with lambda < shift according to the Sylvester's law of inertia.
		///
		/// @param[in]  shiftVec Shift vectors.
		/// @param[out] inertiaVec Negative inertia count, the same size as
		/// shiftVec.
		/// @param[in] HMat Hamiltonian matrix saved in distributed compressed
		/// sparse column format. See DistSparseMatrix.
		/// @param[in] SMat Overlap matrix saved in distributed compressed
		/// sparse column format. See DistSparseMatrix.
		///
		/// **Note**: If SMat.size == 0, SMat is treated as an identity matrix.
		/// 
		/// @param[in] ColPerm   Permutation method used for SuperLU_DIST
		///
		/// @param[in] numProcSymbFact Number of processors used for parallel
		/// symbolic factorization and PARMETIS/PT-SCOTCH.
		void CalculateNegativeInertiaReal(
				const std::vector<Real>&       shiftVec, 
				std::vector<Real>&             inertiaVec,
				const DistSparseMatrix<Real>&  HMat,
				const DistSparseMatrix<Real>&  SMat,
				std::string                    ColPerm,
				Int                            numProcSymbFact );


    // *********************************************************************
    // Access data
    // *********************************************************************

    const GridType*  GridPole() const {return gridPole_;}


	}; // PPEXSINewData


} // namespace PEXSI
#endif
