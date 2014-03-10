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
/// @file ppexsi.hpp
/// @brief Main class for parallel %PEXSI.
/// @date 2012-11-22  Initially started.
/// @date 2014-03-09  Revise for the new interface.
#ifndef _PPEXSI_HPP_
#define _PPEXSI_HPP_
#include  "environment.hpp"
#include  "sparse_matrix_impl.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "pole.hpp"
#include	"mpi_interf.hpp"
#include  "superlu_dist_interf.hpp"
#include	"pselinv.hpp"
#include  "c_pexsi_interface.h"

namespace PEXSI{

	/// @class PPEXSIData
	///
	/// @brief Main class for parallel %PEXSI.
	///
	class PPEXSIData{
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

		const GridType*           gridPole_;          // Outer layer communicator. Also used for distributing the DistSparseMatrix.  Each DistSparseMatrix is replicated in the row (numPoleGroup) direction of gridPole.
		const GridType*           gridSelInv_;        // Inner layer communicator for SelInv
		const SuperLUGrid*    gridSuperLU_;           // Inner layer communicator for SuperLU factorization

		SuperNodeType         super_;                 // Supernode partition

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
		PPEXSIData( const PEXSI::GridType* g, Int nprow, Int npcol );
		~PPEXSIData();

		/// @brief Solve is the main subroutine for PPEXSI.
		///
		/// Compute the single particle density matrix, the Helmholtz free
		/// energy density matrix, and the energy density matrix (for
		/// computing the Pulay force) simultaneously.   These matrices can be
		/// called later via member functions DensityMatrix,
		/// FreeEnergyDensityMatrix, EnergyDensityMatrix.
		///
		///
		/// @param[in] numPole Number of poles for the pole expansion
		///	@param[in] temperature  Temperature (in the unit of K)
		///	@param[in] numElectronExact Exact number of electrons
		/// @param[in] gap Band gap (in the unit of au)
		/// @param[in] deltaE Upperbound of the spectrum width
		/// @param[in,out] mu Input: Initial guess of chemical potential (in the
		/// unit of au). Output: The final update of the chemical potential,
		/// for which the number of electrons has NOT been computed.
		/// **Note**: The output mu is NOT the same as last element in muList.
		/// muList[end] corresponds to numElectronList[end].
		/// @param[in,out] muMin Input: Initial guess for the lower bound of the
		/// chemical potential. Output: Lower bound of the chemical potential.
		/// @param[in,out] muMax Input: Initial guess for the upper bound of the
		/// chemical potential. Output: Upper bound of the chemical potential.
		/// @param[in] HMat Hamiltonian matrix saved in distributed compressed
		/// sparse column format. See DistSparseMatrix.
		/// @param[in] SMat Overlap matrix saved in distributed compressed
		/// sparse column format. See DistSparseMatrix.  **Note**: If
		/// SMat.size == 0, SMat is treated as an identity matrix.
		/// @param[in] muMaxIter Maximum iteration number for chemical
		/// potential
		/// @param[in] numElectronTolerance  Stopping criterion for the mu iteration 
		/// @param[in] ColPerm   Permutation method used for SuperLU_DIST
		/// @param[in] numProcSymbFact Number of processors used for parallel
		/// symbolic factorization and PARMETIS/PT-SCOTCH.
		/// @param[in] isFreeEnergyDensityMatrix Whether to compute the Helmholtz free energy matrix
		/// @param[in] isEnergyDensityMatrix Whether to compute the energy density matrix for force
		/// @param[in] isDerivativeTMatrix Whether to compute the derivative of the
		/// single particle density matrix with respect to the temperature.
		/// @param[out] muList Convergence history of the chemical potential
		/// @param[out] numElectronList Convergence history of the number of
		/// electrons
		/// @param[out] numElectronDrvMuList Convergence history of the
		/// derivative of electrons with respect to the chemical potential.
		/// @param[out] isConverged Whether %PEXSI has converged.
		void Solve( 
				Int  numPole, 
				Real temperature,
				Real numElectronExact,
				Real gap,
				Real deltaE,
				Real& mu,
				Real& muMin,
				Real& muMax,
				const DistSparseMatrix<Real>&  HMat,
				const DistSparseMatrix<Real>&  SMat,
				Int  muMaxIter,
				Real numElectronTolerance,
				std::string         ColPerm,
				Int                 numProcSymbFact,
				bool isFreeEnergyDensityMatrix, 
				bool isEnergyDensityMatrix,
				bool isDerivativeTMatrix,
				std::vector<Real>&	muList,
				std::vector<Real>&  numElectronList,
				std::vector<Real>&  numElectronDrvMuList,
				bool&               isConverged
					);



		/// @brief DensityMatrix returns the single particle density matrix.
		DistSparseMatrix<Real>& DensityMatrix () { return rhoMat_; }

		/// @brief DensityMatrix returns the free energy density matrix.
		DistSparseMatrix<Real>& FreeEnergyDensityMatrix () { return freeEnergyDensityMat_; }

		/// @brief DensityMatrix returns the energy density matrix for
		/// computing the force.
		DistSparseMatrix<Real>& EnergyDensityMatrix() { return energyDensityMat_; }

		/// @brief CalculateNumElectron computes the number of electrons given
		/// the current density matrix.
		///
		/// @param[in] SMat overlap matrix.
		///
		/// @return The number of electrons Tr[\rho S]
		Real CalculateNumElectron( const DistSparseMatrix<Real>& SMat );

		/// @brief CalculateNumElectronDrvMu computes the derivative of the
		/// number of electrons with respect to the chemical potential.
		///
		/// @param[in] SMat overlap matrix.
		///
		/// @return The derivative of the number of electrons Tr[f'(H-\muS) S]
		Real CalculateNumElectronDrvMu( const DistSparseMatrix<Real>& SMat );

		/// @brief CalculateNumElectronDrvT computes the derivative of the
		/// number of electrons with respect to the temperature (1/beta, in
		/// atomic unit).
		///
		/// @param[in] SMat overlap matrix.
		///
		/// @return The derivative of the number of electrons Tr[f'(H-\muS) S]
		Real CalculateNumElectronDrvT( const DistSparseMatrix<Real>& SMat );

		/// @brief CalculateTotalEnergy computes the total energy (band energy
		/// part only).
		///
		/// @param[in] HMat Hamilotian matrix.
		///
		/// @return The total energy Tr[ H \rho ]. 
		Real CalculateTotalEnergy( const DistSparseMatrix<Real>& HMat );

		/// @brief CalculateFreeEnergy computes the total Helmholtz free
		/// energy (band energy part only).  
		///
		/// For more information see 
		/// Alavi, A., Kohanoff, J., Parrinello, M., & Frenkel, D. (1994). Ab
		/// initio molecular dynamics with excited electrons. Physical review
		/// letters, 73(19), 2599–2602. 
		///
		/// @param[in] HMat Hamilotian matrix.
		///
		/// @return The Helmholtz free energy Tr[rho_f H]
		Real CalculateFreeEnergy( const DistSparseMatrix<Real>& SMat );

		/// @brief CalculateFreeEnergy computes the force, including the
		/// Hellman-Feynman force and the Pulay force. 
		///
		/// @param[in] HDerivativeMat Derivative of the Hamilotian matrix with
		/// respect to the atomic position R_{i,j}, i = 1, ..., natoms,
		/// j=1,2,3. 
		///
		/// @param[in] HDerivativeMat Derivative of the overlap matrix with
		/// respect to the atomic position R_{i,j}, i = 1, ..., natoms,
		/// j=1,2,3. 
		///
		/// @return The force f_{i,j} with
		/// respect to the atomic position R_{i,j}, i = 1, ..., natoms,
		/// j=1,2,3. 
		Real CalculateForce( 
				const DistSparseMatrix<Real>& HDerivativeMat,  
				const DistSparseMatrix<Real>& SDerivativeMat ); 


		/// @brief Estimate the chemical potential at zero temperature using
		/// quadratic approximation.
		///
		/// @param[in] temperature Temperature (in the unit of K), usually
		/// high (~3000K)
		/// @param[in] mu          Computed chemical potential at high
		/// temperature
		/// @param[in] SMat        Overlap matrix
		///
		/// @return Estimated chemical potential at zero temperature.
		Real EstimateZeroTemperatureChemicalPotential( 
				Real temperature,
				Real mu,
				const DistSparseMatrix<Real>& SMat );


		/// @brief Compute the negative inertia (the number of eigenvalues
		/// below a shift).
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
		void CalculateNegativeInertia( 
				const std::vector<Real>&       shiftVec, 
				std::vector<Real>&             inertiaVec,
				const DistSparseMatrix<Real>&  HMat,
				const DistSparseMatrix<Real>&  SMat,
				std::string                    ColPerm,
				Int                            numProcSymbFact );


		/// @brief Compute the negative inertia (the number of eigenvalues
		/// below a shift) with real arithmetic which directly using SuperLU.
		void CalculateNegativeInertiaReal( 
				const std::vector<Real>&       shiftVec, 
				std::vector<Real>&             inertiaVec,
				const DistSparseMatrix<Real>&  HMat,
				const DistSparseMatrix<Real>&  SMat,
				std::string                    ColPerm,
				Int                            numProcSymbFact );

    /// @brief Estimate the spectral radius of the matrix stencil (H,S).
    ///
    /// The current method uses the locally optimal conjugate gradient
    /// (LOCG) method.
    ///
    /// @note Rigorously speaking we are not computing the spectral
    /// radius, but **the largest generalized eigenvalue**.  For
    /// standard Kohn-Sham problem with pseudopotential this should be
    /// OK.  But this **should be changed** later, at least as a safe
    /// guard, but also may have practical relavance in the case of e.g.
    /// all electron calculation.
    ///
    /// @param[in] method Method for estimating the spectral radius
    /// - = 0     : Locally optimal conjugate gradient (LOCG) method
    ///             (default).  This works for both S being an identity
    ///             matrix or not.  This is the only implemented one
    ///             currently.
		/// @param[in] HMat Hamiltonian matrix saved in distributed compressed
		/// sparse column format. See DistSparseMatrix.
		/// @param[in] SMat Overlap matrix saved in distributed compressed
		/// sparse column format. See DistSparseMatrix. S can be an identity
    /// matrix.
    /// @param[in] v0   Initial starting vector.  
    ///
    /// If v0.m() == 0, then a random vector is used as the initial
    /// starting vector.
    /// @param[in] tol        Relative tolerance for estimating sigma.
    /// @param[in] maxNumIter Maximum number of iterations.
    /// @param[out] numIter The number of iterations.
    /// @param[out] sigma The estimated spectral radius.
		void EstimateSpectralRadius( 
        Int                            method,
				const DistSparseMatrix<Real>&  HMat,
				const DistSparseMatrix<Real>&  SMat,
        const NumVec<Real>&            v0,
        Real                           tol,
        Int                            maxNumIter,
        Int&                           numIter,
        Real&                          sigma );
	}; // PPEXSIData


} // namespace PEXSI
#endif
