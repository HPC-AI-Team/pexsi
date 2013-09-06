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
/// @file pexsi.hpp
/// @brief Main class of sequential PEXSI.
/// @date 2012-11-02
#ifndef _PEXSI_HPP_
#define _PEXSI_HPP_
#include  "environment_impl.hpp"
#include  "sparse_matrix_impl.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "pole.hpp"
#include  "selinv_interf.hpp"

/// @namespace PEXSI
/// @brief Main namespace for the PEXSI project. 
///
/// Everything used in PEXSI, including the sequential and the parallel
/// version  is under the PEXSI namespace so that PEXSI can be safely
/// included as a subproject.

namespace PEXSI{

/// @class PEXSIData
///
/// @brief Main class for the sequential PEXSI calculation.
class PEXSIData{
private:
	// *********************************************************************
	// Computational variables
	// *********************************************************************
	std::vector<Complex>  zshift_;      // Complex shift for the pole expansion
	std::vector<Complex>  zweightRho_;  // Complex weight for the pole expansion for density
	std::vector<Complex>  zweightHelmholtz_;  // Complex shift for the pole expansion for Helmholtz free energy
	std::vector<Complex>  zweightForce_;  // Complex weight for the pole expansion for force

	SparseMatrix<Real>    rhoMat_;              // Density matrix
	SparseMatrix<Real>    freeEnergyDensityMat_;     // Helmholtz free energy density matrix
	SparseMatrix<Real>    energyDensityMat_;         // Energy density matrix for computing the Pulay force

	/// @brief Calculate the new chemical potential based on the history
	/// using Newton's method.
	Real CalculateChemicalPotentialNewton( 
			const Int iter, 
			const Real numElectronExact, 
			const Real numElectronTolerance,
			const std::vector<Real>& muList,
			const std::vector<Real>& numElectronList );

public:
	PEXSIData(){}

	~PEXSIData() {}

	/// @brief Solve is the main subroutine for PEXSI.
	///
	/// Compute the single particle density matrix, the Helmholtz free
	/// energy density matrix, and the energy density matrix (for
	/// computing the Pulay force) simultaneously.   These matrices can be
	/// called later via member functions DensityMatrix,
	/// FreeEnergyDensityMatrix, EnergyDensityMatrix.
	///
	///
	///	@param[in] numPole Number of poles for the pole expansion
	///	@param[in] temperature  Temperature (in the unit of K)
	///	@param[in] numElectronExact Exact number of electrons
	/// @param[in] gap Band gap (in the unit of au)
	/// @param[in] deltaE Upperbound of the spectrum width
	/// @param[in] mu0 initial guess of chemical potential (in the unit of au)
	/// @param[in] HMat Hamiltonian matrix saved in compressed sparse
	/// column format.  See SparseMatrix.  
	/// @param[in] SMat Overlap matrix saved in compressed sparse column
	/// format. See SparseMatrix.
	/// @param[in] muMaxIter Maximum iteration number for chemical
	/// potential
	/// @param[in] poleTolerance Skip the pole if weight is too small
	/// @param[in] numElectronTolerance  Stopping criterion for the mu iteration 
	/// @param[in] ColPerm   Permutation method used for SelInv. So far
	/// only MMD_AT_PLUS_A is supported.
	/// @param[in] isFreeEnergyDensityMatrix Whether to compute the Helmholtz free energy matrix
	/// @param[in] isEnergyDensityMatrix Whether to compute the energy density matrix for force
	/// @param[out] muList Convergence history of the chemical potential
	/// @param[out] numElectronList Convergence history of the number of
	/// @param[out] isConverged Whether PEXSI has converged.
	/// electrons
	void Solve( 
			Int  numPole, 
			Real temperature,
			Real numElectronExact,
			Real gap,
			Real deltaE,
			Real mu0,
			const SparseMatrix<Real>&  HMat,
		 	const SparseMatrix<Real>&  SMat,
		 	Int  muMaxIter,
			Real poleTolerance,
			Real numElectronTolerance,
			std::string         ColPerm,
			bool isFreeEnergyDensityMatrix, 
			bool isEnergyDensityMatrix,
			std::vector<Real>&	muList,
			std::vector<Real>&  numElectronList,
			bool&               isConverged
			);

	/// @brief DensityMatrix returns the single particle density matrix.
	SparseMatrix<Real>& DensityMatrix () { return rhoMat_; }

	/// @brief CalculateNumElectron computes the number of electrons given
	/// the current density matrix.
	///
	/// @param[in] SMat overlap matrix.
	///
	/// @return The number of electrons Tr[\rho S]
	Real CalculateNumElectron( const SparseMatrix<Real>& SMat );

	/// @brief CalculateTotalEnergy computes the total energy (band energy
	/// part only).
	///
	/// @param[in] HMat Hamilotian matrix.
	///
	/// @return The total energy Tr[ H \rho ]. 
  Real CalculateTotalEnergy( const SparseMatrix<Real>& HMat );

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
	Real CalculateFreeEnergy( const SparseMatrix<Real>& SMat );

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
			const SparseMatrix<Real>& HDerivativeMat,  
			const SparseMatrix<Real>& SDerivativeMat ); 

};



} // namespace PEXSI
#endif
