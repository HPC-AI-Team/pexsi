/// @file ppexsi.hpp
/// @brief Main class for parallel PEXSI.
/// @author Lin Lin
/// @date 2012-11-22
#ifndef _PPEXSI_HPP_
#define _PPEXSI_HPP_
#include  "environment_impl.hpp"
#include  "sparse_matrix_impl.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "pole.hpp"
#include	"mpi_interf.hpp"
#include  "superlu_dist_interf.hpp"
#include	"pselinv.hpp"

namespace PEXSI{

/// @class PPEXSIData
///
/// @brief Main class for parallel PEXSI.
///
/// PPEXSI uses SuperLU_DIST for parallel factorization and PSelInv for
/// the selected inversion.  
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

  const Grid*           gridPole_;              // Outer layer communicator. Also used for distributing the DistSparseMatrix.  Each DistSparseMatrix is replicated in the row (numPoleGroup) direction of gridPole.
	const Grid*           gridSelInv_;            // Inner layer communicator for SelInv
	const SuperLUGrid*    gridSuperLU_;           // Inner layer communicator for SuperLU factorization

	SuperNode             super_;                 // Supernode partition

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
	PPEXSIData( const PEXSI::Grid* g, Int nprow, Int npcol );
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
	///	@param[in] numPole Number of poles for the pole expansion
	///	@param[in] temperature  Temperature (in the unit of K)
	///	@param[in] numElectronExact Exact number of electrons
	/// @param[in] gap Band gap (in the unit of au)
	/// @param[in] deltaE Upperbound of the spectrum width
	/// @param[in] mu0 Initial guess of chemical potential (in the unit of au)
	/// @param[in] muMin Initial guess for the lower bound of the chemical
	/// potential.
	/// @param[in] muMax Initial guess for the upper bound of the chemical
	/// potential.
	/// @param[in] HMat Hamiltonian matrix saved in distributed compressed
	/// sparse column format. See DistSparseMatrix.
	/// @param[in] SMat Overlap matrix saved in distributed compressed
	/// sparse column format. See DistSparseMatrix.
	///
	/// **Note**: If SMat.size == 0, SMat is treated as an identity matrix.
	/// 
	/// @param[in] muMaxIter Maximum iteration number for chemical
	/// potential
	/// @param[in] poleTolerance Skip the pole if weight is too small
	/// @param[in] numElectronTolerance  Stopping criterion for the mu iteration 
	/// @param[in] ColPerm   Permutation method used for SuperLU_DIST
	/// @param[in] isFreeEnergyDensityMatrix Whether to compute the Helmholtz free energy matrix
	/// @param[in] isEnergyDensityMatrix Whether to compute the energy density matrix for force
	/// @param[in] isDerivativeTMatrix Whether to compute the derivative of the
	/// single particle density matrix with respect to the temperature.
	/// @param[out] muList Convergence history of the chemical potential
	/// @param[out] numElectronList Convergence history of the number of
	/// electrons
	/// @param[out] numElectronDrvMuList Convergence history of the
	/// derivative of electrons with respect to the chemical potential.
	/// @param[out] isConverged Whether PEXSI has converged.
	void Solve( 
			Int  numPole, 
			Real temperature,
			Real numElectronExact,
			Real gap,
			Real deltaE,
			Real mu0,
			Real muMin,
			Real muMax,
			const DistSparseMatrix<Real>&  HMat,
		 	const DistSparseMatrix<Real>&  SMat,
		 	Int  muMaxIter,
			Real poleTolerance,
			Real numElectronTolerance,
			std::string         ColPerm,
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
	void PPEXSIData::CalculateNegativeInertia( 
			const std::vector<Real>&       shiftVec, 
			std::vector<Real>&             inertiaVec,
			const DistSparseMatrix<Real>&  HMat,
			const DistSparseMatrix<Real>&  SMat,
			std::string                    ColPerm);

};


} // namespace PEXSI
#endif
