/// @file ppexsi.hpp
/// @brief Main class for parallel PEXSI.
/// @author Lin Lin
/// @date 2012-11-22
#ifndef _PPEXSI_HPP_
#define _PPEXSI_HPP_
#include  "environment_impl.hpp"
#include  "sparse_matrix.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "pole.hpp"
#include	"mpi_interf.hpp"
#include  "superlu_dist_interf.hpp"
#include	"pselinv.hpp"

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
	std::vector<Complex>  zweightFreeEnergy_;  // Complex shift for the pole expansion for Helmholtz free energy
	std::vector<Complex>  zweightForce_;  // Complex weight for the pole expansion for force
	Int                   numPoleUsed_;            // Number of poles after truncation

  const Grid*           gridPole_;              // Outer layer communicator. Also used for distributing the DistSparseMatrix.  Each DistSparseMatrix is replicated in the row (numPoleGroup) direction of gridPole.
	const Grid*           gridSelInv_;            // Inner layer communicator for SelInv
	const SuperLUGrid*    gridSuperLU_;           // Inner layer communicator for SuperLU factorization

	SuperNode             super_;                 // Supernode partition

	std::vector<Real>          muList_;              // chemical potential
	std::vector<Real>          numElectronList_;     // number of electrons

	DistSparseMatrix<Real>     rhoMat_;                   // Density matrix
	DistSparseMatrix<Real>     freeEnergyDensityMat_;     // Helmholtz free energy density matrix
	DistSparseMatrix<Real>     energyDensityMat_;         // Energy density matrix for computing the Pulay force

	// *********************************************************************
	// Member functions
	// *********************************************************************
	Real UpdateChemicalPotential( const Int iter );

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
	/// @param[in] mu0 initial guess of chemical potential (in the unit of au)
	/// @param[in] HMat Hamiltonian matrix saved in distributed compressed
	/// sparse column format. See DistSparseMatrix.
	/// @param[in] SMat Overlap matrix saved in distributed compressed
	/// sparse column format. See DistSparseMatrix.
	/// @param[in] muMaxIter Maximum iteration number for chemical
	/// potential
	/// @param[in] numElectronTolerance  Stopping criterion for the mu iteration 
	/// @param[in] isFreeEnergyDensityMatrix Whether to compute the Helmholtz free energy matrix
	/// @param[in] isEnergyDensityMatrix Whether to compute the energy density matrix for force
	/// @param[out] muList Convergence history of the chemical potential
	/// @param[out] numElectronList Convergence history of the number of
	/// electrons
	void Solve( 
			Int  numPole, 
			Real temperature,
			Real numElectronExact,
			Real gap,
			Real deltaE,
			Real mu0,
			const DistSparseMatrix<Real>&  HMat,
		 	const DistSparseMatrix<Real>&  SMat,
		 	Int  muMaxIter,
			Real numElectronTolerance,
			bool isFreeEnergyDensityMatrix, 
			bool isEnergyDensityMatrix,
			std::vector<Real>&	muList,
			std::vector<Real>&  numElectronList
			);
			


	/// @brief DensityMatrix returns the single particle density matrix.
	DistSparseMatrix<Real>& DensityMatrix () { return rhoMat_; }

	/// @brief DensityMatrix returns the free energy density matrix.
	DistSparseMatrix<Real>& FreeEnergyDensityMatrix () { return freeEnergyDensityMat_; }

	/// @brief DensityMatrix returns the energy density matrix for
	/// computing the force.
	DistSparseMatrix<Real>& EnergyDensityMatrix() { return energyDensityMat_; }

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
	Real CalculateFreeEnergy( const DistSparseMatrix<Real>& HMat );

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

	/// @brief CalculateNumElectron computes the number of electrons given
	/// the current density matrix.
	///
	/// @param[in] SMat overlap matrix.
	///
	/// @return The number of electrons Tr[\rho S]
	Real CalculateNumElectron( const DistSparseMatrix<Real>& SDerivativeMat );
};


} // namespace PEXSI
#endif
