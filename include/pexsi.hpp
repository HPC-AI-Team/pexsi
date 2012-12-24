/// @file pexsi.hpp
/// @brief Main class of sequential PEXSI.
/// @author Lin Lin
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
/// 
/// TODO Helmholtz free energy and force support
class PEXSIData{
private:
	// *********************************************************************
	// Computational variables
	// *********************************************************************
	std::vector<Complex>  zshift_;      // Complex shift for the pole expansion
	std::vector<Complex>  zweightRho_;  // Complex weight for the pole expansion for density

	SparseMatrix<Real>    rhoMat_;              // Density matrix

	/// @brief Calculate the new chemical potential based on the history.
	Real CalculateChemicalPotential( 
			const Int iter, 
			const Real numElectronExact, 
			const Real numElectronTolerance,
			const std::vector<Real>& muList,
			const std::vector<Real>& numElectronList );

public:
	// *********************************************************************
	// Input parameters
	// *********************************************************************
//	Real gap;                    // Band gap (in the unit of au)
//	Real temperature;            // Temperature (in the unit of K)
//	Real deltaE;                 // an upperbound of the spectrum width
//	Int  numPole;                // Number of poles for the pole expansion
//	Real poleTolerance;          // Truncation tolerance for the absolute value of the weight
//	Real mu0;                    // initial guess of chemical potential (in the unit of au)
//	Real numElectronExact;       // Exact number of electrons
//
//	Real numElectronTolerance;   // Stopping criterion for the mu iteration 
//	Int  muMaxIter;              // Maximum iteration number for mu
//	Int  permOrder;              // Order of the permutation

	/* -order     :   Reordering strategy:
		 order = -1 (default) : Multiple Minimum Degree Reordering.
		 If METIS is supported, then the following choices
		 are also available:
		 order = 2 : Node Nested Dissection
		 order = 3 : Edge Nested Dissection  */

//	bool isHelmholtz;                           // Whether to compute the Helmholtz free energy
//	bool isForce;                               // Whether to compute the force


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
			std::vector<Real>&  numElectronList
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

};



} // namespace PEXSI
#endif
