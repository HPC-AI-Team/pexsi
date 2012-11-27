#ifndef _PEXSI_HPP_
#define _PEXSI_HPP_
#include  "environment_impl.hpp"
#include  "sparse_matrix.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "pole.hpp"
#include	"mpi_interf.hpp"
#include  "superlu_dist_interf.hpp"
#include	"pselinv.hpp"

namespace PEXSI{

// *********************************************************************
// Main class for parallel PEXSI
// *********************************************************************
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


public:
	PPEXSIData( const PEXSI::Grid* g, Int nprow, Int npcol );
	~PPEXSIData();

	void Solve( 
			Int  numPole,                // Number of poles for the pole expansion
			Real temperature,            // Temperature (in the unit of K)
			Real numElectronExact,       // Exact number of electrons
			Real gap,                    // Band gap (in the unit of au)
			Real deltaE,                 // an upperbound of the spectrum width
			Real mu0,                    // initial guess of chemical potential (in the unit of au)
			const DistSparseMatrix<Real>&  HMat, // Hamiltonian matrix 
			const DistSparseMatrix<Real>&  SMat, // Overlap matrix
			Int  muMaxIter,              // Maximum iteration number for mu
			Real numElectronTolerance,   // Stopping criterion for the mu iteration 
			bool isFreeEnergyDensityMatrix,           // Whether to compute the Helmholtz free energy matrix
			bool isEnergyDensityMatrix                // Whether to compute the energy density matrix for force
			);

	std::vector<Real>& MuHistory( ) { return muList_; }

	std::vector<Real>& NumElectronHistory( ) { return numElectronList_; }

	DistSparseMatrix<Real>& DensityMatrix () { return rhoMat_; }

	DistSparseMatrix<Real>& FreeEnergyDensityMatrix () { return freeEnergyDensityMat_; }

	DistSparseMatrix<Real>& EnergyDensityMatrix() { return energyDensityMat_; }

  Real CalculateEnergy( const DistSparseMatrix<Real>& HMat );

	Real CalculateFreeEnergy( const DistSparseMatrix<Real>& HMat );

	Real CalculateForce( 
			const DistSparseMatrix<Real>& HDerivativeMat,  
			const DistSparseMatrix<Real>& SDerivativeMat ); 
};


} // namespace PEXSI
#endif
