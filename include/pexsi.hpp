#ifndef _PEXSI_HPP_
#define _PEXSI_HPP_
#include  "environment_impl.hpp"
#include  "sparse_matrix.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "getpolef.h"
#include  "selinv_Interf.h"

// *********************************************************************
// Main class for performing the pole expansion and selected inversion
// *********************************************************************
namespace PEXSI{

class PEXSIData{
public:
	// *********************************************************************
	// Input parameters
	// *********************************************************************
	Real gap;                    // Band gap (in the unit of au)
	Real temperature;            // Temperature (in the unit of K)
	Real deltaE;                 // an upperbound of the spectrum width
	Int  numPole;                // Number of poles for the pole expansion
	Real mu0;                    // initial guess of chemical potential (in the unit of au)
	Real numElectronExact;       // Exact number of electrons

	Real numElectronTolerance;   // Stopping criterion for the mu iteration 
	Int  muMaxIter;              // Maximum iteration number for mu
	Int  permOrder;              // Order of the permutation

	bool isHelmholtz;                           // Whether to compute the Helmholtz free energy
	bool isForce;                               // Whether to compute the force


	// H, S, Rho shares the same sparsity pattern.
	SparseMatrix<Real>  HMat;                    // Hamiltonian matrix
	SparseMatrix<Real>  SMat;                    // Overlap matrix

	Real totalEnergy;             // The total energy (linear part)
	Real totalFreeEnergy;         // The total Helmholtz free energy (linear part)

	// *********************************************************************
	// Output parameters
	// *********************************************************************
	SparseMatrix<Real>  rhoMat;                  // Density matrix

	std::vector<Real> muList;              // chemical potential
	std::vector<Real> numElectronList;     // number of electrons


	// *********************************************************************
	// Computational variables
	// *********************************************************************

	CpxNumVec  zshift;      // Complex shift for the pole expansion
	CpxNumVec  zweightRho;  // Complex weight for the pole expansion for density
	CpxNumVec  zweightFreeEnergy;  // Complex shift for the pole expansion for Helmholtz free energy
	CpxNumVec  zweightForce;  // Complex weight for the pole expansion for force

public:
	PEXSIData(){}
	~PEXSIData() {}

  void Setup();
	void Solve();
	void UpdateMu();
};

////	Real traceprod(Int Ndof, IntNumVec& colptr_H, IntNumVec& rowind_H, 
////			DblNumVec& nzval1, DblNumVec& nzval2);


} // namespace PEXSI
#endif
