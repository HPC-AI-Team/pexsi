#ifndef _PEXSI_HPP_
#define _PEXSI_HPP_
#include  "environment_impl.hpp"
#include  "sparse_matrix.hpp"
#include  "numvec_impl.hpp"
//#include  "getpolef.h"
//#include  "selinv_Interf.h"

// *********************************************************************
// Main class for performing the pole expansion and selected inversion
// *********************************************************************
namespace PEXSI{

class PEXSI{
public:

	// *********************************************************************
	// Input parameters
	// *********************************************************************
	Real gap_;                    // Band gap (in the unit of au)
	Real temperature_;            // Temperature (in the unit of K)
	Real deltaE_;                 // an upperbound of the spectrum width
	Int  numPole_;                // Number of poles for the pole expansion
	Real mu0_;                    // initial guess of chemical potential (in the unit of au)
	Real numElectronExact_;       // Exact number of electrons

	Int  muMaxIter_;              // Maximum iteration number for mu
	Int  permOrder_;            // Order of the permutation

	bool isHelmholtz_;                           // Whether to compute the Helmholtz free energy
	bool isForce_;                               // Whether to compute the force


	// H, S, Rho shares the same sparsity pattern.
	SparseMatrix<Real>  HMat_;                    // Hamiltonian matrix
	SparseMatrix<Real>  SMat_;                    // Overlap matrix
  
	Real totalEnergy_;             // The total energy (linear part)
	Real totalFreeEnergy_;         // The total Helmholtz free energy (linear part)
	
	// *********************************************************************
	// Output parameters
	// *********************************************************************
	SparseMatrix<Real>  rhoMat_;                  // Density matrix

	std::vector<Real> muList_;              // chemical potential
	std::vector<Real> numElectronList_;     // number of electrons


	// *********************************************************************
	// Computational variables
	// *********************************************************************
	SparseMatrix<Complex>  greenMat_;                // Green's function matrix

	CpxNumVec  zshiftRho_;   // Complex shift for the pole expansion for density
	CpxNumVec  zweightRho_;  // Complex weight for the pole expansion for density

//	CpxNumVec _zshift_hmz;   // Complex shift for the pole expansion for Helmholtz free energy
//	CpxNumVec _zweight_hmz;  // Complex weight for the pole expansion for Helmholtz free energy
//
//	CpxNumVec _zshift_frc;   // Complex shift for the pole expansion for force
//	CpxNumVec _zweight_frc;  // Complex weight for the pole expansion for force
public:
	PEXSI();
	~PEXSI();
	void Setup();
	void Solve();
//	Real traceprod(Int Ndof, IntNumVec& colptr_H, IntNumVec& rowind_H, 
//			DblNumVec& nzval1, DblNumVec& nzval2);
};

} // namespace PEXSI
#endif
