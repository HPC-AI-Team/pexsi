#ifndef _PEXSI_HPP_
#define _PEXSI_HPP_
#include  "environment_impl.hpp"
#include  "sparse_matrix.hpp"
#include  "numvec_impl.hpp"
#include  "utility.hpp"
#include  "pole.hpp"
#include	"mpi_interf.hpp"
#include	"pselinv.hpp"

namespace PEXSI{

// *********************************************************************
// Main class for parallel PEXSI
// *********************************************************************
class PPEXSIData{
private:
  
public:
	// *********************************************************************
	// Input parameters
	// *********************************************************************
	Real gap;                    // Band gap (in the unit of au)
	Real temperature;            // Temperature (in the unit of K)
	Real deltaE;                 // an upperbound of the spectrum width
	Int  numPole;                // Number of poles for the pole expansion
	Real poleTolerance;          // Truncation tolerance for the absolute value of the weight
	Real mu0;                    // initial guess of chemical potential (in the unit of au)
	Real numElectronExact;       // Exact number of electrons

	Real numElectronTolerance;   // Stopping criterion for the mu iteration 
	Int  muMaxIter;              // Maximum iteration number for mu

	bool isHelmholtz;                           // Whether to compute the Helmholtz free energy
	bool isForce;                               // Whether to compute the force

	// H, S, Rho shares the same sparsity pattern.
	DistSparseMatrix<Real>  HMat;                    // Hamiltonian matrix
	DistSparseMatrix<Real>  SMat;                    // Overlap matrix

	Real totalEnergy;             // The total energy (linear part)
	Real totalFreeEnergy;         // The total Helmholtz free energy (linear part)

	// *********************************************************************
	// Output parameters
	// *********************************************************************
	DistSparseMatrix<Real>  rhoMat;        // Density matrix

	std::vector<Real> muList;              // chemical potential
	std::vector<Real> numElectronList;     // number of electrons


	// *********************************************************************
	// Computational variables
	// *********************************************************************

	std::vector<Complex>  zshift;      // Complex shift for the pole expansion
	std::vector<Complex>  zweightRho;  // Complex weight for the pole expansion for density
	std::vector<Complex>  zweightFreeEnergy;  // Complex shift for the pole expansion for Helmholtz free energy
	std::vector<Complex>  zweightForce;  // Complex weight for the pole expansion for force
	Int        numPoleUsed;            // Number of poles after truncation

public:
	PPEXSIData(){}
	~PPEXSIData() {}

  void Setup();
	
	void Solve();

	Real ProductTrace	( const DblNumVec& nzval1, const DblNumVec& nzval2 );

	Real UpdateChemicalPotential( const Int iter );
};


} // namespace PEXSI
#endif
