/// @file ppexsi.cpp
/// @brief Implementation of the parallel version of PEXSI.
/// @author Lin Lin
/// @date 2012-11-20
#include "ppexsi.hpp"

namespace PEXSI{

PPEXSIData::PPEXSIData	( const PEXSI::Grid* g, Int nprow, Int npcol ): gridPole_(g)
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::PPEXSIData");
#endif
	if( nprow != npcol ){
		throw std::runtime_error( "PSelInv only allows to use square grid." );
	}
	if( gridPole_->numProcCol != nprow * npcol ){
		throw std::runtime_error( "The number of processors numProcCol do not match nprow * npcol." );
	}
	gridSuperLU_  = new SuperLUGrid( gridPole_->rowComm, nprow, npcol );
	gridSelInv_   = new Grid( gridPole_->rowComm, nprow, npcol );

  
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::PPEXSIData  ----- 


PPEXSIData::~PPEXSIData	(  )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::~PPEXSIData");
#endif
	if( gridSuperLU_ != NULL ){
		delete gridSuperLU_;
	}
	
	if( gridSelInv_ != NULL ){
		delete gridSelInv_;
	}
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::~PPEXSIData  ----- 


// Main subroutine for the electronic structure calculation
void PPEXSIData::Solve( 
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
		){
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::Solve");
#endif


	// *********************************************************************
	// Check the input parameters
	// *********************************************************************
	if( numPole % 2 != 0 ){
		throw std::logic_error( "Must be even number of poles!" );
	}

	// TODO Check H and S have the same pattern

	// TODO Check H and S agree with gridPole_


	// *********************************************************************
	// Initialize
	// *********************************************************************
	muList.clear();
	numElectronList.clear();

	DistSparseMatrix<Complex>  AMat;              // A = H - z * S
	DistSparseMatrix<Real>& rhoMat = rhoMat_;     // rename for convenience
	// Copy the pattern
	CopyPattern( HMat, AMat );
	CopyPattern( HMat, rhoMat );

	SetValue( AMat.nzvalLocal, Z_ZERO );          // Symbolic factorization does not need value

	statusOFS << "AMat.nnzLocal = " << AMat.nnzLocal << std::endl;

	SuperLUMatrix              luMat( *gridSuperLU_ );  // SuperLU matrix.

	// *********************************************************************
	// Symbolic factorization.  
	// Each numPoleGroup perform independently
	// *********************************************************************
	luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
	luMat.SymbolicFactorize();
	luMat.SymbolicToSuperNode( super_ );
	luMat.DestroyAOnly();


	Real muNow = mu0;
	Real numElectronNow;

	Real timeMuSta, timeMuEnd;



	// Iteration with chemical potentials
	for(Int iter = 0; iter < muMaxIter; iter++){
		GetTime( timeMuSta );

		{
			std::ostringstream msg;
			msg << "Iteration " << iter << ", mu = " << muNow;
			PrintBlock( statusOFS, msg.str() );
		}

		// Reinitialize the variables
		SetValue( rhoMat.nzvalLocal, 0.0 );

		// Initialize the number of electrons
		numElectronNow  = 0.0;

		//Initialize the pole expansion
		zshift_.clear(); zshift_.resize( numPole );
		zweightRho_.clear(); zweightRho_.resize( numPole );

		GetPoleDensity( &zshift_[0], &zweightRho_[0],
				numPole, temperature, gap, deltaE, muNow ); 

#if ( _DEBUGlevel_ >= 0 )
		statusOFS << "zshift" << std::endl << zshift_ << std::endl;
		statusOFS << "zweightRho" << std::endl << zweightRho_ << std::endl;
#endif

		// for each pole, perform LDLT factoriation and selected inversion
		Real timeSta, timeEnd;
		Real timePoleSta, timePoleEnd;

		for(Int l = 0; l < numPole; l++){
			if( MYROW( gridPole_ ) == PROW( l, gridPole_ ) ){

				GetTime( timePoleSta );
				statusOFS << "Pole " << l << " processing..." << std::endl;
#if ( _DEBUGlevel_ >= 0 )
				statusOFS << "zshift = " << zshift_[l] << ", " 
					<< "zweightRho = " << zweightRho_[l] << std::endl;
#endif
				// FIXME magic number here
				if( std::abs( zweightRho_[l] ) < 1e-8 ){
					statusOFS << "|zweightRho| < 1e-8, pass this pole" << std::endl;
					continue;
				}


				for( Int i = 0; i < HMat.nnzLocal; i++ ){
					AMat.nzvalLocal(i) = HMat.nzvalLocal(i) - zshift_[l] * SMat.nzvalLocal(i);
				}


				// *********************************************************************
				// Factorization
				// *********************************************************************
				// Important: the distribution in pzsymbfact is going to mess up the
				// A matrix.  Recompute the matrix A here.
				luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );

				Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

				GetTime( timeTotalFactorizationSta );

				// Data redistribution
				luMat.Distribute();

				// Numerical factorization
				luMat.NumericalFactorize();
				luMat.DestroyAOnly();

				GetTime( timeTotalFactorizationEnd );

				statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 

				// *********************************************************************
				// Selected inversion
				// *********************************************************************
				Real timeTotalSelInvSta, timeTotalSelInvEnd;
				GetTime( timeTotalSelInvSta );

				PMatrix PMloc( gridSelInv_, &super_ ); // A^{-1} in PMatrix format

				luMat.LUstructToPMatrix( PMloc );

				PMloc.ConstructCommunicationPattern();

				PMloc.PreSelInv();

				PMloc.SelInv();

				GetTime( timeTotalSelInvEnd );

				statusOFS << "Time for total selected inversion is " <<
					timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;

				// *********************************************************************
				// Postprocessing
				// *********************************************************************

				DistSparseMatrix<Complex>  AinvMat;       // A^{-1} in DistSparseMatrix format

				Real timePostProcessingSta, timePostProcessingEnd;

				GetTime( timePostProcessingSta );

				PMloc.PMatrixToDistSparseMatrix( AMat, AinvMat );
				
#if ( _DEBUGlevel_ >= 0 )
				statusOFS << "rhoMat.nnzLocal = " << rhoMat.nnzLocal << std::endl;
				statusOFS << "AinvMat.nnzLocal = " << AinvMat.nnzLocal << std::endl;
#endif

				// Update the density matrix
				for( Int i = 0; i < rhoMat.nnzLocal; i++ ){
					rhoMat.nzvalLocal(i) += 
						zweightRho_[l].real() * AinvMat.nzvalLocal(i).imag() + 
						zweightRho_[l].imag() * AinvMat.nzvalLocal(i).real();
				}

				GetTime( timePostProcessingEnd );

				statusOFS << "Time for postprocessing is " <<
					timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;

				GetTime( timePoleEnd );

				statusOFS << "Time for pole " << l << " is " <<
					timePoleEnd - timePoleSta << " [s]" << std::endl << std::endl;

			} // if I am in charge of this pole
		} // for(l)

		// Reduce the density matrix across the processor rows in gridPole_

		DblNumVec nzvalRhoMatLocal = rhoMat.nzvalLocal;
		SetValue( rhoMat.nzvalLocal, 0.0 );
		
		mpi::Allreduce( nzvalRhoMatLocal.Data(), rhoMat.nzvalLocal.Data(),
				rhoMat.nnzLocal, MPI_SUM, gridPole_->colComm );

		// All processors groups compute the number of electrons

		numElectronNow = CalculateNumElectron( SMat );

		muList.push_back(muNow);
		numElectronList.push_back( numElectronNow );

		Print( statusOFS, "Computed number of electron = ", numElectronNow );
		Print( statusOFS, "Exact number of electron    = ", numElectronExact );

		// TODO
		// Reduce band energy
		// Reduce Helmholtz free energy
		// Reduce force

		if( std::abs( numElectronExact - numElectronList[iter] ) <
				numElectronTolerance ){
			break;
		}

		muNow = CalculateChemicalPotential( 
				iter, numElectronExact, numElectronTolerance,
				muList, numElectronList );

		GetTime( timeMuEnd );

		statusOFS << std::endl << "Time for mu iteration " << iter << " is " <<
			timeMuEnd - timeMuSta << " [s]" << std::endl << std::endl;
	} // for ( iteration of the chemical potential )

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::Solve----- 


  
Real PPEXSIData::CalculateChemicalPotential	( 
			const Int iter, 
			const Real numElectronExact, 
			const Real numElectronTolerance, 
			const std::vector<Real>& muList,
			const std::vector<Real>& numElectronList )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateChemicalPotential");
#endif
  // FIXME Magic number here
	Real  muMin = -5.0, muMax = 5.0, muMinStep = 0.01;;
	Real  muNew;

	if( iter == 0 ){
		if( numElectronExact > numElectronList[iter] ){
			muNew = muList[iter] + muMinStep;
		}
		else{
			muNew = muList[iter] - muMinStep;
		}
	}
	else{
		if( std::abs(numElectronList[iter] -  numElectronList[iter-1])
		    < numElectronTolerance ){
			statusOFS << "The number of electrons did not change." << std::endl;
			if( numElectronExact > numElectronList[iter] ){
				muNew = muList[iter] + muMinStep;
			}
			else{
				muNew = muList[iter] - muMinStep;
			}
		}
		else {
			muNew = muList[iter] + (muList[iter] - muList[iter-1]) / 
				( numElectronList[iter] - numElectronList[iter-1] ) *
				( numElectronExact - numElectronList[iter] );
			if( muNew < muMin || muNew > muMax ){
				statusOFS << "muNew = " << muNew << " is out of bound ["
					<< muMin << ", " << muMax << "]" << std::endl;
				if( numElectronExact > numElectronList[iter] ){
					muNew = muList[iter] + muMinStep;
				}
				else{
					muNew = muList[iter] - muMinStep;
				}
			}
		} // if ( numElectron changed )
	} // if (iter == 0)

#ifndef _RELEASE_
	PopCallStack();
#endif

	return muNew;
} 		// -----  end of method PPEXSIData::CalculateChemicalPotential  ----- 


Real
PPEXSIData::CalculateNumElectron	( const DistSparseMatrix<Real>& SMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateNumElectron");
#endif
	Real numElecLocal = 0.0, numElec = 0.0;
	
	// TODO Check SMat and rhoMat has the same sparsity

	numElecLocal = blas::Dot( SMat.nnzLocal, SMat.nzvalLocal.Data(),
			1, rhoMat_.nzvalLocal.Data(), 1 );
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "numElecLocal = " << numElecLocal << std::endl;
#endif


	mpi::Allreduce( &numElecLocal, &numElec, 1, MPI_SUM, rhoMat_.comm ); 


#ifndef _RELEASE_
	PopCallStack();
#endif

	return numElec;
} 		// -----  end of method PPEXSIData::CalculateNumElectron  ----- 
} //  namespace PEXSI
