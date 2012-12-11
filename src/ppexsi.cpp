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
		Real poleTolerance,
		Real numElectronTolerance,
		std::string         ColPerm,
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
	// rename for convenience
	DistSparseMatrix<Real>& rhoMat       = rhoMat_;     
	DistSparseMatrix<Real>& hmzMat       = freeEnergyDensityMat_;
	DistSparseMatrix<Real>& frcMat     = energyDensityMat_;

	// Copy the pattern
	CopyPattern( HMat, AMat );
	CopyPattern( HMat, rhoMat );
	if( isFreeEnergyDensityMatrix )
		CopyPattern( HMat, hmzMat );
	if( isEnergyDensityMatrix )
		CopyPattern( HMat, frcMat );

	SetValue( AMat.nzvalLocal, Z_ZERO );          // Symbolic factorization does not need value

	statusOFS << "AMat.nnzLocal = " << AMat.nnzLocal << std::endl;
	statusOFS << "AMat.nnz      = " << AMat.nnz      << std::endl;

	SuperLUOptions   luOpt;

	luOpt.ColPerm = ColPerm;

	SuperLUMatrix    luMat( *gridSuperLU_, luOpt );  // SuperLU matrix.

	// *********************************************************************
	// Symbolic factorization.  
	// Each numPoleGroup perform independently
	// *********************************************************************
	luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
	luMat.SymbolicFactorize();
	luMat.SymbolicToSuperNode( super_ );
	luMat.DestroyAOnly();

#if ( _DEBUGlevel_ >= 0 )
	statusOFS << "perm: "    << std::endl << super_.perm     << std::endl;
	statusOFS << "permInv: " << std::endl << super_.permInv  << std::endl;
	statusOFS << "superIdx:" << std::endl << super_.superIdx << std::endl;
	statusOFS << "superPtr:" << std::endl << super_.superPtr << std::endl; 
#endif

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
		if( isFreeEnergyDensityMatrix )
			SetValue( hmzMat.nzvalLocal, 0.0 );
		if( isEnergyDensityMatrix )
			SetValue( frcMat.nzvalLocal, 0.0 );

		// Initialize the number of electrons
		numElectronNow  = 0.0;

		//Initialize the pole expansion
		zshift_.resize( numPole );
		zweightRho_.resize( numPole );

		GetPoleDensity( &zshift_[0], &zweightRho_[0],
				numPole, temperature, gap, deltaE, muNow ); 

		if( isFreeEnergyDensityMatrix ){
			std::vector<Complex>  zshiftTmp( numPole );
			zweightHelmholtz_.resize( numPole );
			GetPoleHelmholtz( &zshiftTmp[0], &zweightHelmholtz_[0], 
				numPole, temperature, gap, deltaE, muNow ); 
			Real norm = 0.0;
			for( Int i = 0; i < numPole; i++ ){
				norm += pow( std::abs( zshiftTmp[i] - zshift_[i] ), 2.0 );
			}
			if( norm > 1e-12 )
				throw std::runtime_error("The pole shifts for rho and for Helmholtz do not match.");
		}

		if( isEnergyDensityMatrix ){
			std::vector<Complex>  zshiftTmp( numPole );
			zweightForce_.resize( numPole );
			GetPoleForce( &zshiftTmp[0], &zweightForce_[0],
				numPole, temperature, gap, deltaE, muNow ); 
			Real norm = 0.0;
			for( Int i = 0; i < numPole; i++ ){
				norm += pow( std::abs( zshiftTmp[i] - zshift_[i] ), 2.0 );
			}
			if( norm > 1e-12 )
				throw std::runtime_error("The pole shifts for rho and for energy density matrix do not match.");
		}


#if ( _DEBUGlevel_ >= 0 )
		statusOFS << "zshift" << std::endl << zshift_ << std::endl;
		statusOFS << "zweightRho" << std::endl << zweightRho_ << std::endl;
#endif

		// Sort the poles according to their weights
		std::vector<Int> sortIdx( numPole );
		std::vector<Real> absWeightRho( numPole );
		for( Int i = 0; i < sortIdx.size(); i++ ){
			sortIdx[i]      = i;
			absWeightRho[i] = std::abs( zweightRho_[i] );
		}
		// Sort in DESCENDING order
		std::sort( sortIdx.begin(), sortIdx.end(), 
				IndexComp<std::vector<Real>& >( absWeightRho ) ) ;
		std::reverse( sortIdx.begin(), sortIdx.end() );

		std::vector<Complex>  zshiftSort( numPole );
		std::vector<Complex>  zweightRhoSort( numPole );
		for( Int i = 0; i < numPole; i++ ){
			zshiftSort[i]     = zshift_[sortIdx[i]];
			zweightRhoSort[i] = zweightRho_[sortIdx[i]];
		}
		zshift_     = zshiftSort;
		zweightRho_ = zweightRhoSort;

		if( isFreeEnergyDensityMatrix ){
			std::vector<Complex>  zweightHelmholtzSort( numPole );
			for( Int i = 0; i < numPole; i++ ){
				zweightHelmholtzSort[i] = zweightHelmholtz_[sortIdx[i]];
			}
			zweightHelmholtz_ = zweightHelmholtzSort;
		}

		if( isEnergyDensityMatrix ){
			std::vector<Complex>  zweightForceSort( numPole );
			for( Int i = 0; i < numPole; i++ ){
				zweightForceSort[i] = zweightForce_[sortIdx[i]];
			}
			zweightForce_ = zweightForceSort;
		}

#if ( _DEBUGlevel_ >= 0 )
		statusOFS << "sorted indicies (according to |weightRho|) " << std::endl << sortIdx << std::endl;
		statusOFS << "sorted zshift" << std::endl << zshift_ << std::endl;
		statusOFS << "sorted zweightRho" << std::endl << zweightRho_ << std::endl;
		if( isFreeEnergyDensityMatrix )
			statusOFS << "sorted zweightHelmholtz" << std::endl << zweightHelmholtz_ << std::endl;
		if( isEnergyDensityMatrix )
			statusOFS << "sorted zweightForce" << std::endl << zweightForce_ << std::endl;
#endif

		// for each pole, perform LDLT factoriation and selected inversion
		Real timeSta, timeEnd;
		Real timePoleSta, timePoleEnd;

		Int numPoleComputed = 0;
		for(Int l = 0; l < numPole; l++){
			if( MYROW( gridPole_ ) == PROW( l, gridPole_ ) ){

				GetTime( timePoleSta );
				statusOFS << "Pole " << l << " processing..." << std::endl;
#if ( _DEBUGlevel_ >= 0 )
				statusOFS << "zshift           = " << zshift_[l] << std::endl;
				statusOFS	<< "zweightRho       = " << zweightRho_[l] << std::endl;
				if( isFreeEnergyDensityMatrix )
					statusOFS << "zweightHelmholtz = " << zweightHelmholtz_[l] << std::endl;
				if( isEnergyDensityMatrix )
					statusOFS << "zweightForce     = " << zweightForce_[l] << std::endl;
#endif
				if( std::abs( zweightRho_[l] ) < poleTolerance ){
					statusOFS << "|weightRho| < poleTolerance, pass the computation of this pole" << std::endl;
				}
				else
				{
					numPoleComputed++;


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

					Int nnzLocal = PMloc.NnzLocal();
					statusOFS << "Number of local nonzeros (L+U) = " << nnzLocal << std::endl;
					Int nnz      = PMloc.Nnz();
					statusOFS << "Number of nonzeros (L+U)       = " << nnz << std::endl;


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

					// Update the density matrix. The following lines are equivalent to
					//
					//				for( Int i = 0; i < rhoMat.nnzLocal; i++ ){
					//					rhoMat.nzvalLocal(i) += 
					//						zweightRho_[l].real() * AinvMat.nzvalLocal(i).imag() + 
					//						zweightRho_[l].imag() * AinvMat.nzvalLocal(i).real();
					//				}
					// 
					// But done more cache-efficiently with blas.
					Real* AinvMatRealPtr = (Real*)AinvMat.nzvalLocal.Data();
					Real* AinvMatImagPtr = AinvMatRealPtr + 1;
					blas::Axpy( rhoMat.nnzLocal, zweightRho_[l].real(), AinvMatImagPtr, 2, 
							rhoMat.nzvalLocal.Data(), 1 );
					blas::Axpy( rhoMat.nnzLocal, zweightRho_[l].imag(), AinvMatRealPtr, 2,
							rhoMat.nzvalLocal.Data(), 1 );

					if( isFreeEnergyDensityMatrix ){
						blas::Axpy( hmzMat.nnzLocal, zweightHelmholtz_[l].real(), AinvMatImagPtr, 2,
								hmzMat.nzvalLocal.Data(), 1 );
						blas::Axpy( hmzMat.nnzLocal, zweightHelmholtz_[l].imag(), AinvMatRealPtr, 2,
								hmzMat.nzvalLocal.Data(), 1 );
					}

					if( isEnergyDensityMatrix ){
						blas::Axpy( frcMat.nnzLocal, zweightForce_[l].real(), AinvMatImagPtr, 2,
								frcMat.nzvalLocal.Data(), 1 );
						blas::Axpy( frcMat.nnzLocal, zweightForce_[l].imag(), AinvMatRealPtr, 2, 
								frcMat.nzvalLocal.Data(), 1 );
					}

					// Update the free energy density matrix and energy density matrix similarly
					GetTime( timePostProcessingEnd );

					statusOFS << "Time for postprocessing is " <<
						timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
				}

				// Compute the contribution to the number of electrons for each pole
#if ( _DEBUGlevel_ >= 0 )
				Real numElecLocal = 
					numElecLocal = blas::Dot( SMat.nnzLocal, SMat.nzvalLocal.Data(),
							1, rhoMat_.nzvalLocal.Data(), 1 );

				Real numElec;
				mpi::Allreduce( &numElecLocal, &numElec, 1, MPI_SUM, gridPole_->comm ); 

				statusOFS << std::endl << "numElecLocal = " << numElecLocal << std::endl;
				statusOFS << "numElecTotal = " << numElec << std::endl << std::endl;
#endif

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

		if( isFreeEnergyDensityMatrix ){
			DblNumVec nzvalHmzMatLocal = hmzMat.nzvalLocal;
			SetValue( hmzMat.nzvalLocal, 0.0 );

			mpi::Allreduce( nzvalHmzMatLocal.Data(), hmzMat.nzvalLocal.Data(),
					hmzMat.nnzLocal, MPI_SUM, gridPole_->colComm );
		}

		if( isEnergyDensityMatrix ){
			DblNumVec nzvalFrcMatLocal = frcMat.nzvalLocal;
			SetValue( frcMat.nzvalLocal, 0.0 );

			mpi::Allreduce( nzvalFrcMatLocal.Data(), frcMat.nzvalLocal.Data(),
					frcMat.nnzLocal, MPI_SUM, gridPole_->colComm );
		}

		// All processors groups compute the number of electrons, and total
		// energy, and optimally the helmholtz free energy

		numElectronNow = CalculateNumElectron( SMat );

    Real totalEnergy = CalculateTotalEnergy( HMat );
		
		Real totalFreeEnergy;
		if( isFreeEnergyDensityMatrix )
			totalFreeEnergy = CalculateFreeEnergy( HMat );

		muList.push_back(muNow);
		numElectronList.push_back( numElectronNow );

		statusOFS << std::endl;
		Print( statusOFS, "Number of poles computed    = ", numPoleComputed );
		Print( statusOFS, "Computed number of electron = ", numElectronNow );
		Print( statusOFS, "Exact number of electron    = ", numElectronExact );
		Print( statusOFS, "Total energy                = ", totalEnergy );
		if( isFreeEnergyDensityMatrix )
			Print( statusOFS, "Total free energy           = ", totalFreeEnergy );

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
	statusOFS << std::endl << "numElecLocal = " << numElecLocal << std::endl;
#endif

	mpi::Allreduce( &numElecLocal, &numElec, 1, MPI_SUM, rhoMat_.comm ); 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return numElec;
} 		// -----  end of method PPEXSIData::CalculateNumElectron  ----- 


Real
PPEXSIData::CalculateTotalEnergy	( const DistSparseMatrix<Real>& HMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateTotalEnergy");
#endif
	
	Real totalEnergyLocal = 0.0, totalEnergy = 0.0;
	
	// TODO Check HMat and rhoMat has the same sparsity

	totalEnergyLocal = blas::Dot( HMat.nnzLocal, HMat.nzvalLocal.Data(),
			1, rhoMat_.nzvalLocal.Data(), 1 );
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << std::endl << "TotalEnergyLocal = " << totalEnergyLocal << std::endl;
#endif

	mpi::Allreduce( &totalEnergyLocal, &totalEnergy, 1, MPI_SUM, rhoMat_.comm ); 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return totalEnergy;
} 		// -----  end of method PPEXSIData::CalculateTotalEnergy  ----- 

Real 
PPEXSIData::CalculateFreeEnergy	( const DistSparseMatrix<Real>& HMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateFreeEnergy");
#endif
	
	Real totalFreeEnergyLocal = 0.0, totalFreeEnergy = 0.0;
	
	// TODO Check HMat and freeEnergyDensityMat_ has the same sparsity

	totalFreeEnergyLocal = blas::Dot( HMat.nnzLocal, HMat.nzvalLocal.Data(),
			1, freeEnergyDensityMat_.nzvalLocal.Data(), 1 );
#if ( _DEBUGlevel_ >= 0 )
	statusOFS << std::endl << "TotalFreeEnergyLocal = " << totalFreeEnergyLocal << std::endl;
#endif

	mpi::Allreduce( &totalFreeEnergyLocal, &totalFreeEnergy, 1, MPI_SUM, 
			freeEnergyDensityMat_.comm ); 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return totalFreeEnergy;
} 		// -----  end of method PPEXSIData::CalculateFreeEnergy  ----- 


Real
PPEXSIData::CalculateForce	( 
		const DistSparseMatrix<Real>& HDerivativeMat,  
		const DistSparseMatrix<Real>& SDerivativeMat )
{
#ifndef _RELEASE_
	PushCallStack("PPEXSIData::CalculateForce");
#endif

	Real totalForceLocal = 0.0, totalForce = 0.0;

	// TODO Check HDerivativeMat, SDerivativeMat, rhoMat and
	// energyDensityMat_ has the same sparsity pattern

	totalForceLocal = 
		- blas::Dot( HDerivativeMat.nnzLocal, HDerivativeMat.nzvalLocal.Data(),
				1, rhoMat_.nzvalLocal.Data(), 1 )
		+ blas::Dot( SDerivativeMat.nnzLocal, SDerivativeMat.nzvalLocal.Data(),
				1, energyDensityMat_.nzvalLocal.Data(), 1 );

#if ( _DEBUGlevel_ >= 0 )
	statusOFS << std::endl << "TotalForceLocal = " << totalForceLocal << std::endl;
#endif

	mpi::Allreduce( &totalForceLocal, &totalForce, 1, MPI_SUM, rhoMat_.comm );

#ifndef _RELEASE_
	PopCallStack();
#endif

	return totalForce;
} 		// -----  end of method PPEXSIData::CalculateForce  ----- 

} //  namespace PEXSI
