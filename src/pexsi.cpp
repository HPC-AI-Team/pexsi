/// @file pexsi.cpp
/// @brief Implementation of the sequential PEXSI.
/// @author Lin Lin
/// @date 2012-10-15
#include "pexsi.hpp"

namespace PEXSI{


Real PEXSIData::CalculateChemicalPotential	( 
			const Int iter, 
			const Real numElectronExact, 
			const Real numElectronTolerance, 
			const std::vector<Real>& muList,
			const std::vector<Real>& numElectronList )
{
#ifndef _RELEASE_
	PushCallStack("PEXSIData::CalculateChemicalPotential");
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
} 		// -----  end of method PEXSIData::CalculateChemicalPotential  ----- 


// Main subroutine for the electronic structure calculation
void PEXSIData::Solve( 
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
		){
#ifndef _RELEASE_
	PushCallStack("PEXSIData::Solve");
#endif


	// *********************************************************************
	// Check the input parameters
	// *********************************************************************
	if( numPole % 2 != 0 ){
		throw std::logic_error( "Must be even number of poles!" );
	}

	// TODO Check H and S have the same pattern

	// *********************************************************************
	// Initialize
	// *********************************************************************
	muList.clear();
	numElectronList.clear();

	SparseMatrix<Complex>  AMat;              // A = H - z * S
	// rename for convenience
	SparseMatrix<Real>& rhoMat       = rhoMat_;     
//	SparseMatrix<Real>& hmzMat       = freeEnergyDensityMat_;
//	SparseMatrix<Real>& frcMat     = energyDensityMat_;

	// Copy the pattern
	CopyPattern( HMat, AMat );
	CopyPattern( HMat, rhoMat );
//	if( isFreeEnergyDensityMatrix )
//		CopyPattern( HMat, hmzMat );
//	if( isEnergyDensityMatrix )
//		CopyPattern( HMat, frcMat );

	SetValue( AMat.nzval, Z_ZERO );  

	// The number of nonzero factors in L (including diagonal blocks).
	Int Lnnz;
	// Permutation order. Currently only MMD is supported.
	Int permOrder;

	if( ColPerm == "MMD_AT_PLUS_A" ){
		permOrder = -1;
	}
	else{
		std::ostringstream msg;
		msg << ColPerm << " is not a supported ColPerm type. Try (case sensitive) " << std::endl
			  << "MMD_AT_PLUS_A" << std::endl;
		throw std::runtime_error( msg.str().c_str() );
	}


	// *********************************************************************
	// Symbolic factorization.  
	// Each numPoleGroup perform independently
	// *********************************************************************
	SelInvInterface  selinv;
	selinv.SymbolicFactorize( AMat, permOrder, NULL, Lnnz );
	
	Print( statusOFS, "After preprocessing" );
	Print( statusOFS, "Ndof = ", HMat.size );
	Print( statusOFS, "number of nonzeros in H = ", HMat.nnz );
	Print( statusOFS, "number of nonzeros in L = ", Lnnz );
	
	SparseMatrix<Complex>     invAMat;            // inverse(A)
	invAMat.size = HMat.size; 
	invAMat.nnz  = Lnnz;
	invAMat.colptr.Resize( invAMat.size + 1 );
	invAMat.rowind.Resize( Lnnz );
	invAMat.nzval.Resize( Lnnz );

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
		SetValue( rhoMat.nzval, 0.0 );
//		if( isFreeEnergyDensityMatrix )
//			SetValue( hmzMat.nzvalLocal, 0.0 );
//		if( isEnergyDensityMatrix )
//			SetValue( frcMat.nzvalLocal, 0.0 );

		// Initialize the number of electrons
		numElectronNow  = 0.0;

		//Initialize the pole expansion
		zshift_.resize( numPole );
		zweightRho_.resize( numPole );

		GetPoleDensity( &zshift_[0], &zweightRho_[0],
				numPole, temperature, gap, deltaE, muNow ); 

//		if( isFreeEnergyDensityMatrix ){
//			std::vector<Complex>  zshiftTmp( numPole );
//			zweightHelmholtz_.resize( numPole );
//			GetPoleHelmholtz( &zshiftTmp[0], &zweightHelmholtz_[0], 
//				numPole, temperature, gap, deltaE, muNow ); 
//			Real norm = 0.0;
//			for( Int i = 0; i < numPole; i++ ){
//				norm += pow( std::abs( zshiftTmp[i] - zshift_[i] ), 2.0 );
//			}
//			if( norm > 1e-12 )
//				throw std::runtime_error("The pole shifts for rho and for Helmholtz do not match.");
//		}
//
//		if( isEnergyDensityMatrix ){
//			std::vector<Complex>  zshiftTmp( numPole );
//			zweightForce_.resize( numPole );
//			GetPoleForce( &zshiftTmp[0], &zweightForce_[0],
//				numPole, temperature, gap, deltaE, muNow ); 
//			Real norm = 0.0;
//			for( Int i = 0; i < numPole; i++ ){
//				norm += pow( std::abs( zshiftTmp[i] - zshift_[i] ), 2.0 );
//			}
//			if( norm > 1e-12 )
//				throw std::runtime_error("The pole shifts for rho and for energy density matrix do not match.");
//		}


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

//		if( isFreeEnergyDensityMatrix ){
//			std::vector<Complex>  zweightHelmholtzSort( numPole );
//			for( Int i = 0; i < numPole; i++ ){
//				zweightHelmholtzSort[i] = zweightHelmholtz_[sortIdx[i]];
//			}
//			zweightHelmholtz_ = zweightHelmholtzSort;
//		}
//
//		if( isEnergyDensityMatrix ){
//			std::vector<Complex>  zweightForceSort( numPole );
//			for( Int i = 0; i < numPole; i++ ){
//				zweightForceSort[i] = zweightForce_[sortIdx[i]];
//			}
//			zweightForce_ = zweightForceSort;
//		}

#if ( _DEBUGlevel_ >= 0 )
		statusOFS << "sorted indicies (according to |weightRho|) " << std::endl << sortIdx << std::endl;
		statusOFS << "sorted zshift" << std::endl << zshift_ << std::endl;
		statusOFS << "sorted zweightRho" << std::endl << zweightRho_ << std::endl;
//		if( isFreeEnergyDensityMatrix )
//			statusOFS << "sorted zweightHelmholtz" << std::endl << zweightHelmholtz_ << std::endl;
//		if( isEnergyDensityMatrix )
//			statusOFS << "sorted zweightForce" << std::endl << zweightForce_ << std::endl;
#endif

		// for each pole, perform LDLT factoriation and selected inversion
		Real timeSta, timeEnd;
		Real timePoleSta, timePoleEnd;

		Int numPoleComputed = 0;
		for(Int l = 0; l < numPole; l++){

			GetTime( timePoleSta );
			statusOFS << "Pole " << l << " processing..." << std::endl;
#if ( _DEBUGlevel_ >= 0 )
			statusOFS << "zshift           = " << zshift_[l] << std::endl;
			statusOFS	<< "zweightRho       = " << zweightRho_[l] << std::endl;
//			if( isFreeEnergyDensityMatrix )
//				statusOFS << "zweightHelmholtz = " << zweightHelmholtz_[l] << std::endl;
//			if( isEnergyDensityMatrix )
//				statusOFS << "zweightForce     = " << zweightForce_[l] << std::endl;
#endif
			if( std::abs( zweightRho_[l] ) < poleTolerance ){
				statusOFS << "|weightRho| < poleTolerance, pass the computation of this pole" << std::endl;
			}
			else
			{
				numPoleComputed++;

				for( Int i = 0; i < HMat.nnz; i++ ){
					AMat.nzval(i) = HMat.nzval(i) - zshift_[l] * SMat.nzval(i);
				}


				// *********************************************************************
				// Factorization
				// *********************************************************************
				Real timeTotalFactorizationSta, timeTotalFactorizationEnd;

				GetTime( timeTotalFactorizationSta );

				selinv.NumericalFactorize( AMat );

				GetTime( timeTotalFactorizationEnd );

				statusOFS << "Time for total factorization is " << timeTotalFactorizationEnd - timeTotalFactorizationSta<< " [s]" << std::endl; 

				// *********************************************************************
				// Selected inversion
				// *********************************************************************
				Real timeTotalSelInvSta, timeTotalSelInvEnd;
				GetTime( timeTotalSelInvSta );
			
				selinv.SelInv( invAMat );

				GetTime( timeTotalSelInvEnd );

				statusOFS << "Time for total selected inversion is " <<
					timeTotalSelInvEnd  - timeTotalSelInvSta << " [s]" << std::endl;

				// *********************************************************************
				// Postprocessing
				// *********************************************************************

				Real timePostProcessingSta, timePostProcessingEnd;

				GetTime( timePostProcessingSta );

				{
					// Use the pointer format
					Int* colptrHPtr = HMat.colptr.Data();
					Int* rowindHPtr = HMat.rowind.Data();
					Int* colptrInvAPtr = invAMat.colptr.Data();
					Int* rowindInvAPtr = invAMat.rowind.Data();
					Real* nzvalRhoPtr  = rhoMat.nzval.Data();
					Complex* nzvalInvAPtr = invAMat.nzval.Data();
					Complex  zweightl = zweightRho_[l];

					for(Int j = 1; j < HMat.size+1; j++){
						for(Int ii = colptrHPtr[j-1]; ii < colptrHPtr[j]; ii++){
							Int kk;
							for(kk = colptrInvAPtr[j-1]; kk < colptrInvAPtr[j]; kk++){
								if( rowindHPtr[ii-1] == rowindInvAPtr[kk-1] ){
									nzvalRhoPtr[ii-1] += 
										zweightl.real() * nzvalInvAPtr[kk-1].imag() +
										zweightl.imag() * nzvalInvAPtr[kk-1].real();
									break;
								}
							}
							if( kk == colptrInvAPtr[j] ){
								std::ostringstream msg;
								msg << "H(" << HMat.rowind(ii-1) << ", " << j << ") cannot be found in invA" << std::endl;
								throw std::logic_error( msg.str().c_str() );
							}
						}
					} // for (j)
				} 


//				if( isFreeEnergyDensityMatrix ){
//					blas::Axpy( hmzMat.nnzLocal, zweightHelmholtz_[l].real(), AinvMatImagPtr, 2,
//							hmzMat.nzvalLocal.Data(), 1 );
//					blas::Axpy( hmzMat.nnzLocal, zweightHelmholtz_[l].imag(), AinvMatRealPtr, 2,
//							hmzMat.nzvalLocal.Data(), 1 );
//				}

//				if( isEnergyDensityMatrix ){
//					blas::Axpy( frcMat.nnzLocal, zweightForce_[l].real(), AinvMatImagPtr, 2,
//							frcMat.nzvalLocal.Data(), 1 );
//					blas::Axpy( frcMat.nnzLocal, zweightForce_[l].imag(), AinvMatRealPtr, 2, 
//							frcMat.nzvalLocal.Data(), 1 );
//				}

				// Update the free energy density matrix and energy density matrix similarly
				GetTime( timePostProcessingEnd );

				statusOFS << "Time for postprocessing is " <<
					timePostProcessingEnd - timePostProcessingSta << " [s]" << std::endl;
			}

			// Compute the contribution to the number of electrons for each pole
#if ( _DEBUGlevel_ >= 0 )
			numElectronNow = CalculateNumElectron( SMat );

			statusOFS << "numElecTotal = " << numElectronNow << std::endl << std::endl;
#endif

			GetTime( timePoleEnd );

			statusOFS << "Time for pole " << l << " is " <<
				timePoleEnd - timePoleSta << " [s]" << std::endl << std::endl;

		} // for(l)

		// Reduce the density matrix across the processor rows in gridPole_

//		if( isFreeEnergyDensityMatrix ){
//			DblNumVec nzvalHmzMatLocal = hmzMat.nzvalLocal;
//			SetValue( hmzMat.nzvalLocal, 0.0 );
//
//			mpi::Allreduce( nzvalHmzMatLocal.Data(), hmzMat.nzvalLocal.Data(),
//					hmzMat.nnzLocal, MPI_SUM, gridPole_->colComm );
//		}

//		if( isEnergyDensityMatrix ){
//			DblNumVec nzvalFrcMatLocal = frcMat.nzvalLocal;
//			SetValue( frcMat.nzvalLocal, 0.0 );
//
//			mpi::Allreduce( nzvalFrcMatLocal.Data(), frcMat.nzvalLocal.Data(),
//					frcMat.nnzLocal, MPI_SUM, gridPole_->colComm );
//		}

		// All processors groups compute the number of electrons, and total
		// energy, and optimally the helmholtz free energy

		numElectronNow = CalculateNumElectron( SMat );

    Real totalEnergy = CalculateTotalEnergy( HMat );
		
//		Real totalFreeEnergy;
//		if( isFreeEnergyDensityMatrix )
//			totalFreeEnergy = CalculateFreeEnergy( HMat );

		muList.push_back(muNow);
		numElectronList.push_back( numElectronNow );

		statusOFS << std::endl;
		Print( statusOFS, "Number of poles computed    = ", numPoleComputed );
		Print( statusOFS, "Computed number of electron = ", numElectronNow );
		Print( statusOFS, "Exact number of electron    = ", numElectronExact );
		Print( statusOFS, "Total energy                = ", totalEnergy );
//		if( isFreeEnergyDensityMatrix )
//			Print( statusOFS, "Total free energy           = ", totalFreeEnergy );

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
} 		// -----  end of method PEXSIData::Solve----- 

Real
PEXSIData::CalculateNumElectron	( const SparseMatrix<Real>& SMat )
{
#ifndef _RELEASE_
	PushCallStack("PEXSIData::CalculateNumElectron");
#endif
	// TODO Check SMat and rhoMat has the same sparsity

	Real val = 0.0;
	Real *ptr1 = SMat.nzval.Data();
	Real *ptr2 = rhoMat_.nzval.Data();
	for(Int j = 0; j < SMat.size; j++){
		for(Int i = SMat.colptr[j] - 1; i < SMat.colptr[j+1] - 1; i++){
			if( j+1 == SMat.rowind[i] ){ // diagonal
				val += ptr1[i]*ptr2[i];  
			}
			else{
				val += 2.0 * ptr1[i]*ptr2[i];
			}
		}
	}

#ifndef _RELEASE_
	PopCallStack();
#endif

	return val;
} 		// -----  end of method PEXSIData::CalculateNumElectron  ----- 

Real
PEXSIData::CalculateTotalEnergy	( const SparseMatrix<Real>& HMat )
{
#ifndef _RELEASE_
	PushCallStack("PEXSIData::CalculateTotalEnergy");
#endif
	
	// TODO Check HMat and rhoMat has the same sparsity

	Real val = 0.0;
	Real *ptr1 = HMat.nzval.Data();
	Real *ptr2 = rhoMat_.nzval.Data();
	for(Int j = 0; j < HMat.size; j++){
		for(Int i = HMat.colptr[j] - 1; i < HMat.colptr[j+1] - 1; i++){
			if( j+1 == HMat.rowind[i] ){ // diagonal
				val += ptr1[i]*ptr2[i];  
			}
			else{
				val += 2.0 * ptr1[i]*ptr2[i];
			}
		}
	}

#ifndef _RELEASE_
	PopCallStack();
#endif

	return val;
} 		// -----  end of method PEXSIData::CalculateTotalEnergy  ----- 




} //  namespace PEXSI
