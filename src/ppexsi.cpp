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

	DistSparseMatrix<Complex>  AMat;             // A = H - z * S
	// Copy the pattern
	CopyPattern( HMat, AMat );
	SetValue( AMat.nzvalLocal, Z_ZERO );          // Symbolic factorization does not need value

	statusOFS << "AMat.nnzLocal = " << AMat.nnzLocal << std::endl;

	SuperLUMatrix              luMat( *gridSuperLU_ );  // SuperLU matrix.  Can be used for sparsity pattern.



	// *********************************************************************
	// Symbolic factorization.  
  // Each numPoleGroup perform independently
	// *********************************************************************
	luMat.DistSparseMatrixToSuperMatrixNRloc( AMat );
	luMat.SymbolicFactorize();
	luMat.DestroyAOnly();


//	Real muNow = mu0;
//	Real numElectronNow;
//
//	Real timeMuSta, timeMuEnd;
//
//	for(Int iter = 0; iter < muMaxIter; iter++){
//		GetTime( timeMuSta );
//		statusOFS << "Iteration " << iter << ", mu = " << muNow << std::endl;
//		// Reinitialize the variables
//		SetValue( rhoMat.nzval, 0.0 );
//
//		// Initialize the number of electrons
//		numElectronNow  = 0.0;
//
//		//Initialize the pole expansion
//
//		std::vector<Complex> zshiftRaw( numPole );
//		std::vector<Complex> zweightRhoRaw( numPole );
//
//		GetPoleDensity(&zshiftRaw[0], &zweightRhoRaw[0],
//				numPole, temperature, gap, deltaE, muNow); 
//
//		// Sort and truncate the poles according to the weights
//		{
//			std::vector<std::pair<Real,Int> >  weightAbs( numPole );
//			for( Int i = 0; i < numPole; i++ ){
//				weightAbs[i] = std::pair<Real,Int>(abs( zweightRhoRaw[i] ), i);
//			}
//			std::sort( weightAbs.begin(), weightAbs.end(), PairGtComparator );
//
//			zshift.clear();
//			zweightRho.clear();
//
//			numPoleUsed = 0;
//			for( Int i = 0; i < numPole; i++ ){
//				if( weightAbs[i].first > poleTolerance ){
//					zshift.push_back( zshiftRaw[weightAbs[i].second] );
//					zweightRho.push_back( zweightRhoRaw[weightAbs[i].second] );
//					numPoleUsed++;
//				}
//			}
//
//		}
//
//		
//		Print( statusOFS, "Number of poles used = ", numPoleUsed );
//		Print( statusOFS, "zshift" );
//		statusOFS << zshift << std::endl;
//		Print( statusOFS, "zweightRho " );
//		statusOFS << zweightRho << std::endl;
//
//		// for each pole, perform LDLT factoriation and selected inversion
//
//		Real timeSta, timeEnd;
//		Real timePoleSta, timePoleEnd;
//
//		for(Int l = 0; l < numPoleUsed; l++){
//			GetTime( timePoleSta );
//			statusOFS << "Pole " << l << std::endl;
//			statusOFS << "zshift = " << zshift[l] << ", " 
//				<< "zweightRho = " << zweightRho[l] << std::endl;
//
//			for(Int i = 0; i < HMat.nnz; i++){
//				AMat.nzval(i) = HMat.nzval(i) - zshift[l] * SMat.nzval(i);
//			}
//
//
//			GetTime( timeSta );
//
//			SelInvInterface::ldlt_fact__(&token, HMat.colptr.Data(),
//					HMat.rowind.Data(), 
//					reinterpret_cast<doublecomplex*>(AMat.nzval.Data()));
//
//			GetTime( timeEnd );
//
//			Print( statusOFS, "Factorization done" );
//			Print( statusOFS, "Factorization time = ", timeEnd - timeSta, "[s]" );
//			
//			GetTime( timeSta );
//
//			SelInvInterface::ldlt_blkselinv__(&token, invAMat.colptr.Data(),
//					invAMat.rowind.Data(), 
//					reinterpret_cast<doublecomplex*>(invAMat.nzval.Data()), &dumpL);
//			
//			GetTime( timeEnd );
//
//			Print( statusOFS, "Selected inversion done" );
//			Print( statusOFS, "Selected inversion time = ", timeEnd - timeSta, "[s]" );
//
//			// Evaluate the electron density
//			GetTime( timeSta );
//			
//			{
//				// Otherwise the speed is too slow.
//				Int* colptrHPtr = HMat.colptr.Data();
//				Int* rowindHPtr = HMat.rowind.Data();
//				Int* colptrInvAPtr = invAMat.colptr.Data();
//				Int* rowindInvAPtr = invAMat.rowind.Data();
//				Real* nzvalRhoPtr  = rhoMat.nzval.Data();
//				Complex* nzvalInvAPtr = invAMat.nzval.Data();
//				Complex  zweightl = zweightRho[l];
//
//				for(Int j = 1; j < HMat.size+1; j++){
//					for(Int ii = colptrHPtr[j-1]; ii < colptrHPtr[j]; ii++){
//						Int kk;
//						for(kk = colptrInvAPtr[j-1]; kk < colptrInvAPtr[j]; kk++){
//							if( rowindHPtr[ii-1] == rowindInvAPtr[kk-1] ){
//								nzvalRhoPtr[ii-1] += 
//									zweightl.real() * nzvalInvAPtr[kk-1].imag() +
//									zweightl.imag() * nzvalInvAPtr[kk-1].real();
//								break;
//							}
//						}
//						if( kk == colptrInvAPtr[j] ){
//							std::ostringstream msg;
//							msg << "H(" << HMat.rowind(ii-1) << ", " << j << ") cannot be found in invA" << std::endl;
//							throw std::logic_error( msg.str().c_str() );
//						}
//					}
//				} // for (j)
//			} 
//			
//			GetTime( timeEnd );
//
//			Print( statusOFS, "Evaluating density done" );
//			Print( statusOFS, "Evaluating density time = ", timeEnd - timeSta, "[s]" );
//
//
//			// Accumulate the number of electrons
//			GetTime( timeSta );
//			numElectronNow = ProductTrace( SMat.nzval, rhoMat.nzval );
//			GetTime( timeEnd );
//			
//			GetTime( timePoleEnd );
//
//			Print( statusOFS, "Accumulated number of electron = ", numElectronNow );
//			Print( statusOFS, "Evaluating number of electrons done" );
//			Print( statusOFS, "Evaluating number of electrons time = ", timeEnd - timeSta, "[s]" );
//
//			Print( statusOFS, "\nTime for this pole = ", timePoleEnd - timePoleSta, "[s]" );
//
//		} // for(l)
//
//		// Reduce Ne
//		
//		muList.push_back(muNow);
//		numElectronList.push_back( numElectronNow );
//		
//		Print( statusOFS, "Computed number of electron = ", numElectronNow );
//		Print( statusOFS, "Exact number of electron    = ", numElectronExact );
//
//
//		// Reduce band energy
//		// Reduce Helmholtz free energy
//		// Reduce force
//		if( std::abs( numElectronExact - numElectronList[iter] ) <
//				numElectronTolerance ){
//			break;
//		}
//
//		muNow = UpdateChemicalPotential( iter );
//
//		GetTime( timeMuEnd );
//
//		Print( statusOFS, "Total wall clock time for this iteration = ", 
//				timeMuEnd - timeMuSta, "[s]" );
//	}
//
//	// FIXME
//	if ( permOrder == 0) free(perm);
//	SelInvInterface::ldlt_free__(&token); 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PPEXSIData::Solve----- 


  
// Use symmetry, compute the trace of the product of two matrices
// sharing the same structure as H


//Real PPEXSIData::UpdateChemicalPotential	( const Int iter )
//{
//#ifndef _RELEASE_
//	PushCallStack("PPEXSIData::UpdateChemicalPotential");
//#endif
//  // FIXME Magic number here
//	Real  muMin = -5.0, muMax = 5.0, muMinStep = 0.01;;
//	Real  muNew;
//
//	if( iter == 0 ){
//		if( numElectronExact > numElectronList[iter] ){
//			muNew = muList[iter] + muMinStep;
//		}
//		else{
//			muNew = muList[iter] - muMinStep;
//		}
//	}
//	else{
//		if( std::abs(numElectronList[iter] -  numElectronList[iter-1])
//		    < numElectronTolerance ){
//			statusOFS << "The number of electrons did not change." << std::endl;
//			if( numElectronExact > numElectronList[iter] ){
//				muNew = muList[iter] + muMinStep;
//			}
//			else{
//				muNew = muList[iter] - muMinStep;
//			}
//		}
//		else {
//			muNew = muList[iter] + (muList[iter] - muList[iter-1]) / 
//				( numElectronList[iter] - numElectronList[iter-1] ) *
//				( numElectronExact - numElectronList[iter] );
//			if( muNew < muMin || muNew > muMax ){
//				statusOFS << "muNew = " << muNew << " is out of bound ["
//					<< muMin << ", " << muMax << "]" << std::endl;
//				if( numElectronExact > numElectronList[iter] ){
//					muNew = muList[iter] + muMinStep;
//				}
//				else{
//					muNew = muList[iter] - muMinStep;
//				}
//			}
//		} // if ( numElectron changed )
//	} // if (iter == 0)
//
//#ifndef _RELEASE_
//	PopCallStack();
//#endif
//
//	return muNew;
//} 		// -----  end of method PPEXSIData::UpdateChemicalPotential  ----- 



} //  namespace PEXSI
