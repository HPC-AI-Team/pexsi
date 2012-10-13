#include "pexsi.hpp"

namespace PEXSI{

void
PEXSIData::Setup	(  )
{
#ifndef _RELEASE_
	PushCallStack("PEXSIData::Setup");
#endif
	if( numPole % 2 != 0 ){
		throw std::logic_error( "Must be even number of poles!" );
	}

	muList.clear();
	numElectronList.clear();

	rhoMat.size = HMat.size;
	rhoMat.nnz  = HMat.nnz;
	rhoMat.colptr = IntNumVec( HMat.colptr.m(), false, HMat.colptr.Data() );
	rhoMat.rowind = IntNumVec( HMat.rowind.m(), false, HMat.rowind.Data() );
	rhoMat.nzval.Resize( HMat.nnz );
 	SetValue( rhoMat.nzval, 0.0 );

	totalEnergy = 0.0;
	totalFreeEnergy = 0.0;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PEXSIData::Setup  ----- 


// Main subroutine for the electronic structure calculation
void PEXSIData::Solve( )
{
#ifndef _RELEASE_
	PushCallStack("PEXSIData::Solve");
#endif
	int token = 0;                                
	int dumpL = 0;
	int Lnnz;

	int *perm;                             // Permutation vector

	// Since all the poles share the same algebraic structure, the
	// preprocess can be considered as a once-for-all process,The
	// preprocessing procedure is a once-for-all calculation 

	perm = NULL;   // perm is not used here unless order  0
	SelInvInterface::ldlt_preprocess__(&token, &HMat.size, 
			HMat.colptr.Data(), HMat.rowind.Data(), 
			&Lnnz, &permOrder, perm);   
	Print( statusOFS, "After preprocessing" );
	Print( statusOFS, "Ndof = ", HMat.size );
	Print( statusOFS, "number of nonzeros in H = ", HMat.nnz );
	Print( statusOFS, "number of nonzeros in L = ", Lnnz );

	SparseMatrix<Complex>     AMat;               // A = H - z S
	SparseMatrix<Complex>     invAMat;            // inverse(A)

	AMat.size = HMat.size;
	AMat.nnz  = HMat.nnz;
	// Save some memory
	AMat.colptr = IntNumVec( HMat.colptr.m(), false, HMat.colptr.Data() );
	AMat.rowind = IntNumVec( HMat.rowind.m(), false, HMat.rowind.Data() );
	AMat.nzval.Resize( HMat.nnz );

	invAMat.size = HMat.size;
	invAMat.nnz  = Lnnz;
	invAMat.colptr.Resize( invAMat.size + 1 );
	invAMat.rowind.Resize( Lnnz );
	invAMat.nzval.Resize( Lnnz );

	muList.clear();
	numElectronList.clear();
	Real muNow = mu0;
	Real numElectronNow;

	Real timeMuSta, timeMuEnd;

	for(Int iter = 0; iter < muMaxIter; iter++){
		GetTime( timeMuSta );
		{
			std::ostringstream msg;
			msg << "Iteration " << iter << ", mu = " << muNow;
			PrintBlock( statusOFS, msg.str() );
		}
		// Reinitialize the variables
		SetValue( rhoMat.nzval, 0.0 );

		// Initialize the number of electrons
		numElectronNow  = 0.0;

		//Initialize the pole expansion

		std::vector<Complex> zshiftRaw( numPole );
		std::vector<Complex> zweightRhoRaw( numPole );

		getpole_rho(reinterpret_cast<doublecomplex*>(&zshiftRaw[0]),
				reinterpret_cast<doublecomplex*>(&zweightRhoRaw[0]),
				&numPole, &temperature, &gap, &deltaE, &muNow); 

		// Sort and truncate the poles according to the weights
		{
			std::vector<std::pair<Real,Int> >  weightAbs( numPole );
			for( Int i = 0; i < numPole; i++ ){
				weightAbs[i] = std::pair<Real,Int>(abs( zweightRhoRaw[i] ), i);
			}
			std::sort( weightAbs.begin(), weightAbs.end(), PairGtComparator );

			zshift.clear();
			zweightRho.clear();

			numPoleUsed = 0;
			for( Int i = 0; i < numPole; i++ ){
				if( weightAbs[i].first > poleTolerance ){
					zshift.push_back( zshiftRaw[weightAbs[i].second] );
					zweightRho.push_back( zweightRhoRaw[weightAbs[i].second] );
					numPoleUsed++;
				}
			}

		}

		
		Print( statusOFS, "Number of poles used = ", numPoleUsed );
		Print( statusOFS, "zshift" );
		statusOFS << zshift << std::endl;
		Print( statusOFS, "zweightRho " );
		statusOFS << zweightRho << std::endl;

		// for each pole, perform LDLT factoriation and selected inversion

		Real timeSta, timeEnd;

		for(Int l = 0; l < numPoleUsed; l++){
			statusOFS << "Pole " << l << std::endl;
			statusOFS << "zshift = " << zshift[l] << ", " 
				<< "zweightRho = " << zweightRho[l] << std::endl;

			for(Int i = 0; i < HMat.nnz; i++){
				AMat.nzval(i) = HMat.nzval(i) - zshift[l] * SMat.nzval(i);
			}


			GetTime( timeSta );

			SelInvInterface::ldlt_fact__(&token, HMat.colptr.Data(),
					HMat.rowind.Data(), 
					reinterpret_cast<doublecomplex*>(AMat.nzval.Data()));

			GetTime( timeEnd );

			Print( statusOFS, "Factorization done" );
			Print( statusOFS, "Factorization time = ", timeEnd - timeSta, "[s]" );
			
			GetTime( timeSta );

			SelInvInterface::ldlt_blkselinv__(&token, invAMat.colptr.Data(),
					invAMat.rowind.Data(), 
					reinterpret_cast<doublecomplex*>(invAMat.nzval.Data()), &dumpL);
			
			GetTime( timeEnd );

			Print( statusOFS, "Selected inversion done" );
			Print( statusOFS, "Selected inversion time = ", timeEnd - timeSta, "[s]" );

			// Evaluate the electron density
			GetTime( timeSta );
			
			{
				// Otherwise the speed is too slow.
				Int* colptrHPtr = HMat.colptr.Data();
				Int* rowindHPtr = HMat.rowind.Data();
				Int* colptrInvAPtr = invAMat.colptr.Data();
				Int* rowindInvAPtr = invAMat.rowind.Data();
				Real* nzvalRhoPtr  = rhoMat.nzval.Data();
				Complex* nzvalInvAPtr = invAMat.nzval.Data();
				Complex  zweightl = zweightRho[l];

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
			
			GetTime( timeEnd );

			Print( statusOFS, "Evaluating density done" );
			Print( statusOFS, "Evaluating density time = ", timeEnd - timeSta, "[s]" );


			// Accumulate the number of electrons
			GetTime( timeSta );
			numElectronNow = ProductTrace( SMat.nzval, rhoMat.nzval );
			GetTime( timeEnd );

			Print( statusOFS, "Accumulated number of electron = ", numElectronNow );
			Print( statusOFS, "Evaluating number of electrons done" );
			Print( statusOFS, "Evaluating number of electrons time = ", timeEnd - timeSta, "[s]" );


		} // for(l)

		// Reduce Ne
		
		muList.push_back(muNow);
		numElectronList.push_back( numElectronNow );
		
		Print( statusOFS, "Computed number of electron = ", numElectronNow );
		Print( statusOFS, "Exact number of electron    = ", numElectronExact );


		// Reduce band energy
		// Reduce Helmholtz free energy
		// Reduce force
		if( std::abs( numElectronExact - numElectronList[iter] ) <
				numElectronTolerance ){
			break;
		}

		muNow = UpdateChemicalPotential( iter );

		GetTime( timeMuEnd );

		Print( statusOFS, "Total wall clock time for this iteration = ", 
				timeMuEnd - timeMuSta, "[s]" );
	}

	// FIXME
	if ( permOrder == 0) free(perm);
	SelInvInterface::ldlt_free__(&token); 

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method PEXSIData::Solve----- 


  
// Use symmetry, compute the trace of the product of two matrices
// sharing the same structure as H

Real PEXSIData::ProductTrace	( 
		const DblNumVec& nzval1, const DblNumVec& nzval2 )
{
#ifndef _RELEASE_
	PushCallStack("PEXSIData::ProductTrace");
#endif
	Real val = 0.0;
	for(Int j = 0; j < HMat.size; j++){
		for(Int i = HMat.colptr(j) - 1; i < HMat.colptr(j+1) - 1; i++){
			if( j+1 == HMat.rowind(i) ){ // diagonal
				val += nzval1(i)*nzval2(i);  
			}
			else{
				val += 2.0 * nzval1(i)*nzval2(i);
			}
		}
	}

#ifndef _RELEASE_
	PopCallStack();
#endif

	return val;
} 		// -----  end of method PEXSIData::ProductTrace  ----- 


Real PEXSIData::UpdateChemicalPotential	( const Int iter )
{
#ifndef _RELEASE_
	PushCallStack("PEXSIData::UpdateChemicalPotential");
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
} 		// -----  end of method PEXSIData::UpdateChemicalPotential  ----- 

} //  namespace PEXSI
