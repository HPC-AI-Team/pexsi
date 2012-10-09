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

	zshift.Resize( numPole );
	SetValue( zshift, Z_ZERO );
	zweightRho.Resize( numPole );
	SetValue( zweightRho, Z_ZERO );

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

	double muNow = mu0;

	for(int iter = 0; iter < muMaxIter; iter++){
		{
			std::ostringstream msg;
			msg << "Iteration " << iter << ", mu = " << muNow;
			PrintBlock( statusOFS, msg.str() );
		}
		// Reinitialize the variables
		SetValue( rhoMat.nzval, 0.0 );

		//Initialize the pole expansion
		getpole_rho(reinterpret_cast<doublecomplex*>(zshift.Data()),
				reinterpret_cast<doublecomplex*>(zweightRho.Data()),
				&numPole, &temperature, &gap, &deltaE, &muNow); 

		Print( statusOFS, "zshift" );
		statusOFS << zshift << std::endl;
		Print( statusOFS, "zweightRho " );
		statusOFS << zweightRho << std::endl;

		// for each pole, perform LDLT factoriation and selected inversion

		Real timeSta, timeEnd;

		for(Int l = 0; l < numPole; l++){
			statusOFS << "Pole " << l << std::endl;
			statusOFS << "zshift = " << zshift(l) << ", " 
				<< "zweightRho = " << zweightRho(l) << std::endl;

			for(Int i = 0; i < HMat.nnz; i++){
				AMat.nzval(i) = HMat.nzval(i) - zshift(l) * SMat.nzval(i);
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
				Int* colptrHPtr = HMat.colptr.Data();
				Int* rowindHPtr = HMat.rowind.Data();
				Int* colptrInvAPtr = invAMat.colptr.Data();
				Int* rowindInvAPtr = invAMat.rowind.Data();
				Real* nzvalRhoPtr  = rhoMat.nzval.Data();
				Complex* nzvalInvAPtr = invAMat.nzval.Data();
				Complex  ztmp = zweightRho(l);

				for(Int j = 1; j < HMat.size+1; j++){
					statusOFS << j << std::endl;
					for(Int ii = colptrHPtr[j-1]; ii < colptrHPtr[j]; ii++){
//					for(Int ii = HMat.colptr(j-1); ii < HMat.colptr(j); ii++){
						Int kk;
						for(kk = colptrInvAPtr[j-1]; kk < colptrInvAPtr[j]; kk++){
//						for(kk = invAMat.colptr(j-1); kk < invAMat.colptr(j); kk++){
//							if( HMat.rowind(ii-1) == invAMat.rowind(kk-1) ){
							if( rowindHPtr[ii-1] == rowindInvAPtr[kk-1] ){
								nzvalRhoPtr[ii-1] += 
									ztmp.real() * nzvalInvAPtr[kk-1].imag() +
									ztmp.imag() * nzvalInvAPtr[kk-1].real();
//								rhoMat.nzval(ii-1) += 
//									zweightRho(l).real() * invAMat.nzval(kk-1).imag() +
//									zweightRho(l).imag() * invAMat.nzval(kk-1).real();
								break;
							}
						}
						if( kk == colptrInvAPtr[j] ){
//						if( kk == invAMat.colptr(j) ){
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

		} // for(l)

		// Reduce Ne
//		_Ne[iter] = this->traceprod(_Ndof, _colptr_H, _rowind_H, _nzval_S, _nzval_rho);
//		cerr << "Ne["<< iter << "] = " << _Ne[iter] << endl;


		// Reduce band energy
		// Reduce Helmholtz free energy
		// Reduce force
		//    _Ne[iter] = 0.0;
		//
		//    if( abs(_Ne[iter] - _Neexact) < _tol_Ne ) break;

		// update_mu(iter, _maxit_mu, mu, Ne);
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
//double PEXSI::traceprod(int Ndof, IntNumVec& colptr_H, IntNumVec& rowind_H, 
//		DblNumVec& nzval1, DblNumVec& nzval2){
//	double val = 0.0;
//	for(int j = 0; j < Ndof; j++){
//		for(int i = colptr_H[j]-1; i < colptr_H[j+1]-1; i++){
//			if( j+1 == rowind_H[i] ){ // diagonal
//				val += nzval1[i]*nzval2[i];  
//			}
//			else{
//				val += 2.0 * nzval1[i]*nzval2[i];
//			}
//		}
//	}
//	return val;
//}

} //  namespace PEXSI
