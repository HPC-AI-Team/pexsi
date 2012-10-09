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

	/* -order     :   Reordering strategy:
		 order = -1 (default) : Multiple Minimum Degree Reordering.
		 If METIS is supported, then the following choices
		 are also available:
		 order = 2 : Node Nested Dissection
		 order = 3 : Edge Nested Dissection  */

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
	//	FIXME
//	AMat.Resize( HMat.nnz, HMat.size );
//	invAMat.Resize( HMat.nnz, Lnnz );

	double muNow = mu0;
//
//	for(int iter = 0; iter < muMaxIter; iter++){
//		// Reinitialize the variables
//		SetValue( rhoMat.nzval, 0.0 );
//
//		//Initialize the pole expansion
//		getpole_rho(reinterpret_cast<doublecomplex*>(zshift.Data()),
//				reinterpret_cast<doublecomplex*>(zweightRho.Data()),
//				&numPole, &temperature, &gap, &deltaE, &muNow); 
//
//		// for each pole, perform LDLT factoriation and selected inversion
//
//		for(int l = 0; l < numPole; l++){
//			cerr << _zshift_rho[l] << endl;
//			for(int i = 0; i < _Hnnz; i++){
//				nzval_A[i] = _nzval_H[i] - _zshift_rho[l] * _nzval_S[i];
//			}
//
//			SelInvInterface::ldlt_fact__(&token, _colptr_H.data(),
//					_rowind_H.data(), 
//					reinterpret_cast<doublecomplex*>(nzval_A.data()));
//
//			cerr << "Factorization done" << endl;
//
//			SelInvInterface::ldlt_blkselinv__(&token, colptr_invA.data(),
//					rowind_invA.data(),
//					reinterpret_cast<doublecomplex*>(nzval_invA.data()),
//					&dumpL);
//
//			cerr << "Selected inversion done" << endl;
//
//
//			// Evaluate the electron density
//			for(int j = 1; j < _Ndof+1; j++){
//				for(int ii = _colptr_H[j-1]; ii < _colptr_H[j]; ii++){
//					int k;
//					for(k = colptr_invA[j-1]; k < colptr_invA[j]; k++){
//						if( _rowind_H[ii-1] == rowind_invA[k-1] ){
//							_nzval_rho[ii-1] += _zweight_rho[l].real() * nzval_invA[k-1].imag() + 
//								_zweight_rho[l].imag() * nzval_invA[k-1].real();
//							break;
//						}
//					}
//					iA( k != colptr_invA[j] );
//				}
//			}
//
//		} // for(l)
//
//		// Reduce Ne
//		_Ne[iter] = this->traceprod(_Ndof, _colptr_H, _rowind_H, _nzval_S, _nzval_rho);
//		cerr << "Ne["<< iter << "] = " << _Ne[iter] << endl;
//
//
//		// Reduce band energy
//		// Reduce Helmholtz free energy
//		// Reduce force
//		//    _Ne[iter] = 0.0;
//		//
//		//    if( abs(_Ne[iter] - _Neexact) < _tol_Ne ) break;
//
//		// update_mu(iter, _maxit_mu, mu, Ne);
//	}

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
