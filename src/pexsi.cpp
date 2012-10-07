#include "pepsi.hpp"

//LLIN: Default parameter setup 
PEXSI::PEXSI(){
}

PEXSI::~PEXSI(){
}

// Use symmetry, compute the trace of the product of two matrices
// sharing the same structure as H
double PEXSI::traceprod(int Ndof, IntNumVec& colptr_H, IntNumVec& rowind_H, 
		     DblNumVec& nzval1, DblNumVec& nzval2){
  double val = 0.0;
  for(int j = 0; j < Ndof; j++){
    for(int i = colptr_H[j]-1; i < colptr_H[j+1]-1; i++){
      if( j+1 == rowind_H[i] ){ // diagonal
	val += nzval1[i]*nzval2[i];  
      }
      else{
	val += 2.0 * nzval1[i]*nzval2[i];
      }
    }
  }
  return val;
}

int PEXSI::setup(){
  cpx czero = cpx(0.0, 0.0);
  /* Make sure that all the input variables have been give outside FIXME */

  iA( _Npole % 2 == 0 ); // Even number of poles

  _mu.resize(_maxit_mu);  _mu[0] = _mu0;
  _Ne.resize(_maxit_mu);

  _nzval_rho.resize(_Hnnz); setvalue(_nzval_rho, 0.0);

  _egy = 0.0;
  _freeegy = 0.0;

  _zshift_rho.resize(_Npole);  setvalue(_zshift_rho, czero);
  _zweight_rho.resize(_Npole); setvalue(_zweight_rho, czero);
  
  if( _hmzflg == 1 ){
    _zshift_hmz.resize(_Npole);  setvalue(_zshift_hmz, czero);
    _zweight_hmz.resize(_Npole); setvalue(_zweight_hmz, czero);
  }

  if( _frcflg == 1 ){
    _zshift_frc.resize(_Npole);  setvalue(_zshift_frc, czero);
    _zweight_frc.resize(_Npole); setvalue(_zweight_frc, czero);
  }

  return 0;
}

// Main subroutine to solve the electronic structure problem 
int PEXSI::solve(){
  
  int token;
  int* perm;
  CpxNumVec nzval_invA, nzval_A;
  IntNumVec colptr_invA, rowind_invA;
  int dumpL = 0;

  token = 0;
  
  _Ne[0] = 0.0;
  _mu[0] = _mu0;

  /* -order     :   Reordering strategy:
       order = -1 (default) : Multiple Minimum Degree Reordering.
       If METIS is supported, then the following choices
       are also available:
       order = 2 : Node Nested Dissection
       order = 3 : Edge Nested Dissection  */

  // Since all the poles share the same algebraic structure, the
  // preprocess can be considered as a once-for-all process,The
  // preprocessing procedure is a once-for-all calculation 
  
  perm = NULL;   // perm is not used here
  SELINV_Interface::ldlt_preprocess__(&token, &_Ndof, 
				      _colptr_H.data(), _rowind_H.data(), 
				      &_Lnnz, &_order, perm);   
  cerr << "Ndof = " << _Ndof << endl;
  cerr << "Hnnz = " << _Hnnz << endl;
  cerr << "Lnnz = " << _Lnnz << endl;

  colptr_invA.resize(_Ndof+1);
  rowind_invA.resize(_Lnnz);
  nzval_invA.resize(_Lnnz);
  nzval_A.resize(_Hnnz);


  for(int iter = 0; iter < _maxit_mu; iter++){
    // Reinitialize the variables
    setvalue(_nzval_rho, 0.0);

    //Initialize the pole expansion
    iC( getpole_rho(reinterpret_cast<doublecomplex*>(_zshift_rho.data()),
		    reinterpret_cast<doublecomplex*>(_zweight_rho.data()),
		    &_Npole, &_temp, &_gap, &_DeltaE, &_mu[iter]) ); 

    // for each pole, perform LDLT factoriation and selected inversion

    for(int l = 0; l < _Npole; l++){
      cerr << _zshift_rho[l] << endl;
      for(int i = 0; i < _Hnnz; i++){
	nzval_A[i] = _nzval_H[i] - _zshift_rho[l] * _nzval_S[i];
      }
      
      SELINV_Interface::ldlt_fact__(&token, _colptr_H.data(),
				    _rowind_H.data(), 
				    reinterpret_cast<doublecomplex*>(nzval_A.data()));

      cerr << "Factorization done" << endl;

      SELINV_Interface::ldlt_blkselinv__(&token, colptr_invA.data(),
					 rowind_invA.data(),
					 reinterpret_cast<doublecomplex*>(nzval_invA.data()),
					 &dumpL);

      cerr << "Selected inversion done" << endl;


      // Evaluate the electron density
      for(int j = 1; j < _Ndof+1; j++){
	for(int ii = _colptr_H[j-1]; ii < _colptr_H[j]; ii++){
	  int k;
	  for(k = colptr_invA[j-1]; k < colptr_invA[j]; k++){
	    if( _rowind_H[ii-1] == rowind_invA[k-1] ){
	      _nzval_rho[ii-1] += _zweight_rho[l].real() * nzval_invA[k-1].imag() + 
		_zweight_rho[l].imag() * nzval_invA[k-1].real();
	      break;
	    }
	  }
	  iA( k != colptr_invA[j] );
	}
      }

    } // for(l)

    // Reduce Ne
    _Ne[iter] = this->traceprod(_Ndof, _colptr_H, _rowind_H, _nzval_S, _nzval_rho);
    cerr << "Ne["<< iter << "] = " << _Ne[iter] << endl;

     
    // Reduce band energy
    // Reduce Helmholtz free energy
    // Reduce force
//    _Ne[iter] = 0.0;
//
//    if( abs(_Ne[iter] - _Neexact) < _tol_Ne ) break;

    // update_mu(iter, _maxit_mu, mu, Ne);
  }

  if (_order == 0) free(perm);
  SELINV_Interface::ldlt_free__(&token); 

  return 0;
}
