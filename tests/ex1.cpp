#include "pexsi.hpp"

using namespace PEXSI;

int main(int argc, char **argv) 
{
  PEXSI pexsi;

//  pepsi._order = -1;
//  pepsi._Hnnz = 0;
//  pepsi._Npole = 2;
//  pepsi._temp = 300;
//  pepsi._gap = 0.0;
//  pepsi._DeltaE = 1.0;
//  pepsi._mu0 = -1.0;
//  pepsi._Neexact = 1;
//  pepsi._tol_Ne = 1e-6;
//  pepsi._hmzflg = 0;
//  pepsi._frcflg = 0;
//  pepsi._maxit_mu = 1;
//
//  {
//    char* filename = "Hc.ccf";
//    int Ndof, Hnnz;
//    readmatrixheader_(filename, &Ndof, &Hnnz);
//
//    IntNumVec colptr_H(Ndof+1);
//    IntNumVec rowind_H(Hnnz);
//    CpxNumVec tmp(Hnnz);
//    readcmatrix_(filename, colptr_H.data(), 
//		 rowind_H.data(), 
//		 reinterpret_cast<doublecomplex*>(tmp.data()));
//    pepsi._Hnnz = Hnnz;
//    pepsi._Ndof = Ndof;
//    pepsi._colptr_H = colptr_H;
//    pepsi._rowind_H = rowind_H;
//    pepsi._nzval_H.resize(Hnnz);
//    pepsi._nzval_S.resize(Hnnz);
//    cerr << colptr_H << endl;
//    for(int i = 0; i < Hnnz; i++){
//      pepsi._nzval_H[i] = tmp[i].real();
//    }
//
//    // Orthogonal overlap matrix. Note that CCS format here uses one-based
//    // notation.
//    for(int i = 1; i < Ndof+1; i++){
//      for(int j = colptr_H[i-1]; j < colptr_H[i]; j++){
//	if(i == rowind_H[j-1]){
//	  pepsi._nzval_S[j-1] = 1.0;
//	}
//	else{
//	  pepsi._nzval_S[j-1] = 0.0;
//	}
//      }
//    }
//  }
//
//  iC(pepsi.setup());
//  iC(pepsi.solve());
  
  return 0;
}
