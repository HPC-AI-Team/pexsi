/*
	 Copyright (c) 2012 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Authors: Lexing Ying and Lin Lin
	 
   This file is part of PEXSI. All rights reserved.

	 Redistribution and use in source and binary forms, with or without
	 modification, are permitted provided that the following conditions are met:

	 (1) Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.
	 (2) Redistributions in binary form must reproduce the above copyright notice,
	 this list of conditions and the following disclaimer in the documentation
	 and/or other materials provided with the distribution.
	 (3) Neither the name of the University of California, Lawrence Berkeley
	 National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
	 be used to endorse or promote products derived from this software without
	 specific prior written permission.

	 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
	 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	 DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
	 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
	 ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	 You are under no obligation whatsoever to provide any bug fixes, patches, or
	 upgrades to the features, functionality or performance of the source code
	 ("Enhancements") to anyone; however, if you choose to make your Enhancements
	 available either publicly, or directly to Lawrence Berkeley National
	 Laboratory, without imposing a separate written license agreement for such
	 Enhancements, then you hereby grant the following license: a non-exclusive,
	 royalty-free perpetual license to install, use, modify, prepare derivative
	 works, incorporate into other computer software, distribute, and sublicense
	 such enhancements or derivative works thereof, in binary and source code form.
*/
/// @file nummat_impl.hpp
/// @brief Implementation of numerical matrix.
/// @date 2010-09-27
#ifndef _NUMMAT_IMPL_HPP_
#define _NUMMAT_IMPL_HPP_

#include  "nummat_decl.hpp"

namespace  PEXSI{

// TODO Move the things from decl to impl


template <class F> inline void SetValue(NumMat<F>& M, F val)
{
	F *ptr = M.data_;
	for (Int i=0; i < M.m()*M.n(); i++) *(ptr++) = val;
}

template <class F> inline Real Energy(const NumMat<F>& M)
{
  Real sum = 0;
	F *ptr = M.data_;
	for (Int i=0; i < M.m()*M.n(); i++) 
		sum += abs(ptr[i]) * abs(ptr[i]);
  return sum;
}


template <class F> inline void
Transpose ( const NumMat<F>& A, NumMat<F>& B )
{
#ifndef _RELEASE_
	PushCallStack("Transpose");
#endif
	if( A.m() != B.n() || A.n() != B.m() ){
		B.Resize( A.n(), A.m() );
	}

	F* Adata = A.Data();
	F* Bdata = B.Data();
	Int m = A.m(), n = A.n();

	for( Int i = 0; i < m; i++ ){
		for( Int j = 0; j < n; j++ ){
			Bdata[ j + n*i ] = Adata[ i + j*m ];
		}
	}

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Transpose  ----- 

template <class F> inline void
Symmetrize( NumMat<F>& A )
{
#ifndef _RELEASE_
	PushCallStack("Symmetrize");
#endif
	if( A.m() != A.n() ){
		throw std::logic_error( "The matrix to be symmetrized should be a square matrix." );
	}

	NumMat<F> B;
	Transpose( A, B );

	F* Adata = A.Data();
	F* Bdata = B.Data();

	F  half = (F) 0.5;

	for( Int i = 0; i < A.m() * A.n(); i++ ){
		*Adata = half * (*Adata + *Bdata);
		Adata++; Bdata++;
	}

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
}		// -----  end of function Symmetrize ----- 


} // namespace PEXSI

#endif // _NUMMAT_IMPL_HPP_
