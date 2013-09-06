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
/// @file numtns_impl.hpp
/// @brief Implementation of numerical tensor.
/// @date 2010-09-27
#ifndef _NUMTNS_IMPL_HPP_
#define _NUMTNS_IMPL_HPP_

#include  "numtns_decl.hpp"

namespace  PEXSI{

// TODO Move the things from decl to impl


template <class F> inline void SetValue(NumTns<F>& T, F val)
{
	F *ptr = T.data_;
  for(Int i=0; i < T.m() * T.n() * T.p(); i++) *(ptr++) = val; 

	return;
}



template <class F> inline Real Energy(const NumTns<F>& T)
{
  Real sum = 0;

	F *ptr = T.Data();
  for(Int i=0; i < T.m() * T.n() * T.p(); i++) 
		sum += abs(ptr[i]) * abs(ptr[i]);

	return sum;
}

} // namespace PEXSI

#endif // _NUMTNS_IMPL_HPP_
