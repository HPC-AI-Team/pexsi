/*
	 Copyright (c) 2012 The Regents of the University of California,
	 through Lawrence Berkeley National Laboratory.  

   Authors: Lexing Ying, Mathias Jacquelin and Lin Lin
	 
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
/// @file tinyvec_impl.hpp
/// @brief Implementation of tiny vectors
/// @date 2010-09-20
#ifndef  _PEXSI_TINYVEC_IMPL_HPP_
#define  _PEXSI_TINYVEC_IMPL_HPP_


namespace PEXSI{

	// Template for tiny vectors of dimension 3.

	template <class F> 
		inline F&
		Vec3T<F>::operator() ( Int i ) 
		{
#ifndef _RELEASE_
			PushCallStack("Vec3T::operator()");
#endif  // ifndef _RELEASE_
			if( i < 0 || i > 2 ){
				#ifdef USE_ABORT
abort();
#endif
ErrorHandling( "Index is out of bound." );
			}
#ifndef _RELEASE_
			PopCallStack();
#endif  // ifndef _RELEASE_
			return v_[i];
		} 		// -----  end of method Vec3T::operator()  ----- 

	template <class F> 
		inline const F&
		Vec3T<F>::operator() ( Int i ) const
		{
#ifndef _RELEASE_
			PushCallStack("Vec3T::operator()");
#endif  // ifndef _RELEASE_
			if( i < 0 || i > 2 ){
				#ifdef USE_ABORT
abort();
#endif
ErrorHandling( "Index is out of bound." );
			}
#ifndef _RELEASE_
			PopCallStack();
#endif  // ifndef _RELEASE_
			return v_[i];
		} 		// -----  end of method Vec3T::operator()  ----- 

	template <class F> 
		inline F&
		Vec3T<F>::operator[] ( Int i ) 
		{
#ifndef _RELEASE_
			PushCallStack("Vec3T::operator[]");
#endif  // ifndef _RELEASE_
			if( i < 0 || i > 2 ){
				#ifdef USE_ABORT
abort();
#endif
ErrorHandling( "Index is out of bound." );
			}
#ifndef _RELEASE_
			PopCallStack();
#endif  // ifndef _RELEASE_
			return v_[i];
		} 		// -----  end of method Vec3T::operator[]  ----- 

	template <class F> 
		inline const F&
		Vec3T<F>::operator[] ( Int i ) const
		{
#ifndef _RELEASE_
			PushCallStack("Vec3T::operator[]");
#endif  // ifndef _RELEASE_
			if( i < 0 || i > 2 ){
				#ifdef USE_ABORT
abort();
#endif
ErrorHandling( "Index is out of bound." );
			}
#ifndef _RELEASE_
			PopCallStack();
#endif  // ifndef _RELEASE_
			return v_[i];
		} 		// -----  end of method Vec3T::operator[]  ----- 


} // namespace PEXSI

#endif // _PEXSI_TINYVEC_IMPL_HPP_
