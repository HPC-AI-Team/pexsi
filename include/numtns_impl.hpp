/// @file numtns_impl.hpp
/// @brief Implementation of numerical tensor.
/// @author Lexing Ying and Lin Lin
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
