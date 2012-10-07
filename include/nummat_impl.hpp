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

} // namespace PEXSI

#endif // _NUMMAT_IMPL_HPP_
