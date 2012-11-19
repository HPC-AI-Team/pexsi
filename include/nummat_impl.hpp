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


} // namespace PEXSI

#endif // _NUMMAT_IMPL_HPP_
