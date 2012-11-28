#ifndef  _TINYVEC_IMPL_HPP_
#define  _TINYVEC_IMPL_HPP_

#include "tinyvec_decl.hpp"

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
			throw std::logic_error( "Index is out of bound." );
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
			throw std::logic_error( "Index is out of bound." );
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
			throw std::logic_error( "Index is out of bound." );
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
			throw std::logic_error( "Index is out of bound." );
		}
#ifndef _RELEASE_
		PopCallStack();
#endif  // ifndef _RELEASE_
		return v_[i];
	} 		// -----  end of method Vec3T::operator[]  ----- 


} // namespace PEXSI

#endif // _TINYVEC_IMPL_HPP_
