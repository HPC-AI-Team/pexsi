#ifndef _NUMMAT_DECL_HPP_
#define _NUMMAT_DECL_HPP_

#include "environment_impl.hpp"

// TODO Move the things from decl to impl

namespace  PEXSI{

// Templated form of numerical matrix
//
// The main advantage of this portable NumVec structure is that it can
// either own (owndata == true) or view (owndata == false) a piece of
// data.

template <class F>
class NumMat
{
public:
  Int m_, n_;
  bool owndata_;
  F* data_;
public:
  NumMat(Int m=0, Int n=0): m_(m), n_(n), owndata_(true) {
    if(m_>0 && n_>0) { data_ = new F[m_*n_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
  }
  NumMat(Int m, Int n, bool owndata, F* data): m_(m), n_(n), owndata_(owndata) {
    if(owndata_) {
      if(m_>0 && n_>0) { data_ = new F[m_*n_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
      if(m_>0 && n_>0) { for(Int i=0; i<m_*n_; i++) data_[i] = data[i]; }
    } else {
      data_ = data;
    }
  }
  NumMat(const NumMat& C): m_(C.m_), n_(C.n_), owndata_(C.owndata_) {
    if(owndata_) {
      if(m_>0 && n_>0) { data_ = new F[m_*n_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
      if(m_>0 && n_>0) { for(Int i=0; i<m_*n_; i++) data_[i] = C.data_[i]; }
    } else {
      data_ = C.data_;
    }
  }
  ~NumMat() {
    if(owndata_) {
      if(m_>0 && n_>0) { delete[] data_; data_ = NULL; }
    }
  }
  NumMat& operator=(const NumMat& C) {
    if(owndata_) {
      if(m_>0 && n_>0) { delete[] data_; data_ = NULL; }
    }
    m_ = C.m_; n_=C.n_; owndata_=C.owndata_;
    if(owndata_) {
      if(m_>0 && n_>0) { data_ = new F[m_*n_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
      if(m_>0 && n_>0) { for(Int i=0; i<m_*n_; i++) data_[i] = C.data_[i]; }
    } else {
      data_ = C.data_;
    }
    return *this;
  }
  void Resize(Int m, Int n)  {
		if( owndata_ == false ){
			throw std::logic_error("Matrix being resized must own data.");
		}
    if(m_!=m || n_!=n) {
      if(m_>0 && n_>0) { delete[] data_; data_ = NULL; }
      m_ = m; n_ = n;
      if(m_>0 && n_>0) { data_ = new F[m_*n_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
    }
		else {
			std::cerr << "NumMat<F>::Resize is doing nothing" << std::endl;
		}
  }
  const F& operator()(Int i, Int j) const  { 
		if( i < 0 || i >= m_ ||
				j < 0 || j >= n_ ) {
			throw std::logic_error( "Index is out of bound." );
		}
    return data_[i+j*m_];
  }
  F& operator()(Int i, Int j)  { 
		if( i < 0 || i >= m_ ||
				j < 0 || j >= n_ ) {
			throw std::logic_error( "Index is out of bound." );
		}
    return data_[i+j*m_];
  }
  
  F* Data() const { return data_; }
  F* VecData(Int j)  const 
	{ 
		if( j < 0 || j >= n_ ) {
			throw std::logic_error( "Index is out of bound." );
		}
		return &(data_[j*m_]); 
	}
  Int m() const { return m_; }
  Int n() const { return n_; }
};

// Commonly used
typedef NumMat<bool>     BolNumMat;
typedef NumMat<Int>      IntNumMat;
typedef NumMat<Real>     DblNumMat;
typedef NumMat<Complex>  CpxNumMat;

 
} // namespace PEXSI

#endif // _NUMMAT_DECL_HPP_
