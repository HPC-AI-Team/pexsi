#ifndef _NUMTNS_DECL_HPP_
#define _NUMTNS_DECL_HPP_

#include "environment_impl.hpp"
#include "nummat_impl.hpp"

// TODO Move the things from decl to impl

namespace  PEXSI{

// Templated form of numerical tensor
//
// The main advantage of this portable NumVec structure is that it can
// either own (owndata == true) or view (owndata == false) a piece of
// data.

template <class F>
class NumTns
{
public:
  Int m_, n_, p_;
  bool owndata_;
  F* data_;
public:
  NumTns(Int m=0, Int n=0, Int p=0): m_(m), n_(n), p_(p), owndata_(true) {
    if(m_>0 && n_>0 && p_>0) { data_ = new F[m_*n_*p_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
  }
  NumTns(Int m, Int n, Int p, bool owndata, F* data): m_(m), n_(n), p_(p), owndata_(owndata) {
    if(owndata_) {
      if(m_>0 && n_>0 && p_>0) { data_ = new F[m_*n_*p_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
      if(m_>0 && n_>0 && p_>0) { for(Int i=0; i<m_*n_*p_; i++) data_[i] = data[i]; }
    } else {
      data_ = data;
    }
  }
  NumTns(const NumTns& C): m_(C.m_), n_(C.n_), p_(C.p_), owndata_(C.owndata_) {
    if(owndata_) {
      if(m_>0 && n_>0 && p_>0) { data_ = new F[m_*n_*p_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
      if(m_>0 && n_>0 && p_>0) { for(Int i=0; i<m_*n_*p_; i++) data_[i] = C.data_[i]; }
    } else {
      data_ = C.data_;
    }
  }
  ~NumTns() { 
    if(owndata_) { 
      if(m_>0 && n_>0 && p_>0) { delete[] data_; data_ = NULL; } 
    }
  }
  NumTns& operator=(const NumTns& C) {
    if(owndata_) { 
      if(m_>0 && n_>0 && p_>0) { delete[] data_; data_ = NULL; } 
    }
    m_ = C.m_; n_=C.n_; p_=C.p_; owndata_=C.owndata_;
    if(owndata_) {
      if(m_>0 && n_>0 && p_>0) { data_ = new F[m_*n_*p_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
      if(m_>0 && n_>0 && p_>0) { for(Int i=0; i<m_*n_*p_; i++) data_[i] = C.data_[i]; }
    } else {
      data_ = C.data_;
    }
    return *this;
  }
  void Resize(Int m, Int n, Int p)  {
		if( owndata_ == false ){
			throw std::logic_error("Tensor being resized must own data.");
		}
    if(m_!=m || n_!=n || p_!=p) {
      if(m_>0 && n_>0 && p_>0) { delete[] data_; data_ = NULL; } 
      m_ = m; n_ = n; p_=p;
      if(m_>0 && n_>0 && p_>0) { data_ = new F[m_*n_*p_]; if( data_ == NULL ) throw std::runtime_error("Cannot allocate memory."); } else data_=NULL;
    }
  }
  const F& operator()(Int i, Int j, Int k) const  {
		if( i < 0 || i >= m_ ||
				j < 0 || j >= n_ ||
				k < 0 || k >= p_ ) {
			throw std::logic_error( "Index is out of bound." );
		}
    return data_[i+j*m_+k*m_*n_];
  }
  F& operator()(Int i, Int j, Int k)  {
		if( i < 0 || i >= m_ ||
				j < 0 || j >= n_ ||
				k < 0 || k >= p_ ) {
			throw std::logic_error( "Index is out of bound." );
		}
    return data_[i+j*m_+k*m_*n_];
  }
  
//  Int isempty() const {return (m_==0) && (n_==0) && (p_==0);}
  
  F* Data() const { return data_; }

	F* MatData (Int j) const {
		if( j < 0 || j >= p_ ) {
			throw std::logic_error( "Index is out of bound." );
		}
		return &(data_[j*m_*n_]);
	};

	F* VecData (Int j, Int k) const {
		if( j < 0 || j >= n_ ||
				k < 0 || k >= p_ ) {
			throw std::logic_error( "Index is out of bound." );
		}

		return &(data_[k*m_*n_+j*m_]);
	};

  Int m() const { return m_; }
  Int n() const { return n_; }
  Int p() const { return p_; }

};


typedef NumTns<bool>       BolNumTns;
typedef NumTns<Int>        IntNumTns;
typedef NumTns<Real>       DblNumTns;
typedef NumTns<Complex>    CpxNumTns;

} // namespace PEXSI

#endif // _NUMTNS_DECL_HPP_
