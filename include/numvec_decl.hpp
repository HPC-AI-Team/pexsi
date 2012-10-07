#ifndef _NUMVEC_DECL_HPP_
#define _NUMVEC_DECL_HPP_

#include "environment_impl.hpp"

namespace  PEXSI{

// Templated form of numerical vectors
//
// The main advantage of this portable NumVec structure is that it can
// either own (owndata == true) or view (owndata == false) a piece of
// data.

template <class F> class NumVec
{
public:
	Int  m_;                                // The number of elements 
	bool owndata_;                          // Whether it owns the data
	F* data_;                               // The pointer for the actual data
public:
	NumVec(Int m = 0);
	NumVec(Int m, bool owndata, F* data);
	NumVec(const NumVec& C);
	~NumVec();

	NumVec& operator=(const NumVec& C);

	void Resize ( Int m );

	const F& operator()(Int i) const;  
	F& operator()(Int i);  
	const F& operator[](Int i) const;
	F& operator[](Int i);

	bool IsOwnData() const { return owndata_; }
	F*   Data() const { return data_; }
	Int  m () const { return m_; }
};

// Commonly used
typedef NumVec<bool>       BolNumVec;
typedef NumVec<Int>        IntNumVec;
typedef NumVec<Real>       DblNumVec;
typedef NumVec<Complex>    CpxNumVec;

} // namespace PEXSI

#endif // _NUMVEC_DECL_HPP_
