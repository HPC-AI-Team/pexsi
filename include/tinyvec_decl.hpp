/// @file tinyvec_decl.hpp
/// @brief Tiny vectors of dimension 3.
/// @author Lexing Ying and Lin Lin
/// @date 2010-09-20
#ifndef  _TINYVEC_DECL_HPP_
#define  _TINYVEC_DECL_HPP_

#include "environment_impl.hpp"

namespace PEXSI{

/// @class Vec3T
/// 
/// @brief Tiny vectors of dimension 3.
///
/// The main use of tiny vectors is to represent a triplet of three
/// indices (Index3) and three real numbers (Point3). 
template <class F>
	class Vec3T {
	private:
		F v_[3];
	public:
		enum{ X=0, Y=1, Z=2 };
		//------------CONSTRUCTOR AND DESTRUCTOR 
		Vec3T()              { v_[0]=F(0);    v_[1]=F(0);    v_[2]=F(0); }
		Vec3T(const F* f)    { v_[0]=f[0];    v_[1]=f[1];    v_[2]=f[2]; }
		Vec3T(const F a, const F b, const F c)   { v_[0]=a;       v_[1]=b;       v_[2]=c; }
		Vec3T(const Vec3T& c){ v_[0]=c.v_[0]; v_[1]=c.v_[1]; v_[2]=c.v_[2]; }
		~Vec3T() {}
		//------------POINTER and ACCESS
		operator F*()             { return &v_[0]; }
		operator const F*() const { return &v_[0]; }
		F* Data()                 { return &v_[0]; }  //access array
		F& operator()(Int i);
		const F& operator()(Int i) const;
		F& operator[](Int i);
		const F& operator[](Int i) const;
		//------------ASSIGN
		Vec3T& operator= ( const Vec3T& c ) { v_[0] =c.v_[0]; v_[1] =c.v_[1]; v_[2] =c.v_[2]; return *this; }
		Vec3T& operator+=( const Vec3T& c ) { v_[0]+=c.v_[0]; v_[1]+=c.v_[1]; v_[2]+=c.v_[2]; return *this; }
		Vec3T& operator-=( const Vec3T& c ) { v_[0]-=c.v_[0]; v_[1]-=c.v_[1]; v_[2]-=c.v_[2]; return *this; }
		Vec3T& operator*=( const F& s )     { v_[0]*=s;       v_[1]*=s;       v_[2]*=s;       return *this; }
		Vec3T& operator/=( const F& s )     { v_[0]/=s;       v_[1]/=s;       v_[2]/=s;       return *this; }
		//-----------LENGTH
		F l1( void )     const  { F sum=F(0); for(Int i=0; i<3; i++) sum=sum+std::abs(v_[i]); return sum; }
		F linfty( void ) const  { F cur=F(0); for(Int i=0; i<3; i++) cur=std::max(cur,std::abs(v_[i])); return cur; }
		F l2( void )     const  { F sum=F(0); for(Int i=0; i<3; i++) sum=sum+v_[i]*v_[i]; return sqrt(sum); }
	};

// *********************************************************************
// Most commonly used Vec3T thypes
// *********************************************************************
typedef Vec3T<Real>   Point3;
typedef Vec3T<Int>    Index3;

// *********************************************************************
// Compare
// *********************************************************************
template <class F> inline bool operator==(const Vec3T<F>& a, const Vec3T<F>& b) {
	return (a[0]==b[0] && a[1]==b[1] && a[2]==b[2]);
}
template <class F> inline bool operator!=(const Vec3T<F>& a, const Vec3T<F>& b) {
	return !(a==b);
}
template <class F> inline bool operator> (const Vec3T<F>& a, const Vec3T<F>& b) {
	for(Int i=0; i<3; i++) {
		if(     a[i]>b[i])	  return true;
		else if(a[i]<b[i])	  return false;
	}
	return false;
}
template <class F> inline bool operator< (const Vec3T<F>& a, const Vec3T<F>& b) {
	for(Int i=0; i<3; i++) {
		if(     a[i]<b[i])	  return true;
		else if(a[i]>b[i])	  return false;
	}
	return false;
}
template <class F> inline bool operator>=(const Vec3T<F>& a, const Vec3T<F>& b) {
	for(Int i=0; i<3; i++) {
		if(     a[i]>b[i])	  return true;
		else if(a[i]<b[i])	  return false;
	}
	return true;
}
template <class F> inline bool operator<=(const Vec3T<F>& a, const Vec3T<F>& b) {
	for(Int i=0; i<3; i++) {
		if(     a[i]<b[i])	  return true;
		else if(a[i]>b[i])	  return false;
	}
	return true;
}

// *********************************************************************
// Numerical operations
// *********************************************************************
template <class F> inline Vec3T<F> operator- (const Vec3T<F>& a) {
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = -a[i]; return r;
}
template <class F> inline Vec3T<F> operator+ (const Vec3T<F>& a, const Vec3T<F>& b) {
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = a[i]+b[i]; return r; 
}
template <class F> inline Vec3T<F> operator- (const Vec3T<F>& a, const Vec3T<F>& b) {
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = a[i]-b[i]; return r;
}
template <class F> inline Vec3T<F> operator* (F scl, const Vec3T<F>& a) {
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = scl*a[i];  return r;
}
template <class F> inline Vec3T<F> operator* (const Vec3T<F>& a, F scl) {
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = scl*a[i];  return r;
}
template <class F> inline Vec3T<F> operator/ (const Vec3T<F>& a, F scl) {
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = a[i]/scl;  return r;
}
template <class F> inline F operator* (const Vec3T<F>& a, const Vec3T<F>& b) {
	F sum=F(0); for(Int i=0; i<3; i++) sum=sum+a(i)*b(i); return sum;
}
template <class F> inline F dot       (const Vec3T<F>& a, const Vec3T<F>& b) {
	return a*b;
}
template <class F> inline Vec3T<F> operator^ (const Vec3T<F>& a, const Vec3T<F>& b) {
	return Vec3T<F>(a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0));
}
template <class F> inline Vec3T<F> cross     (const Vec3T<F>& a, const Vec3T<F>& b) { 
	return a^b; 
}

// *********************************************************************
// Element wise numerical operations
// *********************************************************************
template <class F> inline Vec3T<F> ewmin(const Vec3T<F>& a, const Vec3T<F>& b) {
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = std::min(a[i], b[i]); return r;
}
template <class F> inline Vec3T<F> ewmax(const Vec3T<F>& a, const Vec3T<F>& b) {
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = std::max(a[i], b[i]); return r;
}
template <class F> inline Vec3T<F> ewabs(const Vec3T<F>& a) {
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = std::abs(a[i]); return r;
}
template <class F> inline Vec3T<F> ewmul(const Vec3T<F>&a, const Vec3T<F>& b) {
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = a[i]*b[i]; return r;
}
template <class F> inline Vec3T<F> ewdiv(const Vec3T<F>&a, const Vec3T<F>& b) { 
	Vec3T<F> r;  for(Int i=0; i<3; i++) r[i] = a[i]/b[i]; return r;
}
template <class F> inline Vec3T<F> ewrnd(const Vec3T<F>&a) { //round
	Vec3T<F> r;  for(Int i=0; i<3; i++)	r[i] = round(a[i]);  return r;
}

// *********************************************************************
// Accumulative boolean operations
// *********************************************************************
template <class F> inline bool allequ(const Vec3T<F>& a, const Vec3T<F>& b) {
	bool res = true;  for(Int i=0; i<3; i++)   res = res && (a(i)==b(i));  return res;
}
template <class F> inline bool allneq(const Vec3T<F>& a, const Vec3T<F>& b) {
	return !(a==b);
}
template <class F> inline bool allgtt(const Vec3T<F>& a, const Vec3T<F>& b) {
	bool res = true;  for(Int i=0; i<3; i++)   res = res && (a(i)> b(i));  return res; 
}
template <class F> inline bool alllst(const Vec3T<F>& a, const Vec3T<F>& b) {
	bool res = true;  for(Int i=0; i<3; i++)   res = res && (a(i)< b(i));  return res; 
}
template <class F> inline bool allgoe(const Vec3T<F>& a, const Vec3T<F>& b) {
	bool res = true;  for(Int i=0; i<3; i++)	res = res && (a(i)>=b(i));  return res; 
}
template <class F> inline bool allloe(const Vec3T<F>& a, const Vec3T<F>& b) {
	bool res = true;  for(Int i=0; i<3; i++)   res = res && (a(i)<=b(i));  return res; 
}


// *********************************************************************
// Input and output
// *********************************************************************
template <class F> std::istream& operator>>(std::istream& is, Vec3T<F>& a) {
	for(Int i=0; i<3; i++) is>>a[i]; return is;
}
template <class F> std::ostream& operator<<(std::ostream& os, const Vec3T<F>& a) { 
	for(Int i=0; i<3; i++) os<<a[i]<<" "; return os;
}



} // namespace PEXSI

#endif // _TINYVEC_DECL_HPP_
