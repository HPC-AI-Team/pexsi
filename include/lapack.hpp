#include  "environment_impl.hpp"

namespace PEXSI {
namespace lapack {

typedef  int                    Int; 
typedef  std::complex<float>    scomplex;
typedef  std::complex<double>   dcomplex;

//
// Machine constants
//

// Relative machine precision
template<typename R> R MachineEpsilon();
template<> float MachineEpsilon<float>();
template<> double MachineEpsilon<double>();

// Minimum number which can be inverted without overflow
template<typename R> R MachineSafeMin();
template<> float MachineSafeMin<float>();
template<> double MachineSafeMin<double>();

// Base of the machine, where the number is represented as 
//   (mantissa) x (base)^(exponent)
template<typename R> R MachineBase();
template<> float MachineBase<float>();
template<> double MachineBase<double>();

// Return the relative machine precision multiplied by the base
template<typename R> R MachinePrecision();
template<> float MachinePrecision<float>();
template<> double MachinePrecision<double>();

// Return the minimum exponent before (gradual) underflow occurs
template<typename R> R MachineUnderflowExponent();
template<> float MachineUnderflowExponent<float>();
template<> double MachineUnderflowExponent<double>();

// Return the underflow threshold: (base)^((underflow exponent)-1)
template<typename R> R MachineUnderflowThreshold();
template<> float MachineUnderflowThreshold<float>();
template<> double MachineUnderflowThreshold<double>();

// Return the largest exponent before overflow
template<typename R> R MachineOverflowExponent();
template<> float MachineOverflowExponent<float>();
template<> double MachineOverflowExponent<double>();

// Return the overflow threshold: (1-(rel. prec.)) * (base)^(overflow exponent)
template<typename R> R MachineOverflowThreshold();
template<> float MachineOverflowThreshold<float>();
template<> double MachineOverflowThreshold<double>();

//
// For safely computing norms without overflow/underflow
//

float SafeNorm( float alpha, float beta );
double SafeNorm( double alpha, double beta );
float SafeNorm( float alpha, float beta, float gamma );
double SafeNorm( double alpha, double beta, double gamma );


//
//
// Given phi and gamma, compute a Givens rotation such that
//
//  |       cs   sn | |   phi |  = | rho |, where cs^2 + |sn|^2 = 1
//  | -conj(sn)  cs | | gamma |    |  0  |
//
// This routine should use the stable approach suggested by Kahan and Demmel
//

void ComputeGivens
( float phi, float gamma,
  float* cs, float* sn, float* rho );

void ComputeGivens
( double phi, double gamma,
  double* cs, double* sn, double* rho );

void ComputeGivens
( scomplex phi, scomplex gamma,
  float* cs, scomplex* sn, scomplex* rho );

void ComputeGivens
( dcomplex phi, dcomplex gamma,
  double* cs, dcomplex* sn, dcomplex* rho );

//
// Cholesky factorization
//

void Cholesky( char uplo, Int n, const float* A, Int lda );
void Cholesky( char uplo, Int n, const double* A, Int lda );
void Cholesky( char uplo, Int n, const scomplex* A, Int lda );
void Cholesky( char uplo, Int n, const dcomplex* A, Int lda );

//
// LU factorization (with partial pivoting)
//

void LU( Int m, Int n, float* A, Int lda, Int* p );
void LU( Int m, Int n, double* A, Int lda, Int* p );
void LU( Int m, Int n, scomplex* A, Int lda, Int* p );
void LU( Int m, Int n, dcomplex* A, Int lda, Int* p );

//
// For reducing well-conditioned Hermitian generalized-definite EVP's
// to standard form.
//

void Hegst
( Int itype, char uplo, 
  Int n, float* A, Int lda, const float* B, Int ldb );
void Hegst
( Int itype, char uplo,
  Int n, double* A, Int lda, const double* B, Int ldb );
void Hegst
( Int itype, char uplo,
  Int n, scomplex* A, Int lda, const scomplex* B, Int ldb );
void Hegst
( Int itype, char uplo,
  Int n, dcomplex* A, Int lda, const dcomplex* B, Int ldb );

//
// For computing the inverse of a triangular matrix
//

void TriangularInverse
( char uplo, char diag, Int n, const float* A, Int lda );
void TriangularInverse
( char uplo, char diag, Int n, const double* A, Int lda );
void TriangularInverse
( char uplo, char diag, Int n, const scomplex* A, Int lda );
void TriangularInverse
( char uplo, char diag, Int n, const dcomplex* A, Int lda );

//
// Compute the SVD of a general matrix using a divide and conquer algorithm
//

void DivideAndConquerSVD
( Int m, Int n, float* A, Int lda, 
  float* s, float* U, Int ldu, float* VTrans, Int ldvt );
void DivideAndConquerSVD
( Int m, Int n, double* A, Int lda, 
  double* s, double* U, Int ldu, double* VTrans, Int ldvt );
void DivideAndConquerSVD
( Int m, Int n, scomplex* A, Int lda, 
  float* s, scomplex* U, Int ldu, scomplex* VAdj, Int ldva );
void DivideAndConquerSVD
( Int m, Int n, dcomplex* A, Int lda, 
  double* s, dcomplex* U, Int ldu, dcomplex* VAdj, Int ldva );

//
// Compute the SVD of a general matrix using the QR algorithm
//

void QRSVD
( Int m, Int n, float* A, Int lda, 
  float* s, float* U, Int ldu, float* VTrans, Int ldvt );
void QRSVD
( Int m, Int n, double* A, Int lda, 
  double* s, double* U, Int ldu, double* VTrans, Int ldvt );
void QRSVD
( Int m, Int n, scomplex* A, Int lda, 
  float* s, scomplex* U, Int ldu, scomplex* VAdj, Int ldva );
void QRSVD
( Int m, Int n, dcomplex* A, Int lda, 
  double* s, dcomplex* U, Int ldu, dcomplex* VAdj, Int ldva );

//
// Compute the singular values of a general matrix using the QR algorithm
//

void SingularValues( Int m, Int n, float* A, Int lda, float* s );
void SingularValues( Int m, Int n, double* A, Int lda, double* s );
void SingularValues( Int m, Int n, scomplex* A, Int lda, float* s );
void SingularValues( Int m, Int n, dcomplex* A, Int lda, double* s );

//
// Compute the SVD of a bidiagonal matrix using the QR algorithm
//

void BidiagQRAlg
( char uplo, Int n, Int numColsVTrans, Int numRowsU,
  float* d, float* e, float* VTrans, Int ldVTrans, float* U, Int ldU );
void BidiagQRAlg
( char uplo, Int n, Int numColsVTrans, Int numRowsU, 
  double* d, double* e, double* VTrans, Int ldVTrans, double* U, Int ldU );
void BidiagQRAlg
( char uplo, Int n, Int numColsVAdj, Int numRowsU,
  float* d, float* e, scomplex* VAdj, Int ldVAdj, scomplex* U, Int ldU );
void BidiagQRAlg
( char uplo, Int n, Int numColsVAdj, Int numRowsU, 
  double* d, double* e, dcomplex* VAdj, Int ldVAdj, dcomplex* U, Int ldU );

// *********************************************************************
// Compute the linear least square problem using SVD
// *********************************************************************
void SVDLeastSquare( Int m, Int n, Int nrhs, float * A, Int lda,
		float * B, Int ldb, float * S, float rcond,
		Int* rank );
void SVDLeastSquare( Int m, Int n, Int nrhs, double * A, Int lda,
		double * B, Int ldb, double * S, double rcond,
		Int* rank );
void SVDLeastSquare( Int m, Int n, Int nrhs, scomplex * A, Int lda,
		scomplex * B, Int ldb, float * S, float rcond,
		Int* rank );
void SVDLeastSquare( Int m, Int n, Int nrhs, dcomplex * A, Int lda,
		dcomplex * B, Int ldb, double * S, double rcond,
		Int* rank );

// *********************************************************************
// Copy
// *********************************************************************

void Copy( char uplo, Int m, Int n, const dcomplex* A, Int lda,
	dcomplex* B, Int ldb	);



} // namespace lapack
} // namespace PEXSI
