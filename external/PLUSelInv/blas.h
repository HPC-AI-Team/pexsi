#ifndef _BLAS_H_
#define _BLAS_H_

//#include "bltypes.h"
#include <complex>
typedef std::complex<float> cpx8;
typedef std::complex<double> cpx16;
typedef int lint;

// LLIN: Commented because of conflict with SuperLU

 extern "C"
{
//  float sasum_(lint *n,float *x,lint *incx);
//  void  saxpy_(lint *n,float *alpha,float *x,lint *incx,float *y,lint *incy);
//  void  saxpyi_(lint *nz,float *a,float *x,lint *indx,float *y);
//  float scasum_(lint *n,cpx8 *x,lint *incx); 
//  float scnrm2_(lint *n,cpx8 *x,lint *incx); 
//  void  scopy_(lint *n,float *x,lint *incx,float *y,lint *incy);
//  float sdot_(lint *n,float *x,lint *incx,float *y,lint *incy);
//  float sdoti_(lint *nz,float *x,lint *indx,float *y);
//  void  sgthr_(lint *nz,float *y,float *x,lint *indx);
//  void  sgthrz_(lint *nz,float *y,float *x,lint *indx);
//  float snrm2_(lint *n,float *x,lint *incx);
//  void  srot_(lint *n,float *x,lint *incx,float *y,lint *incy,float *c,float *s);
//  void  srotg_(float *a,float *b,float *c,float *s);
//  void  sroti_(lint *nz,float *x,lint *indx,float *y,float *c,float *s);
//  void  srotm_(lint *n,float *x,lint *incx,float *y,lint *incy,float *param);
//  void  srotmg_(float *d1,float *d2,float *x1,float *y1,float *param);
//  void  sscal_(lint *n,float *a,float *x,lint *incx);
//  void  ssctr_(lint *nz,float *x,lint *indx,float *y);
//  void  sswap_(lint *n,float *x,lint *incx,float *y,lint *incy);
//  lint   isamax_(lint *n,float *x,lint *incx);
//  lint   isamin_(lint *n,float *x,lint *incx);
//  
//  void caxpy_(lint *n,cpx8 *alpha,cpx8 *x,lint *incx,cpx8 *y,lint *incy); 
//  void caxpyi_(lint *nz,cpx8 *a,cpx8 *x,lint *indx,cpx8 *y); 
//  void ccopy_(lint *n,cpx8 *x,lint *incx,cpx8 *y,lint *incy); 
//  void cdotc_(cpx8 *pres,lint *n,cpx8 *x,lint *incx,cpx8 *y,lint *incy); 
//  void cdotci_(cpx8 *pres,lint *nz,cpx8 *x,lint *indx,cpx8 *y); 
//  void cdotu_(cpx8 *pres,lint *n,cpx8 *x,lint *incx,cpx8 *y,lint *incy); 
//  void cdotui_(cpx8 *pres,lint *nz,cpx8 *x,lint *indx,cpx8 *y); 
//  void cgthr_(lint *nz,cpx8 *y,cpx8 *x,lint *indx); 
//  void cgthrz_(lint *nz,cpx8 *y,cpx8 *x,lint *indx); 
//  void crotg_(cpx8 *a,cpx8 *b,float *c,cpx8 *s); 
//  void cscal_(lint *n,cpx8 *a,cpx8 *x,lint *incx); 
//  void csctr_(lint *nz,cpx8 *x,lint *indx,cpx8 *y); 
//  void csrot_(lint *n,cpx8 *x,lint *incx,cpx8 *y,lint *incy,float *c,float *s); 
//  void csscal_(lint *n,float *a,cpx8 *x,lint *incx); 
//  void cswap_(lint *n,cpx8 *x,lint *incx,cpx8 *y,lint *incy); 
//  lint  icamax_(lint *n,cpx8 *x,lint *incx); 
//  lint  icamin_(lint *n,cpx8 *x,lint *incx); 
//  
//  double dasum_(lint *n,double *x,lint *incx);
//  void   daxpy_(lint *n,double *alpha,double *x,lint *incx,double *y,lint *incy);
//  void   daxpyi_(lint *nz,double *a,double *x,lint *indx,double *y);
//  void   dcopy_(lint *n,double *x,lint *incx,double *y,lint *incy);
//  double ddot_(lint *n,double *x,lint *incx,double *y,lint *incy);
//  double ddoti_(lint *nz,double *x,lint *indx,double *y);
//  void   dgthr_(lint *nz,double *y,double *x,lint *indx);
//  void   dgthrz_(lint *nz,double *y,double *x,lint *indx);
//  double dnrm2_(lint *n,double *x,lint *incx);
//  void   drot_(lint *n,double *x,lint *incx,double *y,lint *incy,double *c,double *s);
//  void   drotg_(double *a,double *b,double *c,double *s);
//  void   droti_(lint *nz,double *x,lint *indx,double *y,double *c,double *s);
//  void   drotm_(lint *n,double *x,lint *incx,double *y,lint *incy,double *param);
//  void   drotmg_(double *d1,double *d2,double *x1,double *y1,double *param);
//  void   dscal_(lint *n,double *a,double *x,lint *incx);
//  void   dsctr_(lint *nz,double *x,lint *indx,double *y);
//  void   dswap_(lint *n,double *x,lint *incx,double *y,lint *incy);
//  double dzasum_(lint *n,cpx16 *x,lint *incx); 
//  double dznrm2_(lint *n,cpx16 *x,lint *incx); 
//  lint    idamax_(lint *n,double *x,lint *incx);
//  lint    idamin_(lint *n,double *x,lint *incx);
//  
//  void zaxpy_(lint *n,cpx16 *alpha,cpx16 *x,lint *incx,cpx16 *y,lint *incy); 
//  void zaxpyi_(lint *nz,cpx16 *a,cpx16 *x,lint *indx,cpx16 *y); 
//  void zcopy_(lint *n,cpx16 *x,lint *incx,cpx16 *y,lint *incy); 
//  void zdotc_(cpx16 *pres,lint *n,cpx16 *x,lint *incx,cpx16 *y,lint *incy); 
//  void zdotci_(cpx16 *pres,lint *nz,cpx16 *x,lint *indx,cpx16 *y); 
//  void zdotu_(cpx16 *pres,lint *n,cpx16 *x,lint *incx,cpx16 *y,lint *incy); 
//  void zdotui_(cpx16 *pres,lint *nz,cpx16 *x,lint *indx,cpx16 *y); 
//  void zdrot_(lint *n,cpx16 *x,lint *incx,cpx16 *y,lint *incy,double *c,double *s); 
//  void zdscal_(lint *n,double *a,cpx16 *x,lint *incx); 
//  void zgthr_(lint *nz,cpx16 *y,cpx16 *x,lint *indx); 
//  void zgthrz_(lint *nz,cpx16 *y,cpx16 *x,lint *indx); 
//  void zrotg_(cpx16 *a,cpx16 *b,double *c,cpx16 *s); 
//  void zscal_(lint *n,cpx16 *a,cpx16 *x,lint *incx); 
//  void zsctr_(lint *nz,cpx16 *x,lint *indx,cpx16 *y); 
//  void zswap_(lint *n,cpx16 *x,lint *incx,cpx16 *y,lint *incy); 
//  lint  izamax_(lint *n,cpx16 *x,lint *incx); 
//  lint  izamin_(lint *n,cpx16 *x,lint *incx); 
//  
//  /* blas level2 */
//  
//  void sgbmv_(char *trans,lint *m,lint *n,lint *kl,lint *ku,float *alpha,float *a,lint *lda,float *x,lint *incx,float *beta,float *y,lint *incy);
//  void sgemv_(char *trans,lint *m,lint *n,float *alpha,float *a,lint *lda,float *x,lint *incx,float *beta,float *y,lint *incy);
//  void sger_(lint *m,lint *n,float *alpha,float *x,lint *incx,float *y,lint *incy,float *a,lint *lda);
//  void ssbmv_(char *uplo,lint *n,lint *k,float *alpha,float *a,lint *lda,float *x,lint *incx,float *beta,float *y,lint *incy);
//  void sspmv_(char *uplo,lint *n,float *alpha,float *ap,float *x,lint *incx,float *beta,float *y,lint *incy);
//  void sspr_(char *uplo,lint *n,float *alpha,float *x,lint *incx,float *ap);
//  void sspr2_(char *uplo,lint *n,float *alpha,float *x,lint *incx,float *y,lint *incy,float *ap);
//  void ssymv_(char *uplo,lint *n,float *alpha,float *a,lint *lda,float *x,lint *incx,float *beta,float *y,lint *incy);
//  void ssyr_(char *uplo,lint *n,float *alpha,float *x,lint *incx,float *a,lint *lda);
//  void ssyr2_(char *uplo,lint *n,float *alpha,float *x,lint *incx,float *y,lint *incy,float *a,lint *lda);
//  void stbmv_(char *uplo,char *trans,char *diag,lint *n,lint *k,float *a,lint *lda,float *x,lint *incx);
//  void stbsv_(char *uplo,char *trans,char *diag,lint *n,lint *k,float *a,lint *lda,float *x,lint *incx);
//  void stpmv_(char *uplo,char *trans,char *diag,lint *n,float *ap,float *x,lint *incx);
//  void stpsv_(char *uplo,char *trans,char *diag,lint *n,float *ap,float *x,lint *incx);
//  void strmv_(char *uplo,char *transa,char *diag,lint *n,float *a,lint *lda,float *b,lint *incx);
//  void strsv_(char *uplo,char *trans,char *diag,lint *n,float *a,lint *lda,float *x,lint *incx);
//  
//  void cgbmv_(char *trans,lint *m,lint *n,lint *kl,lint *ku,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *x,lint *incx,cpx8 *beta,cpx8 *y,lint *incy); 
//  void cgemv_(char *trans,lint *m,lint *n,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *x,lint *incx,cpx8 *beta,cpx8 *y,lint *incy); 
//  void cgerc_(lint *m,lint *n,cpx8 *alpha,cpx8 *x,lint *incx,cpx8 *y,lint *incy,cpx8 *a,lint *lda); 
//  void cgeru_(lint *m,lint *n,cpx8 *alpha,cpx8 *x,lint *incx,cpx8 *y,lint *incy,cpx8 *a,lint *lda); 
//  void chbmv_(char *uplo,lint *n,lint *k,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *x,lint *incx,cpx8 *beta,cpx8 *y,lint *incy); 
//  void chemv_(char *uplo,lint *n,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *x,lint *incx,cpx8 *beta,cpx8 *y,lint *incy); 
//  void cher_(char *uplo,lint *n,float *alpha,cpx8 *x,lint *incx,cpx8 *a,lint *lda); 
//  void cher2_(char *uplo,lint *n,cpx8 *alpha,cpx8 *x,lint *incx,cpx8 *y,lint *incy,cpx8 *a,lint *lda); 
//  void chpmv_(char *uplo,lint *n,cpx8 *alpha,cpx8 *ap,cpx8 *x,lint *incx,cpx8 *beta,cpx8 *y,lint *incy); 
//  void chpr_(char *uplo,lint *n,float *alpha,cpx8 *x,lint *incx,cpx8 *ap); 
//  void chpr2_(char *uplo,lint *n,cpx8 *alpha,cpx8 *x,lint *incx,cpx8 *y,lint *incy,cpx8 *ap); 
//  void ctbmv_(char *uplo,char *trans,char *diag,lint *n,lint *k,cpx8 *a,lint *lda,cpx8 *x,lint *incx); 
//  void ctbsv_(char *uplo,char *trans,char *diag,lint *n,lint *k,cpx8 *a,lint *lda,cpx8 *x,lint *incx); 
//  void ctpmv_(char *uplo,char *trans,char *diag,lint *n,cpx8 *ap,cpx8 *x,lint *incx); 
//  void ctpsv_(char *uplo,char *trans,char *diag,lint *n,cpx8 *ap,cpx8 *x,lint *incx); 
//  void ctrmv_(char *uplo,char *transa,char *diag,lint *n,cpx8 *a,lint *lda,cpx8 *b,lint *incx); 
//  void ctrsv_(char *uplo,char *trans,char *diag,lint *n,cpx8 *a,lint *lda,cpx8 *x,lint *incx); 
//  
//  void dgbmv_(char *trans,lint *m,lint *n,lint *kl,lint *ku,double *alpha,double *a,lint *lda,double *x,lint *incx,double *beta,double *y,lint *incy);
//  void dgemv_(char *trans,lint *m,lint *n,double *alpha,double *a,lint *lda,double *x,lint *incx,double *beta,double *y,lint *incy);
//  void dger_(lint *m,lint *n,double *alpha,double *x,lint *incx,double *y,lint *incy,double *a,lint *lda);
//  void dsbmv_(char *uplo,lint *n,lint *k,double *alpha,double *a,lint *lda,double *x,lint *incx,double *beta,double *y,lint *incy);
//  void dspmv_(char *uplo,lint *n,double *alpha,double *ap,double *x,lint *incx,double *beta,double *y,lint *incy);
//  void dspr_(char *uplo,lint *n,double *alpha,double *x,lint *incx,double *ap);
//  void dspr2_(char *uplo,lint *n,double *alpha,double *x,lint *incx,double *y,lint *incy,double *ap);
//  void dsymv_(char *uplo,lint *n,double *alpha,double *a,lint *lda,double *x,lint *incx,double *beta,double *y,lint *incy);
//  void dsyr_(char *uplo,lint *n,double *alpha,double *x,lint *incx,double *a,lint *lda);
//  void dsyr2_(char *uplo,lint *n,double *alpha,double *x,lint *incx,double *y,lint *incy,double *a,lint *lda);
//  void dtbmv_(char *uplo,char *trans,char *diag,lint *n,lint *k,double *a,lint *lda,double *x,lint *incx);
//  void dtbsv_(char *uplo,char *trans,char *diag,lint *n,lint *k,double *a,lint *lda,double *x,lint *incx);
//  void dtpmv_(char *uplo,char *trans,char *diag,lint *n,double *ap,double *x,lint *incx);
//  void dtpsv_(char *uplo,char *trans,char *diag,lint *n,double *ap,double *x,lint *incx);
//  void dtrmv_(char *uplo,char *transa,char *diag,lint *n,double *a,lint *lda,double *b,lint *incx);
//  void dtrsv_(char *uplo,char *trans,char *diag,lint *n,double *a,lint *lda,double *x,lint *incx);
//  
//  void zgbmv_(char *trans,lint *m,lint *n,lint *kl,lint *ku,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *x,lint *incx,cpx16 *beta,cpx16 *y,lint *incy); 
//  void zgemv_(char *trans,lint *m,lint *n,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *x,lint *incx,cpx16 *beta,cpx16 *y,lint *incy); 
//  void zgerc_(lint *m,lint *n,cpx16 *alpha,cpx16 *x,lint *incx,cpx16 *y,lint *incy,cpx16 *a,lint *lda); 
//  void zgeru_(lint *m,lint *n,cpx16 *alpha,cpx16 *x,lint *incx,cpx16 *y,lint *incy,cpx16 *a,lint *lda); 
//  void zhbmv_(char *uplo,lint *n,lint *k,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *x,lint *incx,cpx16 *beta,cpx16 *y,lint *incy); 
//  void zhemv_(char *uplo,lint *n,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *x,lint *incx,cpx16 *beta,cpx16 *y,lint *incy); 
//  void zher_(char *uplo,lint *n,double *alpha,cpx16 *x,lint *incx,cpx16 *a,lint *lda); 
//  void zher2_(char *uplo,lint *n,cpx16 *alpha,cpx16 *x,lint *incx,cpx16 *y,lint *incy,cpx16 *a,lint *lda); 
//  void zhpmv_(char *uplo,lint *n,cpx16 *alpha,cpx16 *ap,cpx16 *x,lint *incx,cpx16 *beta,cpx16 *y,lint *incy); 
//  void zhpr_(char *uplo,lint *n,double *alpha,cpx16 *x,lint *incx,cpx16 *ap); 
//  void zhpr2_(char *uplo,lint *n,cpx16 *alpha,cpx16 *x,lint *incx,cpx16 *y,lint *incy,cpx16 *ap); 
//  void ztbmv_(char *uplo,char *trans,char *diag,lint *n,lint *k,cpx16 *a,lint *lda,cpx16 *x,lint *incx); 
//  void ztbsv_(char *uplo,char *trans,char *diag,lint *n,lint *k,cpx16 *a,lint *lda,cpx16 *x,lint *incx); 
//  void ztpmv_(char *uplo,char *trans,char *diag,lint *n,cpx16 *ap,cpx16 *x,lint *incx); 
//  void ztpsv_(char *uplo,char *trans,char *diag,lint *n,cpx16 *ap,cpx16 *x,lint *incx); 
//  void ztrmv_(char *uplo,char *transa,char *diag,lint *n,cpx16 *a,lint *lda,cpx16 *b,lint *incx); 
//  void ztrsv_(char *uplo,char *trans,char *diag,lint *n,cpx16 *a,lint *lda,cpx16 *x,lint *incx); 
//  
//  /* blas level3 */
//  
//  void sgemm_(char *transa,char *transb,lint *m,lint *n,lint *k,float *alpha,float *a,lint *lda,float *b,lint *ldb,float *beta,float *c,lint *ldc);
//  void ssymm_(char *side,char *uplo,lint *m,lint *n,float *alpha,float *a,lint *lda,float *b,lint *ldb,float *beta,float *c,lint *ldc);
//  void ssyr2k_(char *uplo,char *trans,lint *n,lint *k,float *alpha,float *a,lint *lda,float *b,lint *ldb,float *beta,float *c,lint *ldc);
//  void ssyrk_(char *uplo,char *trans,lint *n,lint *k,float *alpha,float *a,lint *lda,float *beta,float *c,lint *ldc);
//  void strmm_(char *side,char *uplo,char *transa,char *diag,lint *m,lint *n,float *alpha,float *a,lint *lda,float *b,lint *ldb);
//  void strsm_(char *side,char *uplo,char *transa,char *diag,lint *m,lint *n,float *alpha,float *a,lint *lda,float *b,lint *ldb);
//  
//  void cgemm_(char *transa,char *transb,lint *m,lint *n,lint *k,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *b,lint *ldb,cpx8 *beta,cpx8 *c,lint *ldc); 
//  void chemm_(char *side,char *uplo,lint *m,lint *n,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *b,lint *ldb,cpx8 *beta,cpx8 *c,lint *ldc); 
//  void cher2k_(char *uplo,char *trans,lint *n,lint *k,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *b,lint *ldb,float *beta,cpx8 *c,lint *ldc); 
//  void cherk_(char *uplo,char *trans,lint *n,lint *k,float *alpha,cpx8 *a,lint *lda,float *beta,cpx8 *c,lint *ldc); 
//  void csymm_(char *side,char *uplo,lint *m,lint *n,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *b,lint *ldb,cpx8 *beta,cpx8 *c,lint *ldc); 
//  void csyr2k_(char *uplo,char *trans,lint *n,lint *k,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *b,lint *ldb,cpx8 *beta,cpx8 *c,lint *ldc); 
//  void csyrk_(char *uplo,char *trans,lint *n,lint *k,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *beta,cpx8 *c,lint *ldc); 
//  void ctrmm_(char *side,char *uplo,char *transa,char *diag,lint *m,lint *n,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *b,lint *ldb); 
//  void ctrsm_(char *side,char *uplo,char *transa,char *diag,lint *m,lint *n,cpx8 *alpha,cpx8 *a,lint *lda,cpx8 *b,lint *ldb); 
//  
//  void dgemm_(char *transa,char *transb,lint *m,lint *n,lint *k,double *alpha,double *a,lint *lda,double *b,lint *ldb,double *beta,double *c,lint *ldc);
//  void dsymm_(char *side,char *uplo,lint *m,lint *n,double *alpha,double *a,lint *lda,double *b,lint *ldb,double *beta,double *c,lint *ldc);
//  void dsyr2k_(char *uplo,char *trans,lint *n,lint *k,double *alpha,double *a,lint *lda,double *b,lint *ldb,double *beta,double *c,lint *ldc);
//  void dsyrk_(char *uplo,char *trans,lint *n,lint *k,double *alpha,double *a,lint *lda,double *beta,double *c,lint *ldc);
//  void dtrmm_(char *side,char *uplo,char *transa,char *diag,lint *m,lint *n,double *alpha,double *a,lint *lda,double *b,lint *ldb);
  void dtrsm_(char *side,char *uplo,char *transa,char *diag,lint *m,lint *n,double *alpha,double *a,lint *lda,double *b,lint *ldb);
//  
//  void zgemm_(char *transa,char *transb,lint *m,lint *n,lint *k,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *b,lint *ldb,cpx16 *beta,cpx16 *c,lint *ldc); 
//  void zgemm_(char *transa,char *transb,lint *m,lint *n,lint *k,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *b,lint *ldb,cpx16 *beta,cpx16 *c,lint *ldc); 
//  void zhemm_(char *side,char *uplo,lint *m,lint *n,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *b,lint *ldb,cpx16 *beta,cpx16 *c,lint *ldc); 
//  void zher2k_(char *uplo,char *trans,lint *n,lint *k,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *b,lint *ldb,double *beta,cpx16 *c,lint *ldc); 
//  void zherk_(char *uplo,char *trans,lint *n,lint *k,double *alpha,cpx16 *a,lint *lda,double *beta,cpx16 *c,lint *ldc); 
//  void zsymm_(char *side,char *uplo,lint *m,lint *n,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *b,lint *ldb,cpx16 *beta,cpx16 *c,lint *ldc); 
//  void zsyr2k_(char *uplo,char *trans,lint *n,lint *k,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *b,lint *ldb,cpx16 *beta,cpx16 *c,lint *ldc); 
//  void zsyrk_(char *uplo,char *trans,lint *n,lint *k,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *beta,cpx16 *c,lint *ldc); 
//  void ztrmm_(char *side,char *uplo,char *transa,char *diag,lint *m,lint *n,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *b,lint *ldb); 
//  void ztrsm_(char *side,char *uplo,char *transa,char *diag,lint *m,lint *n,cpx16 *alpha,cpx16 *a,lint *lda,cpx16 *b,lint *ldb); 
}

#endif

