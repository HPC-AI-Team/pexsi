#ifndef SELINV_INTERF_H
#define SELINV_INTERF_H
#include "getpolef.h"

extern "C"
{
  int readmatrixheader_(char *, int *, int *);
 
  int readcmatrix_(char *, int *, int *, doublecomplex *);
  
  void sfinit_(int* outunt, int *neqns, int *nnza, int *xadj, int *adjncy,
	       int* perm, int* invp, int* colcnt, int* nnzl, int *nsub,
	       int* nsuper, int* snode, int* xsuper, int* iwsiz, int* iwork,
	       int* iflag );

  void bfinit_(int *neqns, int *nsuper, int *xsuper, int *snode, 
	       int *xlindx, int *lindx, int *cachsz, int *tmpsiz, int *split);	

  void inpnv_(int *outunt, int *neqns, int *xadjf, int *adjf,
	      doublecomplex *anzf, int *perm, int *invp, int *nsuper, int *xsuper,
	      int *xlindx, int *lindx, int *xlnz, doublecomplex *lnz,
	      int *iwsiz, int *offset, int *iflag);

  void blkfct_(int *outunt, int *nequns, int *nsuper, int *nunrol, int *xsuper,
	       int* snode, int* split, int* xlindx, int *lindx, int *xlnz,
	       doublecomplex *lnz, doublecomplex *diag, int *iwsiz, int *iwork,
	       int* tmpsiz, doublecomplex *tmpvec, int* iflag);

  void blkslv_(int *nsuper, int* xsuper, int *xlindx, int *lindx, int *xlnz,
	       doublecomplex *lnz, doublecomplex *rhs);

  void exdiag_(int *nsuper, int* xsuper, int *xlindx, int *lindx, int *xlnz,
	       doublecomplex *lnz, int *snodes, doublecomplex *diag, doublecomplex *y, int* perm,
	       int* neqns);

  void exdiagblk_(int *nsuper, int *xsuper, int *xlindx, int *lindx, int *xlnz,
		  doublecomplex *lnz, int *snodes, doublecomplex *diag, int *perm,
		  int *neqns, int *dumpl);

  void selinvblk_(int *nsuper, int *xsuper, int *xlindx, int *lindx, int *xlnz,
		  doublecomplex *lnz, int *snodes, int *perm, int *neqns,
		  int *colptr, int *rowind, int *dumpl);

  // in CSUPLDLT
  void ilo2ho_(int *n, int *Hnnz, int *colptr, int *rowind, int *newptr, int *newind, int *ip);

  void flo2ho_(int *n,          int *colptr, int *rowind,
	       doublecomplex *nzvals,  int *xadj,   int *adj,
	       doublecomplex *anz,     int *iwork);

  extern void ordmmd_(int *logfil, int *n,      int *xlindx, int *lindx, 
		      int *invp,   int *perm,   int *iwmax,  int *iwork,
		      int *nnzl,   int *nsub,   int *colcnt, int *nsuper,
		      int *xsuper, int *snodes, int *sfiflg, int *iflag); 

  extern void symfct_(int *logfil, int *n,      int *nnza,   int *xadj,
		      int *adj,    int *perm,   int *invp,   int *colcnt,
		      int *nsuper, int *xsuper, int *snodes, int *nsub,
		      int *xlindx, int *lindx,  int *xlnz,   int *iwsiz,
		      int *iwork,  int *iflag);

#ifdef METIS
  extern void METIS_EdgeND(int *n,       int *xadj, int *adj, int *numflag,
			   int *options, int *perm, int *invp);

  extern void METIS_NodeND(int *n,       int *xadj, int *adj, int *numflag,
			   int *options, int *perm, int *invp);
#endif

}; 


namespace SELINV_Interface
{
  // -----------------------
  //   Function prototypes
  // -----------------------
  void ldlt_preprocess__(int* token,  int *n,    int *colptr,
			 int* rowind, int *Lnnz, int *order, 
			 int *perm);

  void ldlt_fact__(int *token, int *colptr, int *rowind, doublecomplex *nzvals);

  void ldlt_solve__(int *token, doublecomplex *x, doublecomplex *rhs);

  void ldlt_getdiag__(int *token, doublecomplex *diag);

  void ldlt_extract__(int *token, doublecomplex *diag);

  void ldlt_blkextract__(int *token, doublecomplex *diag, int *dumpL);

  void ldlt_blkselinv__(int *token, int* colptr, int* rowind, 
			doublecomplex *inva, int *dumpL);

  void ldlt_free__(int *token);

}

#endif
