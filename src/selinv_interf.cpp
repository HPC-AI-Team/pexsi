/// @file selinv_interf.cpp
/// @brief Implementation of the interface with sequential SelInv.
/// @author Chao Yang and Lin Lin
/// @date 2012-10-16
#include <stdio.h>
#include <stdlib.h>
#include "selinv_interf.hpp"
// TODO Real arithmetic

// METIS is not used by default 
//#define METIS 1

namespace PEXSI{

#ifdef _USE_COMPLEX_
typedef struct { double r, i; } doublecomplex;
#endif


// *********************************************************************
// FORTRAN subroutines for factorization and selected inversion
// *********************************************************************
#ifdef _USE_COMPLEX_
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
#endif


// *********************************************************************
// Definition of interface subroutines with FORTRAN selinv
// *********************************************************************
#ifdef _USE_COMPLEX_
void ldlt_preprocess__(int* token,  int *n,    int *colptr,
		int* rowind, int *Lnnz, int *order, 
		int *perm);

void ldlt_fact__(int *token, int *colptr, int *rowind, doublecomplex *nzvals);

void ldlt_solve__(int *token, doublecomplex *x, doublecomplex *rhs);

void ldlt_getdiag__(int *token, doublecomplex *diag);

void ldlt_blkselinv__(int *token, int* colptr, int* rowind, 
		doublecomplex *inva, int *dumpL);

void ldlt_free__(int *token);


// *********************************************************************
// Global Data structure and variables shared by the interface
// subroutines
// *********************************************************************

typedef struct Anode {
	int    n;
	int    nsuper;
	int    nsub;
	int    nnzl;
	int    *xadj;
	int    *adj;
	doublecomplex *anz;
	doublecomplex *adiag;
	int    *xsuper;
	int    *snodes;
	int    *xlindx;
	int    *lindx;
	int    *xlnz;
	doublecomplex *lnz;
	doublecomplex *diag;
	doublecomplex *tmat;
	int    *perm;
	int    *invp;
	int    *colcnt;
	int    *iwork;
	int    tmpsiz;
	int    *split;
	doublecomplex *newrhs;
}  Anode_type;

// Global variables

// Some global parameters to be optimized
// Only use one matrix here.
const int    MAX_NUM_MAT = 1;
int          cachsz=700; 
int          nunroll=4;
int          fullrep = 0; 
Anode_type   mat[MAX_NUM_MAT]; 



// *********************************************************************
// Implementation of interface subroutines with FORTRAN selinv
// *********************************************************************
void ldlt_preprocess__(int *token, int *n,    int *colptr,
		int *rowind, int *Lnnz, int *order,
		int *perm)
{
	int  numflag, options[8];
	int  nnz,     nnza,  neqns,  nsub,    ibegin, iend,
			 i,        j,  nnzl,   Lnsub,   nsuper, iwsiz, 
			 clnz,    snsz,  tmpsze, tmpsiz,  maxlen, maxsup,
			 rnnz,    ierr,  jglb,   iflag,   sfiflg, logfil=6;
	int  *adj2, *xadj2;

	int  pjend, pibeg, notfound, irow, jcol;
	int  MAXORDER;

#ifdef TIMING
	double t0, t1;
#endif
	/*DEBUG begin*/
	FILE *fp;
	/*DEBUG end*/

	neqns  = *n;
	nnz    = colptr[*n]-1;

	// Check to see whether a full representation the symmetric
	//   matrix is used 
	notfound = 1;
	j = 1;
	while (notfound && j < neqns) {
		/* pointer to the last entry of column j */
		pjend = colptr[j]-1;
		if (pjend >= colptr[j-1]) {
			/* there is at least one off-diagonal entry,
				 find its row index */
			irow  = rowind[pjend-1];

			/* get column pointer to the first entry of irow-th column */
			pibeg = colptr[irow-1];

			/* check to see if upper triangular part is present */
			if ( rowind[pibeg-1] != irow ) {
				fullrep = 1;
			}
			notfound = 0;
		}
		j++;
	}

	/* nnza = number of nonzeros off diagonal */
	nnza   = 2*(nnz-neqns);
	if (fullrep) nnza = nnz-neqns;

	if (fullrep) fprintf(stderr, " Full Respresentation of the Matrix?\n");

	iwsiz  = 7*neqns+3;
	mat[*token].n = neqns;

	// FIXME magic number
	if (*token > 5 || *token < 0) { 
		fprintf(stderr," ldlt_preprocess: Invalid token number!\n");
		ierr = 1;
		exit(1);
	}
	else {

		/* --------------------------------------------------
			 Allocate matrix storage required for reordering
			 and symbolic factorization.
			 -------------------------------------------------- */

		mat[*token].n = *n;

		/* adj, xadj contain the structure of the
			 full representation of the original matrix */

		mat[*token].xadj = (int*) malloc((neqns+1)*sizeof(int));
		if (!mat[*token].xadj) {
			fprintf(stderr, 
					"ldlt_preprocess: memory allocation failed for xadj\n");
			exit(1);
		} 

		mat[*token].adj = (int*) malloc((nnza+neqns)*sizeof(int));
		if (!mat[*token].adj) {
			fprintf(stderr,
					"ldlt_preprocess: memory allocation failed for adj\n");
			exit(1);
		}

		/* Since reordering will destroy the original matrix
			 we need to make an extra copy of the indices and pointers
			 of the full representation of the original matrix.  */

		xadj2 = (int *)malloc((neqns+1)*sizeof(int));
		if (!xadj2) {
			fprintf(stderr,
					"ldlt_preprocess: memory allocation failed for adj2\n");
			exit(1);
		}

		adj2 = (int *)malloc(nnza*sizeof(int));
		if (!adj2) {
			fprintf(stderr,
					"ldlt_preprocess: memory allocation failed for xadj2\n");
			exit(1);
		}

		/* ----------------------------------------------------------
			 convert the indices and pointers of the lower triangular 
			 representation of a symmetric HB matrix to its full 
			 representation. 
			 ---------------------------------------------------------- */

		mat[*token].iwork = (int*) malloc(iwsiz*sizeof(int));
		if (!mat[*token].iwork) { 
			fprintf(stderr, "ldlt_preprocess: Fail to allocate iwork\n");
			exit(1);
		}

		if (fullrep) {
			/* make a copy of the non-zero structure, take out the 
				 diagonals */
			mat[*token].xadj[0] = 1;
			for (i=0;i<neqns;i++) {
				mat[*token].xadj[i+1] = mat[*token].xadj[i] 
					+ (colptr[i+1] - colptr[i] - 1);
			}
			if (mat[*token].xadj[neqns]-1 != nnza) {
				fprintf(stderr, " Something is wrong with the matrix!\n");
			}

			i=0;
			for (jcol = 1; jcol <= neqns; jcol++) {
				ibegin = colptr[jcol-1];
				iend   = colptr[jcol]-1;
				for (irow = ibegin; irow <=iend; irow++) {
					if (rowind[irow-1] != jcol) {
						mat[*token].adj[i] = rowind[irow-1];
						i++;
					}
				}
			}
		}
		else { 
			// convert 
			ilo2ho_(n,                 &nnza,            colptr,
					rowind,            mat[*token].xadj, mat[*token].adj,
					mat[*token].iwork);
		} 

		/* DEBUG */    
		/* printf(" preprocess: adj[nnz-1] = %d\n",
			 mat[*token].adj[nnz-1]); */

		/* -------------------------------------------------
			 Make a copy of the convertd matrix for reordering 
			 ------------------------------------------------- */

		for (i=0;i<=neqns;i++)   xadj2[i] = mat[*token].xadj[i];
		for (i=0;i<nnza;i++)     adj2[i]  = mat[*token].adj[i];

		/* ----------------------------------------
			 Allocate storage for supernode partition 
			 ----------------------------------------*/

		mat[*token].snodes = (int*)malloc(neqns*sizeof(int));
		if (!mat[*token].snodes) {
			fprintf(stderr, 
					"ldlt_preprocess: memory allocation failed for snodes\n");
			exit(1);
		}

		mat[*token].xsuper = (int*)malloc((neqns+1)*sizeof(int)); 
		if (!mat[*token].xsuper) {
			fprintf(stderr,
					"ldlt_preprocess: memory allocation failed for xsuper\n");
			exit(1);
		}

		mat[*token].colcnt = (int*)malloc(neqns*sizeof(int));
		if (!mat[*token].colcnt) {
			fprintf(stderr,
					"ldlt_preprocess: memory allocation failed for colcnt\n");
			exit(1);
		} 

		mat[*token].perm = (int*)malloc(neqns*sizeof(int));
		if (!mat[*token].perm) {
			fprintf(stderr, 
					"ldlt_preprocess: memory allocation failed for perm\n");
			exit(1);
		}

		mat[*token].invp = (int*)malloc(neqns*sizeof(int));
		if (!mat[*token].invp) {
			fprintf(stderr, 
					"ldlt_preprocess: memory allocation failed for invp\n");
			exit(1);
		}

#ifdef TIMING
		t0 = getime();
#endif
		//debug
#ifdef __DEBUG
		printf("order = %d", *order);
#endif

		// FIXME Magic number


#ifdef METIS
		MAXORDER = 3;
#else
		MAXORDER = 1;
#endif

		if (*order < 0 || *order > MAXORDER) {
			/* ------------------------
				 Multiple Minimum Degree
				 ------------------------ */

			ordmmd_(&logfil,            &neqns,             xadj2,
					adj2,               mat[*token].invp,   mat[*token].perm, 
					&iwsiz,             mat[*token].iwork,  &mat[*token].nnzl, 
					&mat[*token].nsub,  mat[*token].colcnt, &mat[*token].nsuper,
					mat[*token].xsuper, mat[*token].snodes, &sfiflg, 
					&iflag);

		}
		else if (*order == 1) {
			printf(" AMD not implemented yet, use MMD !\n");

			ordmmd_(&logfil,            &neqns,             xadj2,
					adj2,               mat[*token].invp,   mat[*token].perm,
					&iwsiz,             mat[*token].iwork,  &mat[*token].nnzl,
					&mat[*token].nsub,  mat[*token].colcnt, &mat[*token].nsuper,
					mat[*token].xsuper, mat[*token].snodes, &sfiflg,
					&iflag);
		}
#ifdef METIS
		else if (*order == 2) {

			/* ----------------------
				 Node Nested Dissection 
				 ---------------------- */

			numflag = 1;
			options[0] = 0;

			/* use default options now
				 options[1] = 3;
				 options[2] = 1;
				 options[3] = 2;
				 options[4] = 0;
				 options[5] = 3;
				 options[6] = 0;
				 options[7] = 1;
				 */

			METIS_NodeND(n,                 xadj2,   adj2, 
					&numflag,          options, mat[*token].perm, 
					mat[*token].invp);
		}
		else if (*order == 3) {
			/* ----------------------
				 Edge Nested Dissection
				 ---------------------- */
			numflag = 1;
			options[0] = 0;
			METIS_EdgeND(n,                 xadj2,    adj2,
					&numflag,          options,  mat[*token].perm, 
					mat[*token].invp);
		}
#endif
		else if (*order == 0) {
			/* use the input permutation vector */
			for (i = 0; i < neqns; i++) {
				mat[*token].perm[i] = perm[i];
				mat[*token].invp[perm[i]-1] = i+1;
			}
		}
		free(xadj2);
		free(adj2);

#ifdef TIMING
		t1 = getime();
		t1 = t1 - t0;
		printf("\n");

		if (*order < 0 || *order > 4) {
			printf(" MMD   :");
		}
		else if (*order == 1)
			printf(" AMD         :");
		else if (*order == 2) {
			printf(" NodeND      :");
		}
		else if (*order == 3) {
			printf(" EdgeND      :");
		}
		else if (*order == 0) {
			printf(" User defined:");
		}
		printf(" Time reordering    = %9.3e\n", t1);
		/*debug*/
		/*
			 fp = fopen("perm.vec","w");
			 for (i = 0; i < *n; i++)
			 fprintf(fp,"%d\n", mat[*token].perm[i]);
			 fclose(fp);
			 */
#endif

#ifdef TIMING
		t0 = getime();
#endif
		/* ----------------------
			 Symbolic factorization 
			 ---------------------- */
		if (*order >= 0 && *order <= 3) {
			/* not needed when MMD has been called */
			for (i=0;i<iwsiz;i++) mat[*token].iwork[i]=0;
			sfinit_(&logfil,              &neqns,        
					&nnza,                mat[*token].xadj, 
					mat[*token].adj,      mat[*token].perm,
					mat[*token].invp,     mat[*token].colcnt,
					&mat[*token].nnzl,    &mat[*token].nsub,
					&mat[*token].nsuper,  mat[*token].snodes,
					mat[*token].xsuper,   &iwsiz,
					mat[*token].iwork,    &iflag);
		}

		nsub   = mat[*token].nsub;
		nsuper = mat[*token].nsuper;
		nnzl   = mat[*token].nnzl;
		*Lnnz  = nnzl;

		mat[*token].xlindx = (int*)malloc((neqns+1)*sizeof(int));
		if (!mat[*token].xlindx) {
			fprintf(stderr,
					"ldlt_preprocess: memory allocation failed for xlindx\n");
			exit(1);
		}

		mat[*token].lindx = (int*)malloc(nsub*sizeof(int));
		if (!mat[*token].lindx) { 
			fprintf(stderr,
					"ldlt_preprocess: memory allocation failed for lindx\n");
			exit(1);
		} 

		mat[*token].xlnz = (int*)malloc((neqns+1)*sizeof(int));
		if (!mat[*token].xlnz) {
			fprintf(stderr,
					"ldlt_preprocess: memory allocation failed for xlnz\n"); 
			exit(1);
		}

		for (i=0; i<iwsiz; i++) mat[*token].iwork[i]=0;

		symfct_(&logfil,              &neqns,
				&nnza,                mat[*token].xadj, 
				mat[*token].adj,      mat[*token].perm,
				mat[*token].invp,     mat[*token].colcnt,
				&mat[*token].nsuper,  mat[*token].xsuper,
				mat[*token].snodes,   &mat[*token].nsub, 
				mat[*token].xlindx,   mat[*token].lindx,
				mat[*token].xlnz,     &iwsiz,
				mat[*token].iwork,    &iflag);

#ifdef TIMING
		t1 = getime();
		t1 = t1 - t0;
		printf(" Time Symbolic factorization      = %9.3e\n\n", t1);
#endif

		/* ---------------------------
			 Prepare for Numerical factorization 
			 and triangular solve.
			 --------------------------- */

		mat[*token].split = (int*)malloc(neqns*sizeof(int));
		if (!mat[*token].split) {
			fprintf(stderr,
					"memory allocation failed for split\n");
			exit(1);
		}


		bfinit_(&neqns,             &nsuper, 
				mat[*token].xsuper, mat[*token].snodes, 
				mat[*token].xlindx, mat[*token].lindx,  
				&cachsz,            &mat[*token].tmpsiz,
				mat[*token].split);
		ierr = 0;

	}
}

void ldlt_fact__(int *token, int *colptr, int *rowind, doublecomplex *nzvals)
{
	int   i, nnz, nnzl, neqns, nsuper, tmpsiz, logfil=6, iflag, iwsiz;
	int   ibegin, iend, irow, jcol;
#ifdef TIMING
	double t0, t1;
#endif
	/* Aug 30, 2009 */
	int   maxsup, jsup, nnzlplus, supsize;


	neqns  = mat[*token].n;
	nnz    = colptr[neqns]-1;
	nnzl   = mat[*token].nnzl;
	nsuper = mat[*token].nsuper;
	tmpsiz = mat[*token].tmpsiz;
	iwsiz  = 7*neqns+3;

#ifdef TIMING
	t0 = getime();
#endif
	if (mat[*token].anz == NULL) {
		mat[*token].anz = (doublecomplex*)malloc((2*nnz-neqns)*sizeof(doublecomplex));
		if (!mat[*token].anz) {
			fprintf(stderr,"memory allocation failed for anz\n");
			exit(1);
		} 
	}

	/* Aug 30, 2008, allocate extra space for diagal extraction 
		 the extra space is not used in the LDLT factorization */
	maxsup = 0;
	nnzlplus = nnzl;
	for (jsup = 0; jsup < nsuper; jsup++) {
		supsize = mat[*token].xsuper[jsup+1]-mat[*token].xsuper[jsup];
		if (supsize > maxsup) maxsup = supsize;
		nnzlplus = nnzlplus + supsize*(supsize-1)/2;
	}

#ifdef __DEBUG
	printf(" nnzlplus = %d\n", nnzlplus);
#endif

	if (mat[*token].lnz == NULL) {
		mat[*token].lnz = (doublecomplex*)malloc(nnzlplus*sizeof(doublecomplex));
		if (!mat[*token].lnz) {
			fprintf(stderr,"memory allocation failed for lnz.\n");
			exit(1);
		}
	}

	if (mat[*token].tmat == NULL) {
		mat[*token].tmat = (doublecomplex*)malloc(tmpsiz*sizeof(doublecomplex));
		if (!mat[*token].tmat) {
			fprintf(stderr,"memory allocation failed for tmat.\n");
			exit(1);
		}
	}

	if (mat[*token].diag == NULL) {  
		mat[*token].diag = (doublecomplex*)malloc(neqns*sizeof(doublecomplex));
		if (!mat[*token].diag) {
			fprintf(stderr,"memory allocation failed for diag.\n");
			exit(1);
		}
	}

	if (mat[*token].newrhs ==NULL) {
		mat[*token].newrhs = (doublecomplex*)malloc(neqns*sizeof(doublecomplex));
		if (!mat[*token].newrhs) {
			fprintf(stderr,"memory allocation failed for newrhs\n");
			exit(1);
		}
	}

	if (fullrep) {
		for (i=0; i<neqns+1; i++) mat[*token].xadj[i] = colptr[i];
		for (i=0; i<nnz; i++) mat[*token].adj[i] = rowind[i];
		for (i=0; i<nnz; i++) mat[*token].anz[i] = nzvals[i];
	}
	else {

		flo2ho_(&neqns,             colptr,           rowind,   
				nzvals,             mat[*token].xadj, mat[*token].adj, 
				mat[*token].anz,    mat[*token].iwork);

	}


	inpnv_(&logfil,            &neqns,              mat[*token].xadj,
			mat[*token].adj,    mat[*token].anz,     mat[*token].perm,   
			mat[*token].invp,   &mat[*token].nsuper, mat[*token].xsuper,
			mat[*token].xlindx, mat[*token].lindx,   mat[*token].xlnz,
			mat[*token].lnz,    &iwsiz,              mat[*token].iwork,
			&iflag);



	for (i=0; i<tmpsiz; i++) {
		(mat[*token].tmat[i]).r = 0.0;
		(mat[*token].tmat[i]).i = 0.0;
	}


	blkfct_(&logfil,            &neqns, 
			&nsuper,            &nunroll,
			mat[*token].xsuper, mat[*token].snodes, 
			mat[*token].split,  mat[*token].xlindx,
			mat[*token].lindx,  mat[*token].xlnz,
			mat[*token].lnz,    mat[*token].diag,
			&iwsiz,             mat[*token].iwork,
			&tmpsiz,            mat[*token].tmat,
			&iflag);


#ifdef TIMING
	t1 = getime();
	t1 = t1 - t0;
	printf(" Time numerical factorization      = %9.3e\n\n", t1);
#endif
}

void ldlt_solve__(int *token, doublecomplex *x, doublecomplex *rhs)
{
	int     neqns, nsuper, i;
	int     *perm, *invp;
	doublecomplex  *newrhs;

	neqns  = mat[*token].n;
	nsuper = mat[*token].nsuper;
	perm   = mat[*token].perm;
	invp   = mat[*token].invp;
	newrhs = mat[*token].newrhs;

	for (i=0; i<neqns; i++) newrhs[i] = rhs[perm[i]-1];

	blkslv_(&nsuper,              mat[*token].xsuper,   mat[*token].xlindx,
			mat[*token].lindx,    mat[*token].xlnz,   mat[*token].lnz,
			mat[*token].newrhs);

	for (i=0; i<neqns; i++) x[i] = newrhs[invp[i]-1];
}

void ldlt_getdiag__(int *token, doublecomplex *diag)
{
	int i, neqns;

	neqns = mat[*token].n;

	for (i=0;i<neqns;i++) diag[i] =  mat[*token].diag[i];
}


void ldlt_blkselinv__(int *token, int* colptr, int* rowind, 
		doublecomplex *inva, int *dumpL)
{
	int   i, nnzl, neqns, nsuper;
	doublecomplex  *ytemp, *plnz, *pinva;

#ifdef TIMING
	double t0, t1;
#endif
	neqns  = mat[*token].n;
	nnzl   = mat[*token].nnzl;
	nsuper = mat[*token].nsuper;

	ytemp = (doublecomplex*)calloc(neqns,sizeof(doublecomplex));
	selinvblk_(&nsuper,            mat[*token].xsuper,  mat[*token].xlindx,
			mat[*token].lindx,  mat[*token].xlnz  ,  mat[*token].lnz   ,
			mat[*token].snodes, mat[*token].perm  ,  &neqns            , 
			colptr            , rowind            ,  dumpL);
	// for(i = 0 ; i < nnzl; i++){
	// printf("lnz[i].r = %25.15f\n", mat[*token].lnz[i].r);
	// }

	/* Copy lnz back to inva */
	plnz = mat[*token].lnz;
	pinva = inva;
	for(i = 0; i < nnzl; i++){
		pinva->r = plnz->r;  pinva->i = plnz->i;
		pinva++; plnz++;
	}

	free(ytemp);
}


void ldlt_free__(int *token)
{
	free(mat[*token].xadj);
	free(mat[*token].adj);

	free(mat[*token].xlindx);
	free(mat[*token].lindx);
	free(mat[*token].xlnz);
	free(mat[*token].lnz);
	free(mat[*token].anz);
	free(mat[*token].diag);
	free(mat[*token].perm);
	free(mat[*token].invp);
	free(mat[*token].colcnt);
	free(mat[*token].xsuper);
	free(mat[*token].snodes);
	free(mat[*token].split);
	free(mat[*token].iwork);
	free(mat[*token].tmat);
	free(mat[*token].newrhs);

	mat[*token].xadj = NULL;
	mat[*token].adj = NULL;

	mat[*token].xlindx=NULL;
	mat[*token].lindx=NULL;
	mat[*token].xlnz=NULL;
	mat[*token].lnz=NULL;
	mat[*token].anz=NULL;
	mat[*token].diag=NULL;
	mat[*token].perm=NULL;
	mat[*token].invp=NULL;
	mat[*token].colcnt=NULL;
	mat[*token].xsuper=NULL;
	mat[*token].snodes=NULL;
	mat[*token].split=NULL;
	mat[*token].iwork=NULL;
	mat[*token].tmat=NULL;
	mat[*token].newrhs=NULL;
}
#endif


// *********************************************************************
// Implementation of the class SelInvInterface
// *********************************************************************

SelInvInterface::SelInvInterface	(  )
{
#ifndef _RELEASE_
	PushCallStack("SelInvInterface::SelInvInterface");
#endif
	isSelInvInitialized_ = false;
	// Only one matrix in the structure mat is used. If more than one
	// matrix is used, change parameter MAX_NUM_MAT
	token_ = 0; 

#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SelInvInterface::SelInvInterface  ----- 

SelInvInterface::~SelInvInterface	(  )
{
#ifndef _RELEASE_
	PushCallStack("SelInvInterface::~SelInvInterface");
#endif
	if( isSelInvInitialized_ )
		ldlt_free__(&token_);

#ifndef _RELEASE_
	PopCallStack();
#endif
} 		// -----  end of method SelInvInterface::~SelInvInterface  ----- 


void
SelInvInterface::SymbolicFactorize	( SparseMatrix<Scalar>& A, Int order, Int* perm, Int& Lnnz )
{
#ifndef _RELEASE_
	PushCallStack("SelInvInterface::SymbolicFactorize");
#endif
	ldlt_preprocess__( &token_, &A.size, A.colptr.Data(),
			A.rowind.Data(), &Lnnz, &order, perm );
	isSelInvInitialized_ = true;

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SelInvInterface::SymbolicFactorize  ----- 


void
SelInvInterface::NumericalFactorize	( SparseMatrix<Scalar>& A )
{
#ifndef _RELEASE_
	PushCallStack("SelInvInterface::NumericalFactorize");
#endif
#ifdef _USE_COMPLEX_
	ldlt_fact__( &token_, A.colptr.Data(), A.rowind.Data(), reinterpret_cast<doublecomplex*>(A.nzval.Data()) );
#endif

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SelInvInterface::NumericalFactorize  ----- 


void
SelInvInterface::Solve	( Scalar* x, Scalar* rhs )
{
#ifndef _RELEASE_
	PushCallStack("SelInvInterface::Solve");
#endif
#ifdef _USE_COMPLEX_
  ldlt_solve__( &token_, reinterpret_cast<doublecomplex*>(x), reinterpret_cast<doublecomplex*>(rhs) );
#endif

#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SelInvInterface::Solve  ----- 


void
SelInvInterface::SelInv	( SparseMatrix<Scalar>& Ainv )
{
#ifndef _RELEASE_
	PushCallStack("SelInvInterface::SelInv");
#endif
	Int dumpL = 0;
#ifdef _USE_COMPLEX_
	ldlt_blkselinv__( &token_, Ainv.colptr.Data(),
			Ainv.rowind.Data(), 
			reinterpret_cast<doublecomplex*>(Ainv.nzval.Data()), &dumpL);
#endif
#ifndef _RELEASE_
	PopCallStack();
#endif

	return ;
} 		// -----  end of method SelInvInterface::SelInv  ----- 
}
