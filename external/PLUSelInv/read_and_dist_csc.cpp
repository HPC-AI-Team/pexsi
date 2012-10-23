#include <math.h>
#include "superlu_ddefs.h"
#include "commoninc.hpp"

// Read and distribute a CSC format (binary) written by
// WriteSparseMatrix.m
//
// This subroutine only supports real arithmetic and symmetric matrix.

int read_and_dist_csc(SuperMatrix *A, int nrhs, double **rhs,
		int *ldb, double **x, int *ldx,
		FILE *fp, gridinfo_t *grid)
{
	SuperMatrix GA;              /* global A */
	double   *nzval_loc;         /* local */
	int_t    *colind_loc, *rowptr_loc, *rowptr;	 /* local */
	int_t    m, n, nnz, tmp, ncols;
	int_t    m_loc, fst_row, nnz_loc, nnz_loc_save, maxnnzloc;
	int_t    m_loc_fst; /* Record m_loc of the first p-1 processors,
												 when mod(m, p) is not zero. */ 
	int_t    iam, row, col, i, j, relpos;
	int_t    *indbuf, *m_loc_vec;
	int_t    size;
	double   *valbuf;
	char     trans[1];
	MPI_Status mpistatus;

	iam = grid->iam;

#if ( DEBUGlevel>=1 )
	CHECK_MALLOC(iam, "Enter dcreate_matrix()");
#endif

	if ( !iam ) {
		fread(&n, sizeof(int_t), 1, fp);
		fread(&nnz, sizeof(int_t), 1, fp);
	}

	MPI_Bcast( &n, 1, mpi_int_t, 0, grid->comm );
	m = n;
	MPI_Bcast( &nnz, 1, mpi_int_t, 0, grid->comm );

	// read row pointers
	rowptr = (int_t*)malloc((n+1)*sizeof(int_t));
	if ( !iam ) {
		fread(&size, sizeof(int_t), 1, fp);
		iA( size = n+1 );
		fread(rowptr, sizeof(int_t), n+1, fp);
	}
	MPI_Bcast( rowptr, n+1, mpi_int_t,  0, grid->comm );

	/* Compute the number of rows on each processor */
	m_loc_vec = (int_t*)malloc(grid->nprow * grid->npcol*sizeof(int_t));
	m_loc_fst = m / (grid->nprow * grid->npcol);
	for (i = 0; i < grid->nprow * grid->npcol; i++) {
		if (i < grid->nprow * grid->npcol-1 ) {
			m_loc_vec[i] = m_loc_fst;
		}
		else { 
			m_loc_vec[i] = m - m_loc_fst*(grid->nprow * grid->npcol-1);
		} 
	}
	m_loc = m_loc_vec[iam];

	rowptr_loc = (int_t*)intMalloc_dist((m_loc+1)); 

	/* construct local row pointer */
	for (i = 0; i < m_loc+1; i++)
		rowptr_loc[i] = rowptr[iam*m_loc_fst+i] - rowptr[iam*m_loc_fst]; 

	/* calculate nnz_loc on each processor */
	nnz_loc = 0;
	for (int irow = 0; irow < m_loc; irow++) {
		ncols = rowptr_loc[irow+1]-rowptr_loc[irow];
		nnz_loc = nnz_loc + ncols;
	}
	MPI_Allreduce(&nnz_loc,&maxnnzloc,1,mpi_int_t,MPI_MAX,grid->comm);

	colind_loc = (int_t*)intMalloc_dist(nnz_loc); 
	nzval_loc = (double*)doubleMalloc_dist(nnz_loc);

	// read and distribute column indices (has to be read by iam==0)
	if ( !iam ) {
		fread(&size, sizeof(int_t), 1, fp); iA( size = nnz );
		indbuf = (int_t*) calloc(maxnnzloc,sizeof(int_t));
		for (int ip = 0; ip < grid->nprow * grid->npcol; ip++) {
			nnz_loc_save = 0;
			for (int irow = ip*m_loc_fst; irow < ip*m_loc_fst+m_loc_vec[ip]; irow++) {
				ncols = rowptr[irow+1]-rowptr[irow];
				fread(&indbuf[nnz_loc_save], sizeof(int_t), ncols, fp);
				nnz_loc_save = nnz_loc_save + ncols;
			}
			if (ip > 0) {
				MPI_Send(&nnz_loc_save, 1, mpi_int_t, ip, 0, grid->comm);
				MPI_Send(indbuf, nnz_loc_save, mpi_int_t, ip, 1, grid->comm);
			}
			else {
				for (i = 0; i < nnz_loc_save; i++) colind_loc[i] = indbuf[i];
			}
		}
		free(indbuf);
	}
	else {
		MPI_Recv(&nnz_loc_save, 1, mpi_int_t, 0, 0, grid->comm, &mpistatus);
		MPI_Recv(colind_loc, nnz_loc_save, mpi_int_t, 0, 1, grid->comm, &mpistatus);
	}

	// read and distribute nzvals
	if ( !iam ) {
		fread(&size, sizeof(int_t), 1, fp); iA( size = nnz );
		valbuf = (double*) malloc(maxnnzloc*sizeof(double));
		for (int ip = 0; ip < grid->nprow * grid->npcol; ip++) {
			nnz_loc_save = 0;
			for (int irow = ip*m_loc_fst; irow < ip*m_loc_fst+m_loc_vec[ip]; irow++) {
				ncols = rowptr[irow+1]-rowptr[irow];
				fread(&valbuf[nnz_loc_save], sizeof(double), ncols, fp);
				nnz_loc_save = nnz_loc_save + ncols;
			}
			if (ip > 0) {
				MPI_Send(&nnz_loc_save, 1, mpi_int_t, ip, 0, grid->comm);
				MPI_Send(valbuf, nnz_loc_save, MPI_DOUBLE, ip, 1, grid->comm);
			}
			else {
				for (i = 0; i < nnz_loc_save; i++) nzval_loc[i] = valbuf[i];
			}
		}
		free(valbuf);
	}
	else {
		MPI_Recv(&nnz_loc_save, 1, mpi_int_t, 0, 0, grid->comm, &mpistatus);
		MPI_Recv(nzval_loc, nnz_loc_save, MPI_DOUBLE, 0, 1, grid->comm, &mpistatus);
	}

	/* Set up the local A in NR_loc format */
	fst_row = iam*m_loc_fst;
	dCreate_CompRowLoc_Matrix_dist(A, m, n, nnz_loc, m_loc, fst_row,
			nzval_loc, colind_loc, rowptr_loc,
			SLU_NR_loc, SLU_D, SLU_GE);

	/* Get the local B */
	if ( !((*rhs) = doubleMalloc_dist(m_loc*nrhs)) )
		ABORT("Malloc fails for rhs[]");
	for (j =0; j < nrhs; ++j) {
		for (i = 0; i < m_loc; ++i) {
			row = fst_row + i;
			(*rhs)[j*m_loc+i] = 1.0;
		}
	}
	*ldb = m_loc;

	/* Set the true X */    
	*ldx = m_loc;
	if ( !((*x) = doubleMalloc_dist(*ldx * nrhs)) )
		ABORT("Malloc fails for x_loc[]");

	/* Get the local part of xtrue_global */
	for (j = 0; j < nrhs; ++j) {
		for (i = 0; i < m_loc; ++i)
			(*x)[i + j*(*ldx)] = 1.0;
	}

	free(rowptr);
	free(m_loc_vec);

#if ( DEBUGlevel>=1 )
	printf("sizeof(NRforamt_loc) %d\n", sizeof(NRformat_loc));
	CHECK_MALLOC(iam, "Exit dcreate_matrix()");
#endif
	return 0;
}
