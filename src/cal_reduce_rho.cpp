#include "cal_reduce_rho.h"
#include "../src_pw/tools.h"
#include "siao_interf.h"
#include "getpolef.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

// LLIN: Macros below are to be compatible with the FORTRAN convention.
#define rowind_H(i) rowind_H[(i)-1]
#define colptr_H(i) colptr_H[(i)-1]
#define nzval_H(i)  nzval_H[(i)-1]
#define nzval_S(i)  nzval_S[(i)-1]
#define nzval_A(i)  nzval_A[(i)-1]
#define nzval_invA(i)  nzval_invA[(i)-1]
#define rowind_invA(i) rowind_invA[(i)-1]
#define colptr_invA(i) colptr_invA[(i)-1]
#define nzval_rho(i)  nzval_rho[(i)-1]

Cal_Reduce_Rho::Cal_Reduce_Rho(){}
Cal_Reduce_Rho::~Cal_Reduce_Rho(){}

/*********************************************************************
  CALREDUCERHO calculates the selected elements of the reduced density 
  matrix. 

Input:
Npole     :    the number of poles to be used.
temp      :    temperature, unit(K)
Gap       :    Energy gap defined to be min(abs(EV-mu)).
EV is the eigenvalue set of Hamiltonian,
unit(hatree) 
deltaE    :    Spectrum width defined to be
max(EV)-min(EV). EV is the eigenvalue set
of Hamiltonian, unit(hartree) 
mu        :    Chemical potential, unit(hartree)
Nnodes    :    The matrix dimension of H, S, Rho
colptr_H  :    Column partition of H, S, Rho
rowind_H  :    Row partition of H, S, Rho
nzval_H   :    Non-zero values of H.
nzval_S   :    Non-zero values of S.

Output:

nzval_rho :    Values of Rho following sparsity pattern 
(colptr_H, rowind_H)

Author:
Lin Lin and Chao Yang
Computer Research Division, Lawrence Berkeley National Lab
Last modified:  09-17-2011
 */

//#ifdef __SELINV
void Cal_Reduce_Rho::init(int &Npole, double &Temp, 
		double &Gap, double &DeltaE, 
		double &mu, int Nnodes, int* colptr_H, int* rowind_H,
		double* nzval_H, double* nzval_S, double* nzval_rho,
		const int &pole1, const int &pole2)
{
	TITLE("Cal_Reduce_Rho","init");
	// Assumes that H and S have the same sparsity pattern
	doublecomplex *zshift, *zweight;
	int token, Lnnz, Hnnz, order;
	int* perm;
	doublecomplex *nzval_invA, *nzval_A;
	int *colptr_invA, *rowind_invA;
	int ii, j, k, l;
	FILE *fid;
	int dumpL = 0;

	zshift  = (doublecomplex*) malloc(Npole*sizeof(doublecomplex));
	zweight = (doublecomplex*) malloc(Npole*sizeof(doublecomplex));
	// Pole expansion
	// zshift:
	// zweight:
	// Npole:
	// Temp:
	// Gap:
	// DeltaE:
	// mu:

	getpole_rho(zshift, zweight, &Npole, &Temp, &Gap, &DeltaE, &mu);

	//LLIN: complex.h incompatible with this code
//	Get_Pole GP;
//	GP.init_rho(zshift, zweight, &Npole, &Temp, &Gap, &DeltaE, &mu);


	token = 0;
	order = -1;

//	for(int i=0; i<Npole; ++i)
//	{
//		printf("i=%d,Zshift=%f,%f\n",i,zshift[i].r,zshift[i].i);
//	}
	
	perm = NULL;   // perm is not used here

	SIAO_Interface::ldlt_preprocess__(&token, &Nnodes, colptr_H, rowind_H, &Lnnz, &order, perm);   

	// Assume that the space for nzval_rho has been assigned outside this
	// subroutine.
	Hnnz        = colptr_H(Nnodes+1) - 1;
	nzval_A     = (doublecomplex*) malloc(sizeof(doublecomplex)*Hnnz);

	colptr_invA = (int*) malloc(sizeof(int)*(Nnodes+1));
	rowind_invA = (int*) malloc(sizeof(int)*Lnnz);
	nzval_invA  = (doublecomplex*) malloc(sizeof(doublecomplex)*Lnnz);

	// LLIN: Note the size of nzval_rho!
	for(ii = 0; ii < Hnnz; ii++)
	{
		nzval_rho[ii] = 0.0;
	}

	// Loop over all the complex poles.
//	for(l = 0; l < Npole; l++)
	for(l = pole1; l < pole2; l++) //mohan update 2011-10-03
	{
		for(ii = 0; ii < Hnnz; ii++)
		{
			nzval_A[ii].r = nzval_H[ii] - zshift[l].r * nzval_S[ii];
			nzval_A[ii].i = -zshift[l].i * nzval_S[ii];
		}


//		if(0)
//		{
//			fid = fopen("A", "w");
//			for(j = 1; j < Nnodes+1; j++)
//			{
//				for(i = colptr_H(j); i < colptr_H(j+1); i++)
//				{
//					fprintf(fid, "%5d %5d %25.15e %25.15e\n",
//							rowind_H(i), j, nzval_A(i).r, nzval_A(i).i);
//				}
//			}
//			fclose(fid);
//		}


		// factorization 
		SIAO_Interface::ldlt_fact__(&token, colptr_H, rowind_H, nzval_A);

		// selected inversion 
		SIAO_Interface::ldlt_blkselinv__(&token, colptr_invA, rowind_invA, 
				nzval_invA, &dumpL);

		// rho = rho + imag(zweight*invA), applied only for the matrix
		//  * elements of rho corresponding to the nonzero elements of H  
		for(j = 1; j < Nnodes + 1; j++)
		{
			for(ii = colptr_H(j); ii < colptr_H(j+1); ii++)
			{
				for( k = colptr_invA(j); k < colptr_invA(j+1); k++)
				{
					if( rowind_H(ii) == rowind_invA(k) )
					{
						nzval_rho(ii) += zweight[l].r * nzval_invA(k).i + 
							zweight[l].i * nzval_invA(k).r;
						break;
					}
				}
				if( k == colptr_invA(j+1) ) 
				{
					fprintf(stderr, "cannot find the index!\n");
				}
			}// end ii
		} // for (j)
	} // end l

	

	SIAO_Interface::ldlt_free__(&token); 

	if (order == 0) free(perm);

	free(zshift);
	free(zweight);
	free(nzval_invA);
	free(nzval_A);
	free(colptr_invA);
	free(rowind_invA);
}
//#endif
