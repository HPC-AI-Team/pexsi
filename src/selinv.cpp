#include "selinv.h"
#include "../src_pw/global.h"
#include "../src_pw/tools.h"
#include "cal_reduce_rho.h"
#include "../src_algorithms/hs_matrix.h"
#include "trace_rho_hs.h"
#include "chem_pot.h"

int Selinv::Npole=40;
double Selinv::temp=2000.0;	
double Selinv::gap=0.0;
double Selinv::deltaE=0.0;
double Selinv::mu=-1.0;
double Selinv::threshold=1.0e-3;
int    Selinv::niter=30;

double Selinv::eband;
double Selinv::nelec;

Selinv::Selinv(){}
Selinv::~Selinv(){}

void Selinv::using_SELINV(const int &ik, double *H, double *S)
{
	TITLE("Void","using_SELINV");
	Trace_Rho_HS TRHS;

//#ifdef __SELINV
	
	// nonzero need to be update later.
	bool** nonzero = new bool*[NLOCAL];
	for(int col=0; col<NLOCAL; ++col)
	{
		nonzero[col] = new bool[NLOCAL-col];
		for(int row=col; row<NLOCAL; ++row)
		{
			int index = row-col;
			nonzero[col][index] = false;
			//nonzero[col][index] = true;
		}
	}


	int Hnnz=TRHS.cal_Hnnz(nonzero);
	//int Hnnz=NLOCAL*(NLOCAL+1)/2;
	
	int* colptr_H=new int[NLOCAL+1];//column element counter
	int* rowind_H=new int[Hnnz];//row index for each guy in column
	double* nzval_H=new double[Hnnz];//values of H
	double* nzval_S=new double[Hnnz];//values of S
	double* nzval_rho=new double[Hnnz];//value of rho
	ZEROS(colptr_H,NLOCAL+1);
	ZEROS(rowind_H,Hnnz);
	ZEROS(nzval_H,Hnnz);
	ZEROS(nzval_S,Hnnz);

	const double mem1=Memory::record("Selinv","nonzero",NLOCAL*(NLOCAL+1.0)/2.0,"bool");
	const double mem2=Memory::record("Selinv","H_S_Rho",Hnnz*3,"double");
	OUT(ofs_running,"mem_nonzero (MB)",mem1);
	OUT(ofs_running,"mem_HSRho (MB)",mem2);
	
	TRHS.transfer_HS(Hnnz, nonzero, colptr_H, rowind_H, H, S, nzval_H, nzval_S);

	for(int col=0; col<NLOCAL; ++col)
	{
		delete[] nonzero[col];
	}
	delete[] nonzero;
	
	if(Npole<=0 || Npole>1000)
	{
		ofs_running << " Npole=" << Npole << endl;
		WARNING_QUIT("Selinv","using_SELINV");
	}
	if(Npole%2!=0)
	{
		ofs_running << " Npole=" << Npole << endl;
		WARNING_QUIT("Selinv","Npole must be an even number!");
	}
	
	assert(deltaE>0.0);
	assert(temp>0.0);

	OUT(ofs_running,"Hnnz",Hnnz);
	OUT(ofs_running,"Npole",Npole);
	OUT(ofs_running,"temp",temp);
	OUT(ofs_running,"gap",gap);
	OUT(ofs_running,"deltaE",deltaE);
	OUT(ofs_running,"mu",mu);
	OUT(ofs_running,"threshold",threshold);
	OUT(ofs_running,"niter",niter);
	OUT(ofs_running,"NLOCAL*(NLOCAL+1)/2",NLOCAL*(NLOCAL+1.0)/2.0);
	OUT(ofs_running,"NLOCAL*NLOCAL",NLOCAL*NLOCAL);
	

	double *mu0 = new double[niter+1];
	double *ne = new double[niter];
	ZEROS(mu0, niter+1);
	ZEROS(ne, niter);
	mu0[0]=mu;
	ofs_running << " " << setw(5) << "Iter" << setw(15) << "Ne_now" 
		<< setw(15) << "Nelec" << setw(20) << "mu_now" 
		<< setw(20) << "mu_next" << setw(20) << "time_selinv" << setw(20) << "time_reduce" << endl; 
	ofs_running << setiosflags(ios::scientific);
	ofs_running << setprecision(6);
	bool converged = false;

	double start1, start2;
	double end1, end2;

	int pole1=0;
	int pole2=Npole;
#ifdef __MPI
	int* distri = new int[NPROC];
	ZEROS(distri, NPROC);
	const int base = Npole/NPROC; //average first
	const int remain = Npole - base * NPROC;//distribute the remain pole
	for(int i=0; i<NPROC; ++i) distri[i]=base;
	for(int i=0; i<remain; ++i) ++distri[i];
	for(int i=0; i<MY_RANK; ++i)//count the start pole
	{
		pole1+=distri[i];
	}	
	pole2 = pole1 + distri[MY_RANK];//the end pole
#endif
	OUT(ofs_running,"Start_Pole",pole1+1);//start from 1
	OUT(ofs_running,"End_Pole",pole2);


	const double threshold_tot = threshold * ucell.nelec;
	for(int iter=0; iter<niter; ++iter)
	{
		time_t time_start = time(NULL);

#ifdef __MPI
		start1 = MPI_Wtime();
#endif
	
		Cal_Reduce_Rho CRR;
		CRR.init(Npole, temp, gap, deltaE, mu0[iter], NLOCAL,
				colptr_H, rowind_H, nzval_H, nzval_S, nzval_rho, pole1, pole2);

#ifdef __MPI
		end1 = MPI_Wtime();
		start2 = MPI_Wtime();
		Parallel_Reduce::reduce_double_all( nzval_rho, Hnnz );	
#endif

#ifdef __MPI
		end2 = MPI_Wtime();
#endif

		if(GAMMA_ONLY_LOCAL)
		{
//			ne[iter]=UHM.GG.cal_rho(); // this can also check the electron number, but make sure rho is zero at first.
			ne[iter]=TRHS.do_it(colptr_H, rowind_H, nzval_rho, nzval_S);
			this->nelec = ne[iter];
		}
	
		// check if the convergence has been achieved.
		if( abs(ucell.nelec-ne[iter]) < threshold_tot)
		{
			mu = mu0[iter];
			
			ofs_running << " " << setw(5) << iter+1 << setw(15) << ne[iter] 
			<< setw(15) << ucell.nelec << setw(20) << mu0[iter] << setw(20) << mu0[iter]
			<< setw(20) << end1-start1 << setw(20) << end2-start2 << setw(20) << "DONE." << endl;
			converged = true;

			bool bit = false;
			if(ParaO.out_hs==3)
			{
				if(MY_RANK==0)
				HS_Matrix::save_HS_ccf(iter, Hnnz, colptr_H, rowind_H, nzval_H, nzval_S, bit);
			}

//			TRHS.transfer_rho(iter,colptr_H, rowind_H, nzval_rho);

			time_t time_end = time(NULL);
			OUT_TIME("selinv time", time_start, time_end);

			break;
		}

		// update the chemical potential
		Chemical_Potential CP;
		CP.update_mu(iter, niter, mu0, ne);

		ofs_running << " " << setw(5) << iter+1 << setw(15) << ne[iter] 
		<< setw(15) << ucell.nelec << setw(20) << mu0[iter] << setw(20) << mu0[iter+1] 
		<< setw(20) << end1-start1 << setw(20) << end2-start2 << endl;

		// output information
		//------------------------------------------------------------------
		// save H and S matrix to disk.
		bool bit = false;
		if(ParaO.out_hs==3)
		{
			if(MY_RANK==0)
			HS_Matrix::save_HS_ccf(iter, Hnnz, colptr_H, rowind_H, nzval_H, nzval_S, bit);
		}
		//	TRHS.transfer_rho(iter,colptr_H, rowind_H, nzval_rho);
		//------------------------------------------------------------------
	
		// time counting
		time_t time_end = time(NULL);
		OUT_TIME("selinv time", time_start, time_end);

	}


#ifdef __MPI
	delete[] distri;
#endif


	if(!converged)
	{
		ofs_warning << "ne not converged!" << endl;
		ofs_warning << "may be you should stop the run." << endl;
	}

	this->eband = TRHS.do_it(colptr_H, rowind_H, nzval_rho, nzval_H);
//	this->nelec = TRHS.do_it(colptr_H, rowind_H, nzval_rho, nzval_S);
	OUT(ofs_running, "Eband from trace[H rho]", this->eband);
	OUT(ofs_running, "nelec from trace[S rho]", this->nelec);

	TRHS.transfer_rho(0,colptr_H, rowind_H, nzval_rho);

	delete[] mu0;
	delete[] ne;

	delete[] colptr_H;
	delete[] rowind_H;
	delete[] nzval_H;
	delete[] nzval_S;
	delete[] nzval_rho;

//#else
//	WARNING_QUIT("Selinv","using_SELINV need to re compile the code with __SELINV.");
//#endif
	return;
}

