#include "chem_pot.h"
#include "../src_pw/global.h"
#include "../src_pw/tools.h"

Chemical_Potential::Chemical_Potential(){}
Chemical_Potential::~Chemical_Potential(){}

// here mu is chemical potential
void Chemical_Potential::update_mu(const int &iter, const int &niter, double *mu0, const double *ne)
{
	TITLE("Chemical_Potential","update_mu");

	double min = -5;
	double max = 5;

	if(iter==0)
	{
		if( ucell.nelec > ne[iter]   ) 
		{
			mu0[iter+1] = mu0[iter] - 0.01;
		}
		else if(ucell.nelec < ne[iter] )
		{
			mu0[iter+1] = mu0[iter] + 0.01;
		}	
	}
	else
	{
		// mu: y
		// nu: x
		// new_y = y0 + (y1-y0)/(x1-x0) * (new_x-x0)
		if(ne[iter]==ne[iter-1])
		{
			mu0[iter+1]=mu0[iter];
			cout << " same ne ! " << endl;
			if( ucell.nelec > ne[iter]   ) 
			{
				mu0[iter+1] = mu0[iter] - 0.01;
			}
			else if(ucell.nelec < ne[iter] )
			{
				mu0[iter+1] = mu0[iter] + 0.01;
			}	
		}
		else
		{
			mu0[iter+1] = mu0[iter] + (mu0[iter]-mu0[iter-1])/(ne[iter]-ne[iter-1]) * (ucell.nelec-ne[iter]);
			if(mu0[iter+1] < min || mu0[iter+1] > max)
			{
				cout << " suppose to be " << mu0[iter+1] << endl;

				if( ne[iter] > ucell.nelec ) mu0[iter+1] = mu0[iter] - 0.01;
				if( ne[iter] < ucell.nelec ) mu0[iter+1] = mu0[iter] + 0.01;

			}
		}
	}
	return;
}
