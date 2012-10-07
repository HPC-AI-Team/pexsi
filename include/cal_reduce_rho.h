#ifndef CAL_REDUCE_RHO_H
#define CAL_REDUCE_RHO_H

#include "getpolef.h"

class Cal_Reduce_Rho
{
	public:
	Cal_Reduce_Rho();
	~Cal_Reduce_Rho();
	
	void init(int &Npole, double &Temp, 
		double &Gap, double &DeltaE, 
		double &mu, int Nnodes, int* colptr_H, int* rowind_H,
		double* nzval_H, double* nzval_S, double* nzval_rho, const int &pole1, const int &pole2);

};

#endif
