#ifndef SELINV_H
#define SELINV_H

class Selinv
{
	public:

	Selinv();
	~Selinv();

	void using_SELINV(const int &ik, double *H, double *S);
	void out_ccf();

	static int Npole;
	static double temp;	
	static double gap;
	static double deltaE;
	static double mu;
	static double threshold;
	static int    niter;

	static double eband;
	static double nelec; 
	
	private:

};

#endif

