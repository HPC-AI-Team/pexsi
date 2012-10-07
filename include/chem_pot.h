#ifndef CHEM_POT_H
#define CHEM_POT_H

class Chemical_Potential
{
	public:
	Chemical_Potential();
	~Chemical_Potential();

	void update_mu(const int &iter, const int &niter, double *mu0, const double *ne);

	private:



};

#endif
