#ifndef CHARGE_EXTRA_H
#define CHARGE_EXTRA_H

using namespace std;

//--------------------------------
// charge extrapolation method
// written by Xiaohui 2014
//--------------------------------
class Charge_Extra
{
	public:

	Charge_Extra();
	~Charge_Extra();

	void allocate_ions(void);
	void extrapolate_charge(void);

	// use "istep = ions.istep"
	int istep;

	// for the second-order extrapolation
	double* pos_old1;
	double* pos_old2;
	double* pos_now;
	double* pos_next;

	private:

	double*** rho_ion; //(dim, nspin, pw.nrxx)
	bool init_rho;
	int dim;

	// use the first-order extrapolation
	double** delta_rho1;
	double** delta_rho2;
	double** delta_rho;

	// use the second-order extrapolation
	double** delta_rho3;
	int pos_dim;
	double alpha,beta;

	void find_alpha_and_beta(void);

};

#endif
