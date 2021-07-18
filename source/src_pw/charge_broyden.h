//=================================
// Mohan add 2010-02-05
//=================================
#ifndef CHARGE_BROYDEN_H
#define CHARGE_BROYDEN_H

//==================================
// (1) Plain Mixing
// (2) KerKer Mixing
// (3) Pulay Mixing
// (4) Modified Broden Mixing
//===================================
#include "tools.h"
#include "../module_base/matrix.h"
#include "charge_pulay.h"

class Charge_Broyden: public Charge_Pulay
{
	public:
	Charge_Broyden();
	~Charge_Broyden();

    void mix_rho(double &dr2,const double &tr2_min,
                 const double &tr2,const int &iter,
                 bool &converged);// mix rho

	private:

	// Sophisticated mixing method.
	void Modified_Broyden_mixing(void);
	void Simplified_Broyden_mixing(const int &iter); //qianrui created 2021-5-15
	void allocate_Broyden();

	void generate_beta(const int &is);
	void generate_Zmk(const int &totstep, const int &irstep, const int &idstep, const int &is);
	void generate_new_broyden_rho(const int &is, const int &irstep);

	bool initb; // b stands for Broyden algorithms.
	double w0;
	double* w;
	int broyden_type;
	matrix beta; // (dstep, dstep)
	matrix betabar; // (dstep, dstep)
	matrix* Zmk;
	matrix* Zmk_old;
};

#endif
