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
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "charge_pulay.h"

class Charge_Broyden: public Charge_Pulay
{
	public:
	Charge_Broyden();
	~Charge_Broyden();
	double get_drho(double** rho, double** rho_save,
		std::complex<double>** rhog, std::complex<double>** rhog_save, const double nelec);

    void mix_rho(const int &iter,
		double** rho,
		double** rho_save,
		std::complex<double>** rhog,
		std::complex<double>** rhog_save
	);// mix rho

	private:

	// Sophisticated mixing method.
	void Simplified_Broyden_mixing(const int &iter,
		double** rho,
		double** rho_save,
		std::complex<double>** rhog,
		std::complex<double>** rhog_save); //qianrui created 2021-5-15
	void allocate_Broyden();

	bool initb; // b stands for Broyden algorithms.
	ModuleBase::matrix beta; // (dstep, dstep)
};

#endif
