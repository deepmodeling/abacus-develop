//=================================
// Mohan add 2010-02-01
//=================================
#ifndef CHARGE_MIXING_H
#define CHARGE_MIXING_H

//==================================
// (1) Plain Mixing
// (2) KerKer Mixing
// (3) Pulay Mixing
// (4) Modified Broden Mixing
//===================================
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "charge.h"
class Charge_Mixing: public Charge
{
	public:
	Charge_Mixing();
	~Charge_Mixing();

    std::string mixing_mode;
    double mixing_beta;
    int mixing_ndim;
	double mixing_gg0; //mohan add 2014-09-27

    void set_mixing
    (
        const std::string &mixing_mode_in,
        const double &mixing_beta_in,
        const int &mixing_ndim_in,
	const double &mixing_gg0_in
    );//mohan add mixing_gg0_in 2014-09-27
	
	// also be used by namespace: GGA_pw
    void set_rhog(const double *rho_in, std::complex<double> *rhog_in)const;
	
	protected:

	// simple plain mixing method.	
    void plain_mixing( double *rho_in, double *rho_save_in ) const;
	
	// Kerker mixing method.
	// Kerker method must use plane wave basis, change the high energy plane wave 
	// more than low energy plane wave.
	void Kerker_mixing( double *rho, const std::complex<double> *residual_g, double *rho_save);
	
	// tools
	void set_rhor(std::complex<double> *rhog, double *rho)const;
	double rhog_dot_product(const std::complex<double>*const*const rhog1, const std::complex<double>*const*const rhog2) const;		// Peize Lin add const 2019-05-01
	
};

#endif
