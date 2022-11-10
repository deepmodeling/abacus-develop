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
class Charge_Mixing
{
	public:
	Charge_Mixing();
	~Charge_Mixing();

    void set_mixing
    (
        const std::string &mixing_mode_in,
        const double &mixing_beta_in,
        const int &mixing_ndim_in,
		const double &mixing_gg0_in
    );//mohan add mixing_gg0_in 2014-09-27
	
    const std::string &get_mixing_mode() const {return mixing_mode;}
    double get_mixing_beta() const {return mixing_beta;}
    int get_mixing_ndim() const {return mixing_ndim;}
    double get_mixing_gg0() const {return mixing_gg0;}

	protected:

    //General parameters
    std::string mixing_mode;
    double mixing_beta;
    int mixing_ndim;
	double mixing_gg0; //mohan add 2014-09-27

	// simple plain mixing method.	
    void plain_mixing( double *rho_in, double *rho_save_in ) const;
	
	// tools
	double rhog_dot_product(const std::complex<double>*const*const rhog1, const std::complex<double>*const*const rhog2) const;		// Peize Lin add const 2019-05-01
	
};

#endif
