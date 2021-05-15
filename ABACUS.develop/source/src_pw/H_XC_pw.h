#ifndef H_XC_PW_H
#define H_XC_PW_H

#include "tools.h"

class H_XC_pw 
{
	public:

	friend class Stress_Func;
	friend class Stress_PW;
	friend class Forces;
    friend class Force_Stress_LCAO; 
	friend class Potential;
	friend class energy;
	friend class eximport;

    H_XC_pw();
    ~H_XC_pw();

	private:

	// the Hartree energy
    static double etxc;
	static double vtxc;

	// compute the exchange-correlation energy 
	// [etxc, vtxc, v] = v_xc(...)
    static std::tuple<double,double,matrix> v_xc(
		const int &nrxx, // number of real-space grid
		const int &ncxyz, // total number of charge grid
		const double &omega, // volume of cell
		const double*const*const rho_in, 
		const double*const rho_core); // core charge density

};

#endif //Exchange-correlation energy
