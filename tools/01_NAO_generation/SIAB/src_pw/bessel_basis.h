//==========================================================
// AUTHOR : mohan
// DATE : 2009-3-29
// Last Modify : 2009-08-28
//==========================================================
#ifndef BESSEL_BASIS_H
#define BESSEL_BASIS_H
#include "../src_spillage/common.h"

//==========================================================
// CLASS :
// NAME :  Bessel_Basis
//==========================================================
class Bessel_Basis
{
public:
	Bessel_Basis();
	~Bessel_Basis();

	// used for a specific group of C4 coefficients.
	void init( 
		const double &dk = 0.01,
		const double &dr = 0.01);

	// get value from interpolation
	double Polynomial_Interpolation2(const int &l, const int &ie, const double &gnorm)const;

	realArray TableOne;

private:
	int kmesh;
	double Dk;

	// init table, used for output overlap Q.
	void init_TableOne(
		const bool smooth_in,
		const double &sigma_in,
		const double &ecut,
		const double &rcut,
		const double &dr,
		const double &dk,
		const int &lmax,
		const int &ecut_number,
		const double &tolerence);
};

#endif
