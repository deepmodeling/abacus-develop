//==========================================================
// AUTHOR : mohan
// DATE : 2009-3-29
// Last Modify : 2021-01-04
//==========================================================
#ifndef BESSEL_BASIS_H
#define BESSEL_BASIS_H
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/realarray.h"

//==========================================================
// CLASS :
// NAME :  Bessel_Basis
//==========================================================
class Bessel_Basis
{
public:
	Bessel_Basis();
	~Bessel_Basis();

	//--------------------------------------------------------------------------
	// used for a specific group of C4 coefficients.
	// mohan updated 2021-01-04, 
	// mohan added a new input parameter lmax_in,
	// if we only generate numerical atomic orbitals based on
	// spherical Bessel functions, lmax_in = ucell.lmax
	// However, if we want to generate spherical Bessel functions for descriptor,
	// then the lmax_in is controlled by user.
	//--------------------------------------------------------------------------
	void init( 
		const bool start_from_file,
		const double &ecutwfc,
		const int &ntype,
		const int &lmax_in,
		const double &dk = 0.01,
		const double &dr = 0.01);

	const int& get_ecut_number() const { return Ecut_number;}

	// get value from interpolation
	double Polynomial_Interpolation(const int &it, const int &l, const int &ic, const double &gnorm)const;
	double Polynomial_Interpolation2(const int &l, const int &ie, const double &gnorm)const;

	// get ecut : get energy cut off, which is used to truncate spherical Bessel function Jlq. 
	// get rcut : cut off radius(angstrom) of radial spherical Bessel function Jlq.
	const double &get_ecut(void) const {return ecut;}
	const double &get_rcut(void) const {return rcut;}
	const double &get_tolerence(void) const {return tolerence;}

	// get_smooth and get_sigma. mohan add 2009-08-28
	// in this case, the Jlq are not the true Jlq.
	const bool &get_smooth(void) const {return smooth;}
	const double &get_sigma(void) const {return sigma;}

private:
	// the most important array to calculate spillage
	ModuleBase::realArray Faln;// (ntype, max_l+1, max_n, ik)

	// Coefficients to be optimized!
	ModuleBase::realArray C4;

	ModuleBase::realArray TableOne;
	int kmesh;
	double Dk;
	int Ecut_number;
	double rcut;
	double ecut;
	double tolerence;

	bool smooth; // mohan add 2009-01-18
	double sigma;

	//========================================================
	// MEMBER FUNCTIONS :
	// NAME : readin ( read in parameters )
	//========================================================
	void readin(const std::string &name); 
	void bcast(void);

	void allocate_C4(
		const int &ntype,
		const int &lmax, 
		const int &nmax,
		const int &ecut_number);

	void readin_C4(
		const std::string &name,
		const int &ntype,
		const int &ecut,
		const int &rcut,
		const int &ecut_number,
		const double &tolerence);

	void init_TableOne(void);

	void init_Faln( const int &ntype,const int &lmax,const int &nmax,const int &ecut_number);
	int nwfc; // number of localized wave functions

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
