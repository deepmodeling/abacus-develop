#ifndef MATHZONE_ADD1_H
#define MATHZONE_ADD1_H

#include <cmath>
#include <cassert>

using namespace std;

/************************************************************************
LiaoChen add @ 2010/03/09 to add efficient functions in LCAO calculation
************************************************************************/
class Mathzone_Add1
{
public:

    Mathzone_Add1();
    ~Mathzone_Add1();

	static void Sph_Bes (double x, int lmax, double *sb, double *dsb); 

	//calculate jlx and its derivatives
	static void Spherical_Bessel
	(           
	    const int &msh, //number of grid points
		const double *r,//radial grid
		const double &q,    //
		const int &l,   //angular momentum
		double *sj,     //jl(1:msh) = j_l(q*r(i)),spherical bessel function
		double *sjp
	);
	
	/***********************************************
	function : integrate function within [a,b]
	formula : \int func(r) dr
	***********************************************/
	static double uni_simpson (const double* func, const int& mshr, const double& dr);
	
	/********************************************************
	function : radial fourier transform
	formula : F(k)_{n,l} = \int r^{n} dr jl(kr) \phi(r)
	input : 
			int n (exponential index)
			int aml (angular momentum)
			int mshk (number of kpoints)
			double* arr_k (an array of kpoints)
			int mshr (number of grids)
			double* ri (radial points)
			double dr (delta r)
			double* phir (function to be transformed)
	output :
			double* phik (with size of mshk)
	********************************************************/
	static void uni_radfft
	(
		const int& n,
		const int& aml,
		const int& mshk,
		const double* arr_k,
		const int& mshr,
		const double* ri,
		const double& dr,
		const double* phir,
		double* phik
	);

	static void Sbt_new 
	(
		const int& polint_order,
		const int& l,
		const double* k,
		const double& dk,
		const int& mshk,
		const double* r,
		const double& dr,
		const int& mshr,
		const double* fr,
		const int& rpow,
		double* fk
	);

	static double dualfac (const int& l);
	static double factorial (const int& l);
	static double fourier_cosine_transform
	(
		const double* func,
		const double* r,
		const int& mshr,
		const double& dr,
		const double& k
	);
	
	static double fourier_sine_transform
	(
		const double* func,
		const double* r,
		const int& mshr,
		const double& dr,
		const double& k
	);
	
/*
	static void RecursiveRadfft
	(
		const int& l,
		const double** sb,
		const double** dsb,
		const double* arr_k,
		const double* kpoint,
		const int& kmesh,
		const double& dk,
		double& rs,
		double& drs
	);	
*/

	static double Polynomial_Interpolation
	(
 		const double* xa,
		const double* ya,
		const int& n,
		const double& x
	);
	
	static void SplineD2
	(
 		const double *rad,
		const double *rad_f,
		const int& mesh,
		const double &yp1,
		const double &ypn,
		double* y2
	);

	static void Cubic_Spline_Interpolation
	(
		const double * const rad,
		const double * const rad_f,
		const double * const y2,
		const int& mesh,
		const double * const r,
		const int& rsize,
		double * const y,
		double * const dy
	);
	
	static double RadialF
	(
		const double* rad,
		const double* rad_f,
		const int& msh,
		const int& l,
		const double& R
	);
	
	static double Uni_RadialF
	(
		const double* rad_f,
		const int& msh,
		const double& dr,
		const double& R
	);
	
	static void test ();		
	static void test2 ();
	
	
	static void Uni_Deriv_Phi
	(
		const double *radf,
		const int &mesh,
		const double &dr,
		const int &nd,
		double* phind
	);
	
	private:
	const static int sph_lmax = 20;
	static double** c_ln_c;
	static double** c_ln_s;
	static bool flag_jlx_expand_coef;

	static void expand_coef_jlx ();
	static double pol_seg_int
		(
		 const int& polint_order,
		 const double* coef,
		 const int& n,
		 const double* k,
		 const int& ik
		);
};

#endif
