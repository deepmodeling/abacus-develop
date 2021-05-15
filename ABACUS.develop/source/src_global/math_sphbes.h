#ifndef MATH_SPHBES_H
#define MATH_SPHBES_H

//========================================================
// Spherical Bessel functions, mohan 2021-05-06
//========================================================
class Sphbes 
{
	public:

    Sphbes();
    ~Sphbes();

    static void Spherical_Bessel
    (
        const int &msh,	//number of grid points
        const double *r,//radial grid
        const double &q,	//
        const int &l,	//angular momentum
        double *jl	//jl(1:msh) = j_l(q*r(i)),spherical bessel function
    );

	static void Spherical_Bessel
	(           
	    const int &msh, //number of grid points
		const double *r,//radial grid
		const double &q,    //
		const int &l,   //angular momentum
		double *sj,     //jl(1:msh) = j_l(q*r(i)),spherical bessel function
		double *sjp
	);

	
    static void Spherical_Bessel_Roots
    (
        const int &num,
        const int &l,
        const double &epsilon,
        double* eigenvalue,
        const double &rcut
    );

private:

    static double Spherical_Bessel_7(const int n, const double &x);

	// Peize Lin change double to void 2019-05-01
    static void BESSJY(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp);

    static void BESCHB(double x, double *gam1, double *gam2, double *gampl, double *gammi);

    static double CHEBEV(double a, double b, double c[], int m, double x);

    static int IMAX(int a, int b);
};

#endif
