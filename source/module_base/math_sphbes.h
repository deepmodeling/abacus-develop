#ifndef MATH_SPHBES_H
#define MATH_SPHBES_H

//========================================================
// Spherical Bessel functions, mohan 2021-05-06
//========================================================
namespace ModuleBase
{

class Sphbes 
{
	public:

    Sphbes();
    ~Sphbes();

    /**
     * @brief spherical bessel jl(qr)
     * 
     * @param msh [in] number of grid points
     * @param r [in] radial grid
     * @param q [in] k_radial
     * @param l [in] angular momentum
     * @param jl [out] jl spherical bessel function
     */
    static void Spherical_Bessel
    (
        const int &msh,
        const double *r,
        const double &q,
        const int &l,
        double *jl
    );

    /**
     * @brief derivative of spherical bessel djl(qr)/d(qr)
     * 
     * @param msh [in] number of grid points
     * @param r [in] radial grid
     * @param q [in] k_radial
     * @param l [in] angular momentum
     * @param jl [out] jl spherical bessel function
     */
    static void dSpherical_Bessel_dx
    (
        const int &msh,
        const double *r,
        const double &q,
        const int &l,
        double *jl
    );

    /**
     * @brief spherical bessel
     * 
     * @param msh [in] number of grid points
     * @param r [in] radial grid
     * @param q [in] k_radial
     * @param l [in] angular momentum
     * @param jl [out] jl spherical bessel function
     * @param sjp [out] sjp[i] is assigned to be 1.0. i < msh.
     */
	static void Spherical_Bessel
	(           
	    const int &msh, 
		const double *r,
		const double &q, 
		const int &l,  
		double *sj,    
		double *sjp
	);

	/**
	 * @brief return num eigenvalues of spherical bessel function
	 * 
	 * @param num [in] the number of eigenvalues
	 * @param l [in] angular number
	 * @param epsilon [in] the accuracy 
	 * @param eigenvalue [out] the calculated eigenvalues
	 * @param rcut [in] the cutoff the radial function
	 */
    static void Spherical_Bessel_Roots
    (
        const int &num,
        const int &l,
        const double &epsilon,
        double* eigenvalue,
        const double &rcut
    );

    //! spherical Bessel function of the first kind
    /*!
     *  This function computes j_l(x) by series expansion for x < l
     *  and by ascending recurrence for x >= l.
     *                                                              */
    static double sphbesj(const int l,   //!< [in] order
                          const double x //!< [in] argument
    );

    static void sphbesj(const int n,
                        const double* const r,
                        const double q,
                        const int l,
                        double* const jl
    );

private:

    static double Spherical_Bessel_7(const int n, const double &x);

	// Peize Lin change double to void 2019-05-01
    static void BESSJY(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp);

    static void BESCHB(double x, double *gam1, double *gam2, double *gampl, double *gammi);

    static double CHEBEV(double a, double b, double c[], int m, double x);

    static int IMAX(int a, int b);

    // utility functions for sphbesj
    static double _sphbesj_ascending_recurrence(int l, double x);
    static double _sphbesj_series(int l, double x);
};

}

#endif
