#ifndef MATH_INTEGRAL_H
#define MATH_INTEGRAL_H

// mohan add 2021-04-03
namespace ModuleBase
{

class Integral
{

    public:

    Integral();
    ~Integral();

    /**
     * @brief Simpson integration with a radial Jacobian.
     * 
     * @param mesh [in] positive odd number of grid points; mesh=1 returns zero
     * @param func [in] function to be integrated
     * @param rab [in] rab[i]=(dr/di)*Delta_i at each grid point. The index spacing Delta_i is 1.
     * @param asum [out] approximation to \int_{r[0]}^{r[mesh-1]} f(r) dr
     * @author Peize Lin
     * @date 2017-10-02
     */
    static void Simpson_Integral 
    (    
        const int mesh,
        const double * const func,
        const double * const rab,
        double &asum
    );

    /**
     * @brief Simpson integration on a uniform grid.
     * 
     * @param mesh [in] positive odd number of grid points; mesh=1 returns zero
     * @param func [in] function to be integrated
     * @param dr [in] interval
     * @param asum [out] approximation to \int_{r[0]}^{r[mesh-1]} f(r) dr
     * @author Peize Lin
     * @date 2017-10-02
     */
    static void Simpson_Integral 
    (
        const int mesh,
        const double * const func,
        const double dr,
        double &asum
    );

    /**
     * @brief Cumulative integration from the first grid point.
     *
     * @param mesh [in] number of grid points, mesh >= 1 (odd or even)
     * @param func [in] function to be integrated
     * @param rab [in] rab[i]=(dr/di)*Delta_i at each grid point. The index spacing Delta_i is 1.
     * @param asum [out] asum[i] \approx \int_{r[0]}^{r[i]} f(r) dr; asum[0]=0.
     * Odd-index endpoints append one trapezoid to the preceding Simpson prefix.
     * @author Peize Lin
     * @date 2016-02-14
     */
    static void Simpson_Integral_0toall 
    (
        const int mesh,
        const double * const func,
        const double * const rab,
        double * const asum
    );

    /**
     * @brief Tail integrals computed as total minus prefix.
     *
     * @param mesh [in] number of grid points, mesh >= 1 (odd or even)
     * @param func [in] function to be integrated 
     * @param rab [in] rab[i]=(dr/di)*Delta_i at each grid point. The index spacing Delta_i is 1.
     * @param asum [out] asum[i] \approx \int_{r[i]}^{r[mesh-1]} f(r) dr; no infinite-tail correction.
     * @author Peize Lin
     * @date 2016-02-14
     */
    static void Simpson_Integral_alltoinf 
    (
        const int mesh,
        const double * const func,
        const double * const rab,
        double * const asum
    );     

    //! Numerical integration on an evenly-spaced grid using Simpson's rule
    static double simpson(const int n,           //!< number of grid points
                          const double* const f, //!< function values at grid points
                          const double dx        //!< grid spacing
    );

    //! Numerical integration on an irregularly-spaced grid using Simpson's rule
    static double simpson(const int n,           //!< number of grid points
                          const double* const f, //!< function values at grid points
                          const double* const h  //!< grid spacing of length n-1, must be positive
    );

    // Grid points and weights used to generate Gauss_Legendre integrals. Returns x in [-1, 1], the integration area is -1, 1.
    // https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
    // \int_{-1}^{1} f(x) dx \approx \sum_{i = 1}^{n} w_{i} f(x_{i})
    // n is the positive number of sample points used.
    // w_{i} are quadrature weights
    // x_{i} are the roots of the n-th Legendre polynomial.
    static void Gauss_Legendre_grid_and_weight(const int n, double *x, double *weights);

    // Gauss-Legendre nodes between xmin and xmax, with oriented integration weights.
    // \int_{xmin}^{xmax} f(x) dx \approx \sum_{i = 1}^{n} w_{i} f(x_{i})
    // xl = (xmax - xmin) / 2
    // xmean = (xmax + xmin) / 2
    // \int_{xmin}^{xmax} f(x) dx = xl \int_{-1}^{1} f(xl * t + xmean) dt
    static void Gauss_Legendre_grid_and_weight(const double xmin, const double xmax, const int n, double *x, double *weights);

    // Grid points and weights used to generate Gauss_Legendre integrals. can be generated with function Lebedev_laikov_grid in math_lebedev_laikov.cpp
    static const double Lebedev_Laikov_grid110_x[110];
    static const double Lebedev_Laikov_grid110_y[110];
    static const double Lebedev_Laikov_grid110_z[110];
    static const double Lebedev_Laikov_grid110_w[110];
};

}
#endif
