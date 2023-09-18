#ifndef SPHERICAL_BESSEL_TRANSFORMER_H_
#define SPHERICAL_BESSEL_TRANSFORMER_H_

#include <memory>
#include <fftw3.h>

namespace ModuleBase
{

/**
 * @brief A class to perform spherical Bessel transforms.
 *
 * The spherical Bessel transform of a function F(x) is defined as
 *
 *                               / +inf     2
 *          G(y) = sqrt(2/pi) *  |      dx x  F(x) j (y*x)
 *                               /  0               l
 *
 * where
 *
 *          j 
 *           l
 *
 * is the l-th order spherical Bessel function of the first kind.
 *
 * This class interprets the input array as
 *
 *                   p
 *          in[i] = x [i] F(x[i])   (p being an input argument)
 *
 * and, on finish, fills the output array with
 *
 *          out[j] = G(y[j])   or  out[j] = dG(y[j])/dy
 *
 * Usage1:
 *
 *      // FFT-based method
 *
 *      SphericalBesselTransformer sbt;
 *
 *      // Default FFTW planner flag is FFTW_ESTIMATE
 *      // which is suitable for handling tasks of different sizes.
 *      sbt.radrfft(0, 2000, ...);
 *      sbt.radrfft(1, 3000, ...);
 *      sbt.radrfft(2, 4000, ...);
 *
 *      // The following flag leads to optimized FFT algorithms at the cost of
 *      // introducing large overhead during planning the FFTs.
 *      sbt.set_fftw_plan_flag(FFTW_MEASURE)
 *
 *      // FFTW plan is created at the first run
 *      // and reused for consecutive same-sized transforms.
 *      sbt.radrfft(0, 5000, ...);
 *      sbt.radrfft(1, 5000, ...);
 *      sbt.radrfft(2, 5000, ...);
 *
 *
 * Usage2:
 *
 *      // Numerical integration (Simpson's rule)
 *
 *      ModuleBase::SphericalBesselTransformer::direct(
 *          l,         // order of transform
 *          ngrid_in,  // number of input grid points
 *          grid_in,   // input grid
 *          value_in,  // input values
 *          ngrid_out, // number of output grid points
 *          grid_out,  // output grid
 *          value_out  // transformed values on the output grid
 *      );
 */
class SphericalBesselTransformer
{
  public:

    ~SphericalBesselTransformer();
    SphericalBesselTransformer(SphericalBesselTransformer const&) = delete;
    SphericalBesselTransformer& operator=(SphericalBesselTransformer const&) = delete;

    /**
     * @brief Creates a SphericalBesselTransformer object handled by a shared pointer
     *
     * This is the only way to create SphericalBesselTransformer objects.
     */
    static std::shared_ptr<SphericalBesselTransformer> create()
    {
        return std::shared_ptr<SphericalBesselTransformer>(new SphericalBesselTransformer);
    }

    /**
     * @brief Performs an l-th order spherical Bessel transform via real-input fast Fourier transforms.
     *
     * This function computes the spherical Bessel transform F(x) -> G(y) (or G's derivative)
     * with input values
     *
     *                   p
     *          in[i] = x [i] F(x[i])
     *
     * where p <= 2 is an integer, and
     *
     *                     cutoff
     *          x[i] = i * -------          i = 0, 1, 2,..., ngrid-1.
     *                     ngrid-1
     *
     * On finish, out[j] = G(y[j]) or dG(y[j])/dy where
     *
     *                      pi
     *          y[j] = j * ------           j = 0, 1, 2,..., ngrid-1.
     *                     cutoff
     *
     * @param[in]   l       order off the transform
     * @param[in]   ngrid   size of the input array
     * @param[in]   cutoff  cutoff distance of input grid
     * @param[in]   in      input values
     * @param[out]  out     transformed values
     * @param[in]   p       exponent of the extra power term in input values, must not exceed 2
     * @param[in]   deriv   if true, the derivative of the transform is computed
     *
     * @note    This function does not allocate memory for output; it must be pre-allocated.
     * @note    F(x) is supposed to be exactly zero at and after cutoff. Results would make
     *          no sense if the input is truncated at a place where F(x) is still significantly
     *          non-zero.
     * @note    FFT-based algorithm is not accurate for high l at small y. Numerical
     *          integration is automatically invoked to handle this case.
     */
    void radrfft(const int l,
                 const int ngrid,
                 const double cutoff,
                 const double* const in,
                 double* const out,
                 const int p = 0,
                 const bool deriv = false
    );

    /**
     * @brief Performs an l-th order spherical Bessel transform via numerical integration with Simpson's rule.
     *
     * This function computes the spherical Bessel transform F(x) -> G(y) (or G's derivative)
     * with input values
     *
     *                   p
     *          in[i] = x [i] F(x[i])
     *
     * where p <= 2 is an integer. On finish, out[j] = G(y[j]) or dG(y[j])/dy. 
     * x & y are specified by grid_in & grid_out, respectively.
     *
     * @param[in]   l           order of the transform
     * @param[in]   ngrid_in    size of the input array
     * @param[in]   grid_in     input grid
     * @param[in]   in          input values
     * @param[in]   ngrid_out   size of the output array
     * @param[in]   grid_out    output grid
     * @param[out]  out         transformed values on the output grid
     * @param[in]   p           exponent of the extra power term in input values, must not exceed 2
     * @param[in]   deriv       if true, the derivative of the transform is computed
     *
     * @note    This function does not allocate memory for output; it must be pre-allocated.
     * @note    Even if the input grid forms a good sampling of F(x), results would still be
     *          inaccurate for very large y values (y*dx ~ pi) because the oscillation of
     *          j_l(y*x) in this case is poorly sampled, in which case Simpson's 1/3 rule
     *          could be a bad approximation.
     * @note    p is restricted to p <= 2 in order to avoid the situation that one has to
     *          determine x^2*F(x) at x = 0 from x[i]^p*F(x[i]).
     */
    void direct(const int l,
                const int ngrid_in,
                const double* const grid_in,
                const double* const in,
                const int ngrid_out,
                const double* const grid_out,
                double* const out,
                const int p = 0,
                const bool deriv = false
    );

    /**
     * @brief Sets the FFTW planner flag.
     *
     * Recommended flags include FFTW_MEASURE and FFTW_ESTIMATE.
     *
     * FFTW_MEASURE yields optimized FFT algorithm at the cost of large overhead,
     * which is suitable for performing many consecutive transforms of the same size.
     *
     * FFTW_ESTIMATE yields less optimized FFT algorithm with much less overhead.
     *
     * @param[in]   new_flag    FFTW planner flag, usually FFTW_MEASURE or FFTW_ESTIMATE
     *
     * @note    Saved fftw_plan will be immediately destroyed if it was created with
     *          a flag other than new_flag.
     */
    void set_fftw_plan_flag(const unsigned new_flag);

    /// clear cached FFTW plan
    void fft_clear();

  private:
    SphericalBesselTransformer() = default;

    /// core function for FFT-based transform
    void radrfft_base(const int l,
                       const int ngrid,
                       const double cutoff,
                       const double* const in,
                       double* const out,
                       const int p = 0
    );

    /// Internal buffer used for in-place real-input FFT (interpreted as double* on input)
    fftw_complex* f_ = nullptr;

    /// FFTW plan saved for reuse
    fftw_plan rfft_plan_ = nullptr;

    /// Size of the planned FFT
    int sz_planned_ = -1;

    /// Planner flag used to create rfft_plan_
    unsigned fftw_plan_flag_ = FFTW_ESTIMATE;

    /// Applies an in-place 1-d real-input discrete Fourier transform to the internal buffer
    void rfft_in_place();

    /// Buffer allocation and plan creation for a real-input FFT of size N.
    void rfft_prepare(const int N);

    /**
     * @brief Polynomial coefficients in the sin & cos expression of the spherical Bessel function.
     *
     * The l-th order spherical Bessel function of the first kind can be expressed as
     *
     *                          sin(x)*P(x) + cos(x)*Q(x)
     *                  j (x) = -------------------------
     *                   l               l+1
     *                                  x
     *
     * where P(x) and Q(x) are polynomials of degree no more than l. This function
     * returns the coefficients within those polynomials.
     *
     * @param[in]   get_sine    if true, compute the coefficients of polynomials attached to sin
     * @param[in]   l           order of the spherical Bessel function
     * @param[in]   n           degree of the polynomial term whose coefficient is computed
     *
     * @return  The polynomial coefficient of the n-th power term in the sin & cos expression
     *          of the l-th order spherical Bessel functions of the first kind.
     *
     * @note    Coefficients grow very quickly as l increases. Currently l is capped at 17
     *          since some coefficients exceed 2^63-1 for l >= 18.
     *
     */
    long long int spherical_bessel_sincos_polycoef(
        const bool get_sine,
        const int l,
        const int n
    );

    /// Computes & stores the values of spherical Bessel function on the given transform grid
    void cache(const int l,
               const int ngrid_in,
               const double* const grid_in,
               const int ngrid_out,
               const double* const grid_out,
               const bool deriv);

    /**
     * @name Cached function values for direct integration
     *                                                                                  */
    ///@{
    bool is_deriv_ = false; ///< if true, the cached values are derivatives of the spherical Bessel function
    int l_ = -1;            ///< order of the cached spherical Bessel function

    int ngrid_in_ = 0;
    int ngrid_out_ = 0;
    double* grid_in_ = nullptr;
    double* grid_out_ = nullptr;

    /// jl_[i*ngrid_in_ + j] = f(l, grid_out_[i] * grid_in_[j]) where f is sphbesj or dsphbesj
    double* jl_ = nullptr;
    ///@}

}; // class SphericalBesselTransformer

} // namespace ModuleBase

#endif
