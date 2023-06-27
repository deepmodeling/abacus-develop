#ifndef SPHERICAL_BESSEL_TRANSFORMER_H_
#define SPHERICAL_BESSEL_TRANSFORMER_H_

#include <fftw3.h>

namespace ModuleBase
{

//! A class to perform spherical Bessel transforms.
/*!
 *  The spherical Bessel transform of a function F(x) is defined as
 *
 *                               / +inf     2
 *          G(y) = sqrt(2/pi) *  |      dx x  F(x) j (x)
 *                               /  0               l
 *
 *  where
 *
 *          j (x)
 *           l
 *
 *  is the l-th order spherical Bessel function of the first kind.
 *
 *  This class interprets the input array as
 *
 *                   p
 *          in[i] = x [i] F(x[i])   (p being an input argument)
 *
 *  and, on finish, fills the output array with
 *
 *          out[j] = G(y[j])
 *
 *  Usage1:
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
 *  Usage2:
 *
 *      // Quadrature (Simpson's rule)
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
 *
 *                                                                              */
class SphericalBesselTransformer
{

  public:
    SphericalBesselTransformer(){};
    ~SphericalBesselTransformer();

    SphericalBesselTransformer(SphericalBesselTransformer const&) = delete;
    SphericalBesselTransformer& operator=(SphericalBesselTransformer const&) = delete;

    //! Performs an l-th order spherical Bessel transform via real-input fast Fourier transforms.
    /*!
     *  This function computes the spherical Bessel transform F(x) -> G(y) with input values
     *
     *                   p
     *          in[i] = x [i] F(x[i])
     *
     *  where p <= 2 is an integer, and
     *
     *                     cutoff
     *          x[i] = i * -------          i = 0, 1, 2,..., ngrid-1.
     *                     ngrid-1
     *
     *  On finish, out[j] = G(y[j]) where
     *
     *                      pi
     *          y[j] = j * ------           j = 0, 1, 2,..., ngrid-1.
     *                     cutoff
     *
     *
     *  @note   This function does not allocate memory for output; it must be pre-allocated.
     *  @note   F(x) is supposed to be exactly zero at and after cutoff. Results would make
     *          no sense if the input is truncated at a place where F(x) is still significantly
     *          non-zero.
     *  @note   FFT-based algorithm is not accurate for high l at small y. Numerical
     *          integration in used to handle this case.
     *  @note   p must not exceed 2 as required by numerical integration, @see @ref direct
     *                                                                                      */
    void radrfft(const int l,            //!< [in] order of the transform
                 const int ngrid,        //!< [in] size of the input array
                 const double cutoff,    //!< [in] cutoff distance of input values
                 const double* const in, //!< [in] input values
                 double* const out,      //!< [out] transformed values
                 const int p = 0         //!< [in] exponent of the pre-multiplied power term in input values,
                                         //!<      must not exceed 2
    );

    //! Performs an l-th order spherical Bessel transform via numerical integration with Simpson's rule.
    /*!
     *  This function computes the spherical Bessel transform F(x) -> G(y) with input values
     *
     *                   p
     *          in[i] = x [i] F(x[i])
     *
     *  where p <= 2 is an integer. On finish, out[j] = G(y[j]). x & y are specified by
     *  grid_in & grid_out, respectively.
     *
     *  @note   This function does not allocate memory for output; it must be pre-allocated.
     *  @note   Even if the input grid forms a good sampling of F(x), results would still be
     *          inaccurate for very large y values (y*dx ~ pi) because the oscillation of
     *          j_l(y*x) in this case is poorly sampled, in which case Simpson's 1/3 rule
     *          could be a bad approximation.
     *  @note   p is limited to p <= 2 in order to avoid the problem of determining
     *          x^2*F(x) at x = 0 from x[i]^p*F(x[i])
     *                                                                                      */
    static void direct(const int l,                  //!< [in] order of the transform
                       const int ngrid_in,           //!< [in] size of the input array
                       const double* const grid_in,  //!< [in] input grid
                       const double* const in,       //!< [in] input values
                       const int ngrid_out,          //!< [in] size of the output array
                       const double* const grid_out, //!< [in] output grid
                       double* const out,            //!< [out] transformed values
                       const int p = 0               //!< [in] exponent of the pre-multiplied power
                                                     //!<      term in input values, must not exceed 2
    );

    //! Sets the FFTW planner flag.
    /*!
     *  Recommended flags include FFTW_MEASURE and FFTW_ESTIMATE.
     *
     *  FFTW_MEASURE leads to optimized FFT algorithm at the cost of large overhead,
     *  which is suitable for performing many consecutive transforms of the same size.
     *
     *  FFTW_ESTIMATE leads to less optimized FFT algorithm with much less overhead,
     *  which is suitable for isolated tasks.
     *
     *  @note   Saved fftw_plan will be immediately destroyed if it was created with
     *          a flag other than new_flag.
     *                                                                                      */
    void set_fftw_plan_flag(const unsigned new_flag /*!< [in] FFTW planner flag */);

  private:
    //! Internal buffer used for in-place real-input FFT (interpreted as double* on input)
    fftw_complex* f_ = nullptr;

    //! FFTW plan saved for reuse
    fftw_plan rfft_plan_ = nullptr;

    //! Size of the planned FFT
    int sz_planned_ = -1;

    //! Planner flag used to create rfft_plan_
    unsigned fftw_plan_flag_ = FFTW_ESTIMATE;

    //! Applies an in-place 1-d real-input discrete Fourier transform to the internal buffer
    void rfft_in_place();

    //! Buffer allocation and plan creation for a real-input FFT
    void rfft_prepare(const int N /*!< [in] size of the FFT to plan ahead */);

    //! Polynomial coefficients in the sin & cos expression of the spherical Bessel function.
    /*!
     *  The l-th order spherical Bessel function of the first kind can be expressed as
     *
     *                          sin(x)*P(x) + cos(x)*Q(x)
     *                  j (x) = -------------------------
     *                   l               l+1
     *                                  x
     *
     *  where P(x) and Q(x) are polynomials of degree no more than l. This function
     *  returns the coefficients within those polynomials.
     *
     *
     *  @note   Coefficients grow very quickly as l increases. Currently l is capped at 10
     *          since some coefficients exceed INT_MAX (2^31-1) for l >= 11.
     *
     *  @return The polynomial coefficient of the n-th power term in the sin & cos expression
     *          of the l-th order spherical Bessel functions of the first kind.
     *                                                                                      */
    long long int spherical_bessel_sincos_polycoef(
        const bool get_sine, //!< [in] specifies if the returned coefficient is associated with sine
        const int l,         //!< [in] order of the spherical Bessel function
        const int n          //!< [in] degree of the polynomial term whose coefficient is computed
    );

}; // class SphericalBesselTransformer

} // namespace ModuleBase

#endif
