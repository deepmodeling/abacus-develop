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
 *  Currently the supported grids must satisfy
 *
 *          x[i] = i*dx             i = 0, 1, 2, ..., N
 *          y[j] = j*pi/(N*dx)      j = 0, 1, 2, ..., N
 *
 *  That is, the input grid must be uniform and starts from 0, and there is no
 *  freedom to choose the output grid once the input grid is specified.
 *
 *
 *  Usage:
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
     *  where p is an arbitrary integer, and
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
     *                                                                                      */
    void radrfft(const int l,            //!< [in] order of the transform
                 const int ngrid,        //!< [in] size of the input array
                 const double cutoff,    //!< [in] cutoff distance of input values
                 const double* const in, //!< [in] input values
                 double* const out,      //!< [out] transformed values
                 const int p = 0         //!< [in] exponent of the pre-multiplied power term in input values
    );

    //! Performs an l-th order spherical Bessel transform via direct integration with Simpson's rule.
    /*!
     *  To be implemented...
     *                                                                                      */
    void direct(const int l,               //!< [in] order of the transform
                const int nin,             //!< [in] size of the input array
                const double* const igrid, //!< [in] input grid
                const double* const in,    //!< [in] values on the input grid
                const int nout,            //!< [in] size of the output array
                const double* const ogrid, //!< [in] output grid
                double* const out,         //!< [out] transformed values
                const int p = 0            //!< [in] exponent of the pre-multiplied power term in input values
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
    int spherical_bessel_sincos_polycoef(
        const bool get_sine, //!< [in] specifies if the returned coefficient is associated with sine
        const int l,         //!< [in] order of the spherical Bessel function
        const int n          //!< [in] degree of the polynomial term whose coefficient is computed
    );

}; // class SphericalBesselTransformer

} // namespace ModuleBase

#endif
