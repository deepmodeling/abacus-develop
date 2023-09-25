#ifndef SPHERICAL_BESSEL_TRANSFORMER_H_
#define SPHERICAL_BESSEL_TRANSFORMER_H_

#include <cstddef>
#include <memory>
#include <fftw3.h>

namespace ModuleBase
{

/**
 * @brief A class to perform spherical Bessel transforms.
 *
 * @note This class is implemented as an opaque shared pointer.
 *       Copy-constructed or assigned objects share the same
 *       underlying implementation class object.
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
 *
 * Usage:
 *
 *      SphericalBesselTransformer sbt;
 *
 *      // Usage 1: FFT-based method
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
 *      // Usage 2: Numerical integration (Simpson's rule)
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
 */
class SphericalBesselTransformer
{
public:
    SphericalBesselTransformer();
    SphericalBesselTransformer(std::nullptr_t);
    SphericalBesselTransformer(SphericalBesselTransformer const&) = default;
    SphericalBesselTransformer& operator=(const SphericalBesselTransformer&) = default;
    ~SphericalBesselTransformer() = default;

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
    ) const;

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
    ) const;

    /**
     * @brief Sets the FFTW planner flag.
     *
     * Accepted flags include FFTW_MEASURE and FFTW_ESTIMATE:
     *
     * - FFTW_MEASURE  yields optimized FFT algorithm at the cost of large overhead;
     * - FFTW_ESTIMATE yields less optimized FFT algorithm with much less overhead.
     *
     * @param[in]   new_flag    FFTW planner flag, FFTW_MEASURE or FFTW_ESTIMATE
     *
     * @note    Saved fftw_plan will be immediately destroyed if it was created with
     *          a different flag from new_flag.
     */
    void set_fftw_plan_flag(const unsigned new_flag) const;

    /// Clears cached FFTW plan.
    void fft_clear() const;

    /// Initializes impl_ if it is nullptr; does nothing otherwise.
    void init();

    /// Returns whether impl_ is ready to use.
    bool is_ready() const { return impl_.get(); }

    /// Returns whether two objects share the same underlying Impl object.
    inline bool operator==(SphericalBesselTransformer const& rhs) const { return this->impl_ == rhs.impl_; }

private:
    class Impl; // forward declaration for detailed implementation class
    std::shared_ptr<Impl> impl_;
};

} // namespace ModuleBase

#endif
