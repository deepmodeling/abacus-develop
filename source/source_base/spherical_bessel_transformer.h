#ifndef SPHERICAL_BESSEL_TRANSFORMER_H_
#define SPHERICAL_BESSEL_TRANSFORMER_H_

#include <memory>
#include <fftw3.h>

namespace ModuleBase
{

/**
 * @brief A class that provides spherical Bessel transforms.
 *
 * @note This class is implemented as an opaque shared pointer. The underlying
 * object of the implementation class, which may cache some tabulated function
 * values, is shared via copy-construction and copy-assignment.
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
 *              out[j] = G(y[j])
 *
 *
 * Usage:
 *
 *      // cache is disabled by default
 *      SphericalBesselTransformer sbt;
 *
 *      // cache enabled mode
 *      // faster for multiple transforms with the same size & order
 *      SphericalBesselTransformer sbt2(true);
 *
 *      //---------------------------------------------------------------------
 *      //                      FFT-based algorithm
 *      //---------------------------------------------------------------------
 *      // basic usage
 *      sbt.radrfft(l, ngrid, cutoff, in, out);
 *
 *      // perform SBT on F(r) with input values r*F(r)
 *      sbt.radrfft(l, ngrid, cutoff, in, out, 1);
 *
 *      //---------------------------------------------------------------------
 *      //              numerical integration (Simpson's rule)
 *      //---------------------------------------------------------------------
 *      // basic usage
 *      sbt.direct(l, ngrid_in, grid_in, value_in, ngrid_out, grid_out, value_out);
 *
 *      // perform SBT on F(r) with input values r^2*F(r)
 *      sbt.direct(l, ngrid_in, grid_in, value_in, ngrid_out, grid_out, value_out, 2);
 *
 */
class SphericalBesselTransformer
{
public:
    SphericalBesselTransformer(const bool cache_enabled = false);
    ~SphericalBesselTransformer() = default;

    SphericalBesselTransformer(SphericalBesselTransformer const&) = default;
    SphericalBesselTransformer(SphericalBesselTransformer &&) = default;

    SphericalBesselTransformer& operator=(const SphericalBesselTransformer&) = default;
    SphericalBesselTransformer& operator=(SphericalBesselTransformer&&) = default;

    /**
     * @brief Spherical Bessel transform via fast Fourier transforms.
     *
     * This function computes the spherical Bessel transform F(x) -> G(y) with input
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
     * On finish, out[j] = G(y[j]) where
     *
     *                      pi
     *          y[j] = j * ------           j = 0, 1, 2,..., ngrid-1.
     *                     cutoff
     *
     * @param[in]   l           order of the transform
     * @param[in]   ngrid       number of grid points (same for input and output)
     * @param[in]   cutoff      cutoff distance of input grid
     * @param[in]   in          input values
     * @param[out]  out         transformed values
     * @param[in]   p           exponent of the extra power term in input values
     *                          (must not exceed 2)
     *
     * @note    F(x) is supposed to be exactly zero at and after cutoff. Results would
     *          make no sense if the input is truncated at a place where F(x) is still
     *          significantly non-zero.
     * @note    FFT-based algorithm is not accurate for high l at small y. Numerical
     *          integration via Simpson's rule is implicitly invoked to handle this case.
     * @note    p is restricted to p <= 2 in order to avoid the situation that one has to
     *          determine x^2*F(x) at x = 0 from x[i]^p*F(x[i]).
     */
    void radrfft(
        const int l,
        const int ngrid,
        const double cutoff,
        const double* const in,
        double* const out,
        const int p = 0
    ) const;


    /**
     * @brief Spherical Bessel transform via numerical integration with Simpson's rule.
     *
     * This function computes the spherical Bessel transform F(x) -> G(y) with input
     *
     *                   p
     *          in[i] = x [i] F(x[i])
     *
     * where p <= 2 is an integer. On finish, out[j] = G(y[j]).
     *
     * @param[in]   l           order of the transform
     * @param[in]   ngrid_in    number of the input grid points
     * @param[in]   grid_in     input grid
     * @param[in]   in          input values
     * @param[in]   ngrid_out   number of the output grid points
     * @param[in]   grid_out    output grid
     * @param[out]  out         transformed values on the output grid
     * @param[in]   p           exponent of the extra power term in input values
     *                          (must not exceed 2)
     *
     * @note    Even if the input grid forms a good sampling of F(x), results would still be
     *          inaccurate for very large y values (y*dx ~ pi) because the oscillation of
     *          j_l(y*x) in this case is poorly sampled, in which case Simpson's 1/3 rule
     *          could be a bad approximation.
     * @note    p is restricted to p <= 2 in order to avoid the situation that one has to
     *          determine x^2*F(x) at x = 0 from x[i]^p*F(x[i]).
     *
     */
    void direct(
        const int l,
        const int ngrid_in,
        const double* const grid_in,
        const double* const in,
        const int ngrid_out,
        const double* const grid_out,
        double* const out,
        const int p = 0
    ) const;


    /// total heap usage (in bytes) from the FFTW buffer and tabulated jl
    size_t heap_usage() const;

    /// clear the FFTW plan & buffer as well as the tabulated jl
    void clear();

    /// check if two objects share the same underlying implementation object
    bool operator==(const SphericalBesselTransformer& rhs) const { return impl_ == rhs.impl_; }


private:

    class Impl; // forward declaration
    std::shared_ptr<Impl> impl_;
};

} // namespace ModuleBase

#endif
