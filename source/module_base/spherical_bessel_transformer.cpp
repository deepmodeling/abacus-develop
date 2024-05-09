#include "module_base/spherical_bessel_transformer.h"

#include <vector>
#include <fftw3.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <functional>
#include <numeric>

#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"

namespace ModuleBase
{

//======================================================================
//
//                      Implementation Class
//
//======================================================================
class SphericalBesselTransformer::Impl
{

public:

    explicit Impl(bool cache_enabled): cache_enabled_(cache_enabled) {}
    ~Impl() { _rfft_clear(); };

    Impl(Impl const&) = delete;
    Impl(Impl&&) = delete;

    Impl& operator=(Impl const&) = delete;
    Impl& operator=(Impl&&) = delete;


    void radrfft(
        int l,
        int ngrid,
        double cutoff,
        const double* in,
        double* out,
        int p = 0
    );


    void direct(
        int l,
        int ngrid_in,
        const double* grid_in,
        const double* in,
        int ngrid_out,
        const double* grid_out,
        double* out,
        int p = 0
    );


private:

    /// core function for FFT-based transform
    void _radrfft_base(
        int l,
        int ngrid,
        double cutoff,
        const double* in,
        double* out,
        int p = 0
    );

    // NOTE: FFTW suggests using fftw_malloc and fftw_free for memory (de)allocation
    // in order to enforce some strong requirements for memory alignment.

    /// internal buffer used for in-place real-input FFT
    fftw_complex* f_ = nullptr;

    /// FFTW plan
    fftw_plan rfft_plan_ = nullptr;

    /// size of the planned FFT
    int sz_planned_ = -1;

    /// buffer allocation and plan creation for a real-input FFT of size n.
    void _rfft_prepare(int n);

    /// clear the FFTW plan and buffer
    void _rfft_clear();


    /**
     * @brief Polynomial coefficients in the sin & cos expression of spherical Bessel function.
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
     * @param[in]   of_sine     if true, compute the coefficients of polynomials attached to sin
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
    long long int _polycoef(bool of_sine, int l, int n);


    /// tabulate spherical Bessel function on the given transform grid
    void _tabulate(
        int l,
        size_t ngrid_in,
        const double* grid_in,
        size_t ngrid_out,
        const double* grid_out
    );

    /// if true, tabulated j_l(grid_out[j] * grid_in[i]) will be cached
    bool cache_enabled_ = false;

    /// order of the cached spherical Bessel function
    int l_ = -1;

    /// cached input grid
    std::vector<double> grid_in_;

    /// cached output grid
    std::vector<double> grid_out_;

    /// cached spherical Bessel function values
    std::vector<double> jl_;
    // jl_[j*ngrid_in_ + i] = sphbesj(l, grid_out_[j] * grid_in_[i])

}; // class SphericalBesselTransformer::Impl


long long int SphericalBesselTransformer::Impl::_polycoef(bool of_sine, int l, int n)
{
    /*
     * The sin & cos coefficients follow the same recurrence relation
     *
     *        c(l,n) = (2l-1) * c(l-1,n) - c(l-2,n-2)
     *
     * with different initial conditions:
     *
     *          cos             sin
     *
     *      c(0,0) =  0     c(0,0) = 1
     *      c(1,0) =  0     c(1,0) = 1
     *      c(1,1) = -1     c(1,1) = 0
     *                                                              */
    assert(l >= 0 && n >= 0);
    assert(l <= 17 && "Impl::_polycoef: some coefficients exceed LLONG_MAX (2^63-1) for l >= 18");

    if (of_sine)
    {
        if (n % 2 == 1 || n > l)
        {
            return 0;
        }
        if (l == 0 && n == 0)
        {
            return 1;
        }
        if (l == 1 && n == 0)
        {
            return 1;
        }
    }
    else
    {
        if (n % 2 == 0 || n > l)
        {
            return 0;
        }
        if (l == 1 && n == 1)
        {
            return -1;
        }
    }

    return (2 * l - 1) * _polycoef(of_sine, l - 1, n)
        - (n >= 2 ? _polycoef(of_sine, l - 2, n - 2) : 0);
}


void SphericalBesselTransformer::Impl::_radrfft_base(
    int l,
    int ngrid,
    double cutoff,
    const double* in,
    double* out,
    int p
)
{
    // this function does not support in-place transform (but radrfft does)
    assert(in != out);

    double pi = std::acos(-1.0);
    int n = ngrid - 1;
    double dx = cutoff / n;
    double dy = pi / cutoff;
    double pref = dx / std::sqrt(2. * pi);

    std::fill(out, out + ngrid, 0.0);

    // rfft buffer (f_) allocation and plan creation
    _rfft_prepare(2 * n);

    // introduced to help manipulate f_ ( double(*)[2] ) as an array of double
    double* f = &f_[0][0];

    for (int m = 0; m <= l; ++m)
    {
        // m even --> sin; f[2*n-i] = -f[i]; out += -imag(rfft(f)) / y^(l+1-m)
        // m odd  --> cos; f[2*n-i] = +f[i]; out += +real(rfft(f)) / y^(l+1-m)
        bool flag = (m % 2 == 0);

        long long int sincos_coef = _polycoef(flag, l, m);

        f[0] = f[n] = 0.0;
        for (int i = 1; i != n; ++i)
        {
            f[i] = pref * sincos_coef * in[i] * std::pow(i * dx, m + 1 - l - p);
            f[2 * n - i] = flag ? -f[i] : f[i];
        }

        fftw_execute(rfft_plan_); // perform in-place rfft on f_

        // summing up the series by ( ... ( ( g0/y + g1 )/y + g2 )/y + ... + gl )/y
        // out[0] is handled independently
        for (int j = 1; j <= n; ++j)
        {
            out[j] = (out[j] + (flag ? -f_[j][1] : f_[j][0])) / (j * dy);
        }
    }

    // out[0] is done by direct integration
    // note that only the zeroth order spherical Bessel function is nonzero at 0
    if (l == 0)
    {
        for (int i = 0; i <= n; ++i)
        {
            out[0] += 2. * pref * in[i] * std::pow(i * dx, 2 - p); // p <= 2 is required!
        }
    }
}

void SphericalBesselTransformer::Impl::radrfft(
    int l,
    int ngrid,
    double cutoff,
    const double* in,
    double* out,
    int p
)
{
    /*
     * An l-th order spherical Bessel transform F(x) -> G(y) can be expressed in terms of Fourier transforms:
     *
     *           l
     *          ---    1        / +inf            -iyx
     * G(y) =   \   -------  Re |      dr f(n,x) e
     *          /     l+1-n     / -inf
     *          ---  y
     *          n=0
     *
     * where
     *
     *              1                              F(x)
     * f(n,x) = ---------- [ c(l,n) + i*s(l,n) ] --------
     *          sqrt(2*pi)                         l-n-1   ,
     *                                            x
     *
     * c(l,n) / s(l,n) are sin / cos coefficients from polycoef, and
     * the domain of F(x) is extended to negative values by defining F(-x) = pow(-1,l)*F(x).
     *
     * With an appropriate grid, the Fourier transform can be approximated by a discrete Fourier transform.
     * Note that given l & n, c(l,n) and s(l,n) cannot be both non-zero. Therefore, each FFT input
     * array is either purely real or purely imaginary, which suggests the use of real-input FFT.
     *
     * If deriv == true, dG(y)/dy is calculated instead of G(y). This is done by using the recurrence
     * relation of j_l(x) to convert the original transform to two Spherical Bessel transforms of order
     * l-1 and l+1 respectively.
     *                                                                                                  */
    assert(l >= 0);
    assert(ngrid > 1);
    assert(p <= 2);

    std::vector<double> out_tmp(ngrid);
    _radrfft_base(l, ngrid, cutoff, in, out_tmp.data(), p);

    // FFT-based method does not yield accurate results for small y at high l
    // use numerical integration in this case
    int n_direct = (l == 0) ? 0 : static_cast<int>(ngrid * std::pow(1e-8, 1. / l));
    if (n_direct > 0)
    {
        std::vector<double> buffer(ngrid + n_direct);
        double* grid_in = buffer.data();
        double* grid_out = grid_in + ngrid;

        double dx = cutoff / (ngrid - 1);
        double dy = std::acos(-1.0) / cutoff;

        std::for_each(grid_in, grid_in + ngrid,
            [&](double& x) { x = (&x - grid_in) * dx; });
        std::for_each(grid_out, grid_out + n_direct,
            [&](double& y) { y = ((&y - grid_out) + 1) * dy; });

        direct(l, ngrid, grid_in, in, n_direct, grid_out, &out_tmp[1], p);
    }

    std::copy(out_tmp.begin(), out_tmp.end(), out);
}


void SphericalBesselTransformer::Impl::direct(
    int l,
    int ngrid_in,
    const double* grid_in,
    const double* in,
    int ngrid_out,
    const double* grid_out,
    double* out,
    int p
)
{
    assert(p <= 2);
    assert(grid_in[0] >= 0.0 && grid_out[0] >= 0.0);
    assert(std::is_sorted(grid_in, grid_in + ngrid_in, std::less_equal<double>()));
    assert(std::is_sorted(grid_out, grid_out + ngrid_out, std::less_equal<double>()));

    std::vector<double> buffer(3 * ngrid_in);
    double* rab = buffer.data();
    double* tmp = rab + ngrid_in; // integrand without the jl part
    double* integrand = tmp + ngrid_in; // integrand

    std::adjacent_difference(grid_in, grid_in + ngrid_in, rab);

    std::copy(in, in + ngrid_in, tmp);
    std::for_each(tmp, tmp + ngrid_in,
        [&](double& x) { x *= std::pow(grid_in[&x - tmp], 2 - p); });

    // compute spherical Bessel function on the grid and store the results in jl_
    // (will be cleared at the end of this function if cache is disabled)
    _tabulate(l, ngrid_in, grid_in, ngrid_out, grid_out);

    for (int j = 0; j < ngrid_out; ++j)
    {
        std::transform(tmp, tmp + ngrid_in, &jl_[j * grid_in_.size()], integrand,
                       std::multiplies<double>());
        out[j] = ModuleBase::Integral::simpson(ngrid_in, integrand, &rab[1]);
    }

    double pref = std::sqrt(2.0 / std::acos(-1.0));
    std::for_each(out, out + ngrid_out, [pref](double& x) { x *= pref; });

    if (!cache_enabled_)
    {
        std::vector<double>().swap(grid_in_);
        std::vector<double>().swap(grid_out_);
        std::vector<double>().swap(jl_);
    }
}


void SphericalBesselTransformer::Impl::_rfft_prepare(int n)
{
    if (n != sz_planned_)
    {
        fftw_free(f_);
        f_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (n / 2 + 1));
        sz_planned_ = n;
        fftw_destroy_plan(rfft_plan_);
        rfft_plan_ = fftw_plan_dft_r2c_1d(n, &f_[0][0], f_, FFTW_ESTIMATE);
    }
}


void SphericalBesselTransformer::Impl::_tabulate(
    int l,
    size_t ngrid_in,
    const double* grid_in,
    size_t ngrid_out,
    const double* grid_out
)
{
    bool is_cached = cache_enabled_ && l == l_
                     && ngrid_in <= grid_in_.size() && ngrid_out <= grid_out_.size() 
                     && std::equal(grid_in, grid_in + ngrid_in, grid_in_.begin())
                     && std::equal(grid_out, grid_out + ngrid_out, grid_out_.begin());
    if (is_cached)
    {
        return;
    }

    l_ = l;
    grid_in_ = std::vector<double>(grid_in, grid_in + ngrid_in);
    grid_out_ = std::vector<double>(grid_out, grid_out + ngrid_out);
    jl_.resize(grid_in_.size() * grid_out_.size());
    for (size_t j = 0; j != ngrid_out; ++j)
    {
        ModuleBase::Sphbes::sphbesj(ngrid_in, grid_in, grid_out[j], l, &jl_[j * ngrid_in]);
    }
}


void SphericalBesselTransformer::Impl::_rfft_clear()
{
    if (rfft_plan_)
    {
        fftw_destroy_plan(rfft_plan_);
        rfft_plan_ = nullptr;
    }

    if (f_)
    {
        fftw_free(f_);
        f_ = nullptr;
    }

    sz_planned_ = -1;
}

//======================================================================
//
//                      Interface Class
//
//======================================================================
SphericalBesselTransformer::SphericalBesselTransformer(bool cache_enabled)
    : impl_(new Impl(cache_enabled))
{}


void SphericalBesselTransformer::radrfft(
    int l,
    int ngrid,
    double cutoff,
    const double* in,
    double* out,
    int p
) const
{
    impl_->radrfft(l, ngrid, cutoff, in, out, p);
}


void SphericalBesselTransformer::direct(
    int l,
    int ngrid_in,
    const double* grid_in,
    const double* in,
    int ngrid_out,
    const double* grid_out,
    double* out,
    int p
) const
{
    impl_->direct(l, ngrid_in, grid_in, in, ngrid_out, grid_out, out, p);
}


} // namespace ModuleBase
