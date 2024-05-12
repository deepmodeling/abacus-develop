#include "module_base/spherical_bessel_transformer.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cassert>
#include <complex>

#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"

namespace ModuleBase
{

//**********************************************************************
//                          Implementation
//**********************************************************************
class SphericalBesselTransformer::Impl
{

public:

    Impl(const bool cache_enabled = false);
    ~Impl() { _rfft_clear(); };

    Impl(Impl const&) = delete;
    Impl(Impl&&) = delete;

    Impl& operator=(Impl const&) = delete;
    Impl& operator=(Impl&&) = delete;

    // see the interface class for details
    void radrfft(
        const int l,
        const int ngrid,
        const double cutoff,
        const double* const in,
        double* const out,
        const int p = 0
    );

    // see the interface class for details
    void direct(
        const int l,
        const int ngrid_in,
        const double* const grid_in,
        const double* const in,
        const int ngrid_out,
        const double* const grid_out,
        double* const out,
        const int p = 0
    );

    // total heap usage (in bytes) from the FFTW buffer and tabulated jl
    size_t heap_usage() const;

    // clear the FFTW plan & buffer as well as the tabulated jl
    void clear();


private:

    //*****************************************************************
    //                      FFT-based algorithm
    //*****************************************************************
    // NOTE: The reason of using a raw pointer to handle the FFT buffer is because
    // FFTW suggests using its own memory allocation utility, which guarantee that
    // the returned pointer obeys any special alignment restrictions imposed by any
    // algorithm in FFTW.

    /// buffer used for in-place real-input FFT
    double* f_ = nullptr;

    /// FFTW plan
    fftw_plan rfft_plan_ = nullptr;

    /// size of the planned FFT
    int sz_planned_ = 0;

    /// planner flag used to create rfft_plan_
    const unsigned fftw_plan_flag_;

    /// buffer allocation and plan creation for in-place real-input FFT
    void _rfft_prepare(const int n);

    /// clear the FFTW plan and buffer
    void _rfft_clear();


    //*****************************************************************
    //                      numerical integration
    //*****************************************************************
    /// if true, jl values used in numerical integration will be cached
    const bool cache_enabled_ = false;

    /// order of the cached spherical Bessel function
    int l_ = -1;

    /// cached input grid
    std::vector<double> grid_in_;

    /// cached output grid
    std::vector<double> grid_out_;

    /// cached spherical Bessel function values
    std::vector<double> jl_;

    /// tabulate spherical Bessel function on the mesh grid
    void _tabulate(
        const int l,
        const int ngrid_in,
        const double* const grid_in,
        const int ngrid_out,
        const double* const grid_out
    );

    /// clear the tabulated jl
    void _table_clear();

}; // class SphericalBesselTransformer::Impl


SphericalBesselTransformer::Impl::Impl(const bool cache_enabled):
    // NOTE: For the current usage of this class, the performance gain
    // by using FFTW_MEASURE instead of FFTW_ESTIMATE usually does not
    // worth the extra overhead, and may cause a timeout of the integrated
    // test. This might need more investigation and change in the future.
    //fftw_plan_flag_(cache_enabled ? FFTW_MEASURE : FFTW_ESTIMATE),
    fftw_plan_flag_(FFTW_ESTIMATE),
    cache_enabled_(cache_enabled)
{}


void SphericalBesselTransformer::Impl::radrfft(
    const int l,
    const int ngrid,
    const double cutoff,
    const double* const in,
    double* const out,
    const int p
)
{
    /*
     * An l-th order spherical Bessel transform F(x) -> G(y) can be expressed
     * in terms of Fourier transforms:
     *
     *                l
     *               ---    1        / +inf            -iyx
     *      G(y) =   \   -------  Re |      dr f(n,x) e
     *               /     l+1-n     / -inf
     *               ---  y
     *               n=0
     *
     * where
     *
     *                   1                              F(x)
     *      f(n,x) = ---------- [ q(l,n) + i*p(l,n) ] --------
     *               sqrt(2*pi)                         l-n-1
     *                                                 x
     *
     * p,q are polynomial coefficients in Rayleigh's formula (see below),
     * and the domain of F(x) is extended to negative values by letting
     * F(-x) = pow(-1,l)*F(x).
     *
     * Note that given l & n, q(l,n) and p(l,n) cannot be both non-zero.
     * Therefore, each FFT input array is either purely real or purely
     * imaginary, which suggests the use of real-input FFT.
     *
     */
    assert(l >= 0);
    assert(ngrid > 1);
    assert(p <= 2);

    const double pi = std::acos(-1.0);
    const int n = ngrid - 1;
    const double dx = cutoff / n;
    const double dy = pi / cutoff;
    const double pref = dx / std::sqrt(2. * pi);

    // temporary storage for the output (in order to support in-place transform)
    std::vector<double> tmp(ngrid);

    // The l-th order spherical Bessel function of the first kind can be expressed as
    //
    //                          sin(x)*P(l,x) + cos(x)*Q(l,x)
    //                  j (x) = -----------------------------
    //                   l               l+1
    //                                  x
    //
    // where P(l,x) and Q(l,x) are polynomials of degree no more than l. Their polynomial
    // coefficients follow the same recurrence relation
    //
    //        c(l,m) = (2l-1) * c(l-1,m) - c(l-2,m-2)
    //
    // with different initial conditions:
    //
    //          cos             sin
    //
    //      c(0,0) =  0     c(0,0) = 1
    //      c(1,0) =  0     c(1,0) = 1
    //      c(1,1) = -1     c(1,1) = 0
    //
    std::vector<std::complex<double>> c((l+1) * (l+2) / 2); // cos->real; sin->imag
    auto idx = [](int ll, int m) { return (ll+1)*ll/2 + m; };

    c[idx(0, 0)] = {0, 1};
    if (l > 0)
    {
        c[idx(1, 0)] = {0, 1};
        c[idx(1, 1)] = {-1, 0};
    }

    for (int ll = 2; ll <= l; ++ll)
    {
        for (int m = 0; m < ll; ++m)
        {
            c[idx(ll,m)] = (2*ll-1.0) * c[idx(ll-1, m)] - (m >= 2 ? c[idx(ll-2, m-2)] : 0);
        }
        c[idx(ll,ll)] = - c[idx(ll-2, ll-2)];
    }

    _rfft_prepare(2 * n);

    bool is_imag = true;
    int sign = -1;
    for (int m = 0; m <= l; ++m)
    {
        // m even --> sin; f[2*n-i] = -f[i]; out += -imag(rfft(f)) / y^(l+1-m)
        // m odd  --> cos; f[2*n-i] = +f[i]; out += +real(rfft(f)) / y^(l+1-m)
        const double coef = reinterpret_cast<double(&)[2]>(c[idx(l,m)])[is_imag];

        f_[0] = f_[n] = 0.0;
        for (int i = 1; i != n; ++i)
        {
            f_[i] = pref * coef * in[i] * std::pow(i * dx, m + 1 - l - p);
            f_[2 * n - i] = sign * f_[i];
        }

        fftw_execute(rfft_plan_); // perform in-place rfft on f_

        // sum up the series by ( ... ( ( g0/y + g1 )/y + g2 )/y + ... + gl )/y
        // out[0] is handled later by direct integration
        for (int j = 1; j <= n; ++j)
        {
            tmp[j] = (tmp[j] + sign * f_[2*j + is_imag]) / (j * dy);
        }

        is_imag = !is_imag;
        sign = -sign;
    }

    // out[0] is done by direct integration
    // note that only the zeroth order spherical Bessel function is nonzero at 0
    if (l == 0)
    {
        for (int i = 0; i <= n; ++i)
        {
            tmp[0] += 2.0 * pref * in[i] * std::pow(i*dx, 2-p); // p <= 2 is required here
        }
    }

    // FFT-based method does not yield accurate results for small y at high l
    // use numerical integration in this case
    const int n_direct = (l == 0) ? 0 : static_cast<int>(ngrid * std::pow(1e-8, 1.0/l));
    if (n_direct > 0)
    {
        std::vector<double> buffer(ngrid + n_direct);
        double* grid_in = buffer.data();
        double* grid_out = grid_in + ngrid;

        std::for_each(grid_in, grid_in + ngrid,
            [&](double& x) { x = (&x - grid_in) * dx; });
        std::for_each(grid_out, grid_out + n_direct,
            [&](double& y) { y = ((&y - grid_out) + 1) * dy; });

        direct(l, ngrid, grid_in, in, n_direct, grid_out, &tmp[1], p);
    }

    std::copy(tmp.begin(), tmp.end(), out);

    if (!cache_enabled_)
    {
        _rfft_clear();
    }
}


void SphericalBesselTransformer::Impl::direct(
    const int l,
    const int ngrid_in,
    const double* const grid_in,
    const double* const in,
    const int ngrid_out,
    const double* const grid_out,
    double* const out,
    const int p
)
{
    assert(p <= 2);
    assert(ngrid_in > 1 && ngrid_out > 0);
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
        double* jl = &jl_[j * grid_in_.size()];
        std::transform(tmp, tmp + ngrid_in, jl, integrand, std::multiplies<double>());
        out[j] = ModuleBase::Integral::simpson(ngrid_in, integrand, &rab[1]);
    }

    const double pref = std::sqrt(2.0 / std::acos(-1.0));
    std::for_each(out, out + ngrid_out, [pref](double& x) { x *= pref; });

    if (!cache_enabled_)
    {
        _table_clear();
    }
}


void SphericalBesselTransformer::Impl::_rfft_prepare(const int n)
{
    if (n != sz_planned_)
    {
        if (f_)
        {
            fftw_free(f_);
        }
        f_ = fftw_alloc_real(sizeof(double) * 2 * (n/2 + 1));
        // see FFTW documentation "one-dimensional DFTs of real data"

        if (rfft_plan_)
        {
            fftw_destroy_plan(rfft_plan_);
        }
        auto* out = reinterpret_cast<fftw_complex*>(f_); // in-place transform
        rfft_plan_ = fftw_plan_dft_r2c_1d(n, f_, out, fftw_plan_flag_);

        sz_planned_ = n;
    }
}


void SphericalBesselTransformer::Impl::_tabulate(
    const int l,
    const int ngrid_in,
    const double* const grid_in,
    const int ngrid_out,
    const double* const grid_out
)
{
    const bool is_cached =
        cache_enabled_ && l == l_
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

    for (int j = 0; j < ngrid_out; ++j)
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

    sz_planned_ = 0;
}


void SphericalBesselTransformer::Impl::_table_clear()
{
    std::vector<double>().swap(grid_in_);
    std::vector<double>().swap(grid_out_);
    std::vector<double>().swap(jl_);
}


size_t SphericalBesselTransformer::Impl::heap_usage() const
{
    return (grid_in_.capacity() + grid_out_.capacity() + jl_.capacity()
            + sz_planned_) * sizeof(double);
}


void SphericalBesselTransformer::Impl::clear()
{
    _rfft_clear();
    _table_clear();
}


//**********************************************************************
//                          Interface
//**********************************************************************
SphericalBesselTransformer::SphericalBesselTransformer(const bool cache_enabled)
    : impl_(new Impl(cache_enabled))
{}


void SphericalBesselTransformer::radrfft(
    const int l,
    const int ngrid,
    const double cutoff,
    const double* const in,
    double* const out,
    const int p
) const
{
    impl_->radrfft(l, ngrid, cutoff, in, out, p);
}


void SphericalBesselTransformer::direct(
    const int l,
    const int ngrid_in,
    const double* const grid_in,
    const double* const in,
    const int ngrid_out,
    const double* const grid_out,
    double* const out,
    const int p
) const
{
    impl_->direct(l, ngrid_in, grid_in, in, ngrid_out, grid_out, out, p);
}


size_t SphericalBesselTransformer::heap_usage() const
{
    return impl_->heap_usage();
}


void SphericalBesselTransformer::clear()
{
    impl_->clear();
}


} // namespace ModuleBase
