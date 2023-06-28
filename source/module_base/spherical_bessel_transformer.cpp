#include "module_base/spherical_bessel_transformer.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <numeric>

#include "module_base/constants.h"
#include "module_base/math_integral.h"
#include "module_base/math_sphbes.h"

namespace ModuleBase
{

SphericalBesselTransformer::~SphericalBesselTransformer()
{
    fftw_destroy_plan(rfft_plan_);
    fftw_free(f_);
}

long long int SphericalBesselTransformer::spherical_bessel_sincos_polycoef(const bool get_sine,
                                                                           const int l,
                                                                           const int n)
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
    assert(l <= 17 && "SphericalBesselTransformer::spherical_bessel_sincos_polycoef: \
                some coefficients exceed LLONG_MAX (2^63-1) for l >= 18");

    if (get_sine)
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

    return (2 * l - 1) * spherical_bessel_sincos_polycoef(get_sine, l - 1, n)
           - (n >= 2 ? spherical_bessel_sincos_polycoef(get_sine, l - 2, n - 2) : 0);
}

void SphericalBesselTransformer::radrfft(const int l,
                                         const int ngrid,
                                         const double cutoff,
                                         const double* const in,
                                         double* const out,
                                         const int p)
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
     * c(l,n) / s(l,n) are sin / cos coefficients from spherical_bessel_sincos_polycoef, and
     * the domain of F(x) is extended to negative values by defining F(-x) = pow(-1,l)*F(x).
     *
     * With an appropriate grid, the Fourier transform can be approximated by a discrete Fourier transform.
     * Note that given l & n, c(l,n) and s(l,n) cannot be both non-zero. Therefore, each FFT input
     * array is either purely real or purely imaginary, which suggests the use of real-input FFT.
     *                                                                                                  */
    assert(l >= 0);
    assert(ngrid > 1);
    assert(p <= 2);

    int n = ngrid - 1;
    double dx = cutoff / n;
    double dy = PI / cutoff;
    double pref = dx / std::sqrt(2. * PI);

    // use a temporary output array in case the transform is in-place
    double* out_tmp = new double[ngrid];
    std::fill(out_tmp, out_tmp + ngrid, 0.0);

    // rfft buffer (f_) allocation and plan creation
    rfft_prepare(2 * n);

    // introduced to help manipulate f_ ( double(*)[2] ) as an array of double
    double* f = &f_[0][0];

    for (int m = 0; m <= l; ++m)
    {

        // m even --> sin; f[2*n-i] = -f[i]; out += -imag(rfft(f)) / y^(l+1-m)
        // m odd  --> cos; f[2*n-i] = +f[i]; out += +real(rfft(f)) / y^(l+1-m)
        bool flag = (m % 2 == 0);

        int sincos_coef = spherical_bessel_sincos_polycoef(flag, l, m);

        f[0] = f[n] = 0.0;
        for (int i = 1; i != n; ++i)
        {
            f[i] = pref * sincos_coef * in[i] * std::pow(i * dx, m + 1 - l - p);
            f[2 * n - i] = flag ? -f[i] : f[i];
        }

        rfft_in_place();

        // summing up the series by ( ... ( ( g0/y + g1 )/y + g2 )/y + ... + gl )/y
        // out[0] is handled independently
        for (int j = 1; j <= n; ++j)
        {
            out_tmp[j] = (out_tmp[j] + (flag ? -f_[j][1] : f_[j][0])) / (j * dy);
        }
    }

    // out[0] is done by direct integration
    // note that only the zeroth order spherical Bessel function is nonzero at 0
    if (l == 0)
    {
        for (int i = 0; i <= n; ++i)
        {
            out_tmp[0] += 2. * pref * in[i] * std::pow(i * dx, 2 - p); // p <= 2 is required!
        }
    }

    // FFT-based method is not accurate for small k at high l
    // use numerical integration in this case
    int n_direct = (l == 0) ? 0 : static_cast<int>(ngrid * std::pow(1e-8, 1. / l));
    if (n_direct > 0)
    {
        double* grid_in = new double[ngrid];
        double* grid_out = new double[n_direct];
        std::for_each(grid_in, grid_in + ngrid, [&](double& x) { x = (&x - grid_in) * dx; });
        std::for_each(grid_out, grid_out + n_direct, [&](double& y) { y = ((&y - grid_out) + 1) * dy; });
        direct(l, ngrid, grid_in, in, n_direct, grid_out, &out_tmp[1], p);
        delete[] grid_in;
        delete[] grid_out;
    }

    std::memcpy(out, out_tmp, ngrid * sizeof(double));
    delete[] out_tmp;
}

void SphericalBesselTransformer::direct(const int l,
                                        const int ngrid_in,
                                        const double* const grid_in,
                                        const double* const in,
                                        const int ngrid_out,
                                        const double* const grid_out,
                                        double* const out,
                                        const int p)
{
    assert(p <= 2);
    assert(grid_in[0] >= 0.0 && std::is_sorted(grid_in, grid_in + ngrid_in, std::less_equal<double>()));
    assert(grid_out[0] >= 0.0 && std::is_sorted(grid_out, grid_out + ngrid_out, std::less_equal<double>()));

    double pref = std::sqrt(2.0 / ModuleBase::PI);

    // NOTE: std::adjacent_difference returns an array of the same size as the input array
    // with rab[0] = grid_in[0] and rab[i] = grid_in[i] - grid[i-1] for i > 0
    double* rab = new double[ngrid_in];
    std::adjacent_difference(grid_in, grid_in + ngrid_in, rab);

    // r^2*f(r) (integrand without j_l)
    double* r2f = new double[ngrid_in];
    std::memcpy(r2f, in, ngrid_in * sizeof(double));
    std::for_each(r2f, r2f + ngrid_in, [&](double& x) { x *= std::pow(grid_in[&x - r2f], 2 - p); });

    // there's no need to use a temporary output array even if the transform is in-place
    // because "in" is not referenced after the construction of r2f.

    double* integrand = new double[ngrid_in];

    // ModuleBase::Sphbes::Spherical_Bessel has potential risk for small function arguments
    double* jl = new double[ngrid_in];
    for (int i_out = 0; i_out != ngrid_out; ++i_out)
    {
        ModuleBase::Sphbes::Spherical_Bessel(ngrid_in, grid_in, grid_out[i_out], l, jl);
        std::transform(r2f, r2f + ngrid_in, jl, integrand, std::multiplies<double>());
        out[i_out] = ModuleBase::Integral::simpson(ngrid_in, integrand, &rab[1]);
    }
    delete[] jl;

    // for (int i_out = 0; i_out != ngrid_out; ++i_out)
    //{
    //     std::transform(grid_in, grid_in + ngrid_in, r2f, integrand, [&](double x, double v) {
    //         return ModuleBase::Sphbes::sphbesj(l, grid_out[i_out] * x) * v;
    //     });
    //     out[i_out] = ModuleBase::Integral::simpson(ngrid_in, integrand, &rab[1]);
    // }

    std::for_each(out, out + ngrid_out, [pref](double& x) { x *= pref; });

    delete[] r2f;
    delete[] integrand;
    delete[] rab;
}

void SphericalBesselTransformer::rfft_in_place()
{
    fftw_execute(rfft_plan_);
}

void SphericalBesselTransformer::rfft_prepare(const int sz)
{
    if (sz != sz_planned_)
    {
        // if a new FFT size is requested,
        // reallocate the buffer and recreate the plan.
        fftw_free(f_);
        f_ = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (sz / 2 + 1));
        sz_planned_ = sz;
        fftw_destroy_plan(rfft_plan_);
        rfft_plan_ = fftw_plan_dft_r2c_1d(sz, &f_[0][0], f_, fftw_plan_flag_);
    }
    else if (!rfft_plan_)
    {
        // if the existing buffer already has the requested size but no plan is available
        rfft_plan_ = fftw_plan_dft_r2c_1d(sz, &f_[0][0], f_, fftw_plan_flag_);
    }
}

void SphericalBesselTransformer::set_fftw_plan_flag(const unsigned new_flag)
{
    assert(new_flag == FFTW_ESTIMATE || new_flag == FFTW_MEASURE);
    if (new_flag != fftw_plan_flag_)
    {
        fftw_destroy_plan(rfft_plan_);
        rfft_plan_ = nullptr;
        fftw_plan_flag_ = new_flag;
    }
}

} // namespace ModuleBase
