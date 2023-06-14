#include "module_base/spherical_bessel_transformer.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "module_base/constants.h"

namespace ModuleBase
{

SphericalBesselTransformer::~SphericalBesselTransformer()
{
    fftw_destroy_plan(rfft_plan_);
    fftw_free(f_);
}

int SphericalBesselTransformer::spherical_bessel_sincos_polycoef(const bool get_sine, const int l, const int n)
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
    assert(l <= 10 && "SphericalBesselTransformer::spherical_bessel_sincos_polycoef: \
                some coefficients exceed INT_MAX (2147483647) for l >= 11");

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

    int n = ngrid - 1;
    double dx = cutoff / n;
    double dy = PI / cutoff;
    double pref = dx / std::sqrt(2. * PI);

    double* out_tmp = new double[ngrid];
    for (int i = 0; i != n + 1; ++i)
    {
        out_tmp[i] = 0.0;
    }

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
        for (int i = 1; i <= n; ++i)
        {
            out_tmp[0] += 2. * pref * in[i] * std::pow(i * dx, 2 - p);
        }

        if (p == 2)
        {
            out_tmp[0] += 2. * pref * in[0];
        }

        if (p > 2)
        {
            // estimate F(0) from the next three points by a quadratic interpolation
            // FIXME this may not be a good numerical solution, to be studied later
            out_tmp[0] += 3 * in[1] / std::pow(dx, p) - 3 * in[2] / std::pow(2 * dx, p) + in[3] / std::pow(3 * dx, p);
        }
    }

    for (int j = 0; j <= n; ++j)
    {
        out[j] = out_tmp[j];
    }

    delete[] out_tmp;
}

void SphericalBesselTransformer::direct(const int l,
                                        const int nin,
                                        const double* const igrid,
                                        const double* const in,
                                        const int nout,
                                        const double* const ogrid,
                                        double* const out,
                                        const int p)
{
    // TBD
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
