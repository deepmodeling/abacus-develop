#include "module_base/spherical_bessel_transformer.h"

#include <algorithm>
#include <memory>

#include "gtest/gtest.h"
#include "module_base/constants.h"

#ifdef __MPI
#include <mpi.h>
#endif

using ModuleBase::PI;
using ModuleBase::SphericalBesselTransformer;

/***********************************************************
 *      Unit test of class SphericalBesselTransform
 ***********************************************************/

/*! Tested functions:
 *
 *  - radrfft
 *      - Performs a spherical Bessel transform via fast
 *        Fourier transforms.
 *
 *  - direct
 *      - Performs a spherical Bessel transform via quadrature
 *        using the Simpson's rule.
 *
 *  - set_fftw_plan_flag
 *      - Sets the planner flag for FFTW plan creation.
 *
 *                                                          */

class SphericalBesselTransformTest : public ::testing::Test
{

  protected:
    /// Allocates buffers
    void SetUp();

    /// Deallocates buffers
    void TearDown();

    /// Gets the maximum absolute element-wise difference between two arrays
    static double max_diff(int sz, const double* arr1, const double* arr2);

    const int sz_max = 10000;   ///< size of each buffer
    double* buffer = nullptr;   ///< buffer for all arrays below

    double* f = nullptr;        ///< input array
    double* g = nullptr;        ///< output array
    double* g_ref = nullptr;    ///< reference array
    double* grid_in = nullptr;  ///< input grid
    double* grid_out = nullptr; ///< output grid

    const double tol = 1e-9; ///< tolerance for element-wise numerical error
};

void SphericalBesselTransformTest::SetUp()
{
    buffer = new double[sz_max * 5];

    f = buffer;
    g = f + sz_max;
    g_ref = g + sz_max;
    grid_in = g_ref + sz_max;
    grid_out = grid_in + sz_max;
}

void SphericalBesselTransformTest::TearDown()
{
    delete[] buffer;
}

double SphericalBesselTransformTest::max_diff(int sz, const double* arr1, const double* arr2)
{
    double diff = 0.0;
    double tmp = 0.0;
    for (int i = 0; i < sz; ++i)
    {
        tmp = std::abs(arr1[i] - arr2[i]);
        if (tmp > diff)
        {
            diff = tmp;
        }
    }
    return diff;
}

TEST_F(SphericalBesselTransformTest, RadrfftBasic)
{
    /*
     * Computes the zeroth, first and second order spherical Bessel
     * transforms of r*exp(-r) and compares the results with analytic
     * expressions:
     *
     * zeroth:  2*sqrt(2/pi) * (3-k^2) / (k^2+1)^3.
     * first :  2*sqrt(2/pi) *   4k    / (k^2+1)^3.
     * second:  2*sqrt(2/pi) *   4k^2  / (k^2+1)^3.
     *                                                              */
    const int sz = 10000;
    assert(sz <= sz_max);

    const double dr = 0.01;
    const double rcut = dr * (sz - 1);
    const double dk = PI / rcut;
    const double pref = std::sqrt(2. / PI) * 2.;

    SphericalBesselTransformer sbt;

    for (int i = 0; i != sz; ++i)
    {
        double r = i * dr;
        f[i] = r * std::exp(-r);
    }

    // zeroth-order transform
    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * (3.0 - k * k) / std::pow(k * k + 1, 3);
    }
    sbt.radrfft(0, sz, rcut, f, g, 0);
    EXPECT_LT(max_diff(sz, g_ref, g), tol);

    // first-order transform
    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * 4.0 * k / std::pow(k * k + 1, 3);
    }
    sbt.radrfft(1, sz, rcut, f, g, 0);
    EXPECT_LT(max_diff(sz, g_ref, g), tol);

    // second-order transform
    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * 4.0 * k * k / std::pow(k * k + 1, 3);
    }
    sbt.radrfft(2, sz, rcut, f, g, 0);
    EXPECT_LT(max_diff(sz, g_ref, g), tol);
}

TEST_F(SphericalBesselTransformTest, RadrfftImplicitExponent)
{
    /*
     * Computes the second order spherical Bessel transform of
     * r^2*exp(-r) with input given as r^(p+2)*exp(-r) instead of
     * bare r^2*exp(-r). Compares the results with the analytic
     * expression:
     *
     *      48*sqrt(2/pi) * k^2  / (k^2+1)^4.
     *                                                          */
    const int sz = 5000;
    assert(sz <= sz_max);

    const double dr = 0.02;
    const double rcut = dr * (sz - 1);
    const double dk = PI / rcut;
    const double pref = std::sqrt(2. / PI) * 48.;

    SphericalBesselTransformer sbt(true);

    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * k * k / std::pow(k * k + 1, 4);
    }

    for (int p = -2; p <= 2; ++p)
    {
        for (int i = 0; i != sz; ++i)
        {
            double r = i * dr;
            f[i] = std::pow(r, 2 + p) * std::exp(-r);
        }
        sbt.radrfft(2, sz, rcut, f, g, p);
        EXPECT_LT(max_diff(sz, g_ref, g), tol);
    }
}

TEST_F(SphericalBesselTransformTest, RadrfftVariableSize)
{
    /*
     * Computes the second order spherical Bessel transform of
     * r^2*exp(-r) with various input sizes. Compares the results
     * with the analytic expression:
     *
     *      48*sqrt(2/pi) * k^2  / (k^2+1)^4.
     *                                                          */
    const double dr = 0.02;
    const double pref = std::sqrt(2. / PI) * 48.;

    SphericalBesselTransformer sbt;

    for (int sz = 5000; sz <= sz_max; sz += 1000)
    {

        for (int i = 0; i != sz; ++i)
        {
            double r = i * dr;
            f[i] = r * r * std::exp(-r);
        }

        const double rcut = dr * (sz - 1);
        const double dk = PI / rcut;
        for (int i = 0; i != sz; ++i)
        {
            double k = dk * i;
            g_ref[i] = pref * k * k / std::pow(k * k + 1, 4);
        }
        sbt.radrfft(2, sz, rcut, f, g, 0);
        EXPECT_LT(max_diff(sz, g_ref, g), tol);
    }
}

TEST_F(SphericalBesselTransformTest, RadrfftInPlace)
{
    /*
     * Performs an in-place second order spherical Bessel transform
     * on r^2*exp(-r^2). Compares the results with the analytic
     * expression:
     *
     *      sqrt(2)/16 * k^2 * exp(-k^2/4)
     *                                                          */
    const double dr = 0.02;
    const double pref = std::sqrt(2.) / 16.;

    SphericalBesselTransformer sbt;

    const double sz = 10000;
    const double rcut = dr * (sz - 1);
    for (int i = 0; i != sz; ++i)
    {
        double r = i * dr;
        f[i] = r * r * std::exp(-r * r);
    }

    double dk = PI / rcut;
    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * k * k * std::exp(-k * k / 4);
    }
    sbt.radrfft(2, sz, rcut, f, f, 0);
    EXPECT_LT(max_diff(sz, g_ref, f), tol);
}

TEST_F(SphericalBesselTransformTest, DirectBasic)
{
    /*
     * Computes the zeroth, first and second order spherical Bessel
     * transforms of r*exp(-r) and compares the results with analytic
     * expressions:
     *
     * zeroth:  2*sqrt(2/pi) * (3-k^2) / (k^2+1)^3.
     * first :  2*sqrt(2/pi) *   4k    / (k^2+1)^3.
     * second:  2*sqrt(2/pi) *   4k^2  / (k^2+1)^3.
     *                                                              */
    const int sz_in = 7000;
    const int sz_out = 5000;
    assert(sz_in <= sz_max && sz_out <= sz_max);

    const double dr = 0.007;
    const double dk = 0.013;
    std::for_each(grid_in, grid_in + sz_in, [&](double& x) { x = (&x - grid_in) * dr; });
    std::for_each(grid_out, grid_out + sz_out, [&](double& x) { x = (&x - grid_out) * dk; });

    const double pref = std::sqrt(2. / PI) * 2.;

    std::for_each(f, f + sz_in, [&](double& x) {
        double r = (&x - f) * dr;
        x = r * std::exp(-r);
    });

    // zeroth-order transform
    std::for_each(g_ref, g_ref + sz_out, [&](double& y) {
        double k = (&y - g_ref) * dk;
        y = pref * (3.0 - k * k) / std::pow(k * k + 1, 3);
    });

    SphericalBesselTransformer sbt;

    sbt.direct(0, sz_in, grid_in, f, sz_out, grid_out, g);
    EXPECT_LT(max_diff(sz_out, g_ref, g), tol);

    // first-order transform
    std::for_each(g_ref, g_ref + sz_out, [&](double& y) {
        double k = (&y - g_ref) * dk;
        y = pref * 4.0 * k / std::pow(k * k + 1, 3);
    });

    sbt.direct(1, sz_in, grid_in, f, sz_out, grid_out, g);
    EXPECT_LT(max_diff(sz_out, g_ref, g), tol);

    // second-order transform
    std::for_each(g_ref, g_ref + sz_out, [&](double& y) {
        double k = (&y - g_ref) * dk;
        y = pref * 4.0 * k * k / std::pow(k * k + 1, 3);
    });

    sbt.direct(2, sz_in, grid_in, f, sz_out, grid_out, g);
    EXPECT_LT(max_diff(sz_out, g_ref, g), tol);
}

TEST_F(SphericalBesselTransformTest, DirectImplicitExponent)
{
    /*
     * Computes the second order spherical Bessel transform of
     * r^2*exp(-r) with input given as r^(p+2)*exp(-r) instead of
     * bare r^2*exp(-r). Compares the results with the analytic
     * expression:
     *
     *      48*sqrt(2/pi) * k^2  / (k^2+1)^4.
     *                                                          */
    const int sz_in = 7000;
    const int sz_out = 6000;
    assert(sz_in <= sz_max && sz_out <= sz_max);

    const double dr = 0.007;
    const double dk = 0.011;
    std::for_each(grid_in, grid_in + sz_in, [&](double& x) { x = (&x - grid_in) * dr; });
    std::for_each(grid_out, grid_out + sz_out, [&](double& x) { x = (&x - grid_out) * dk; });

    const double pref = std::sqrt(2. / PI) * 48.;
    std::for_each(g_ref, g_ref + sz_out, [&](double& y) {
        double k = (&y - g_ref) * dk;
        y = pref * k * k / std::pow(k * k + 1, 4);
    });

    SphericalBesselTransformer sbt(true);

    for (int p = -2; p <= 2; ++p)
    {
        std::for_each(f, f + sz_in, [&](double& x) {
            double r = (&x - f) * dr;
            x = std::pow(r, 2 + p) * std::exp(-r);
        });

        sbt.direct(2, sz_in, grid_in, f, sz_out, grid_out, g, p);
        EXPECT_LT(max_diff(sz_out, g_ref, g), tol);
    }
}

TEST_F(SphericalBesselTransformTest, DirectInPlace)
{
    /*
     * Performs an in-place second order spherical Bessel transform
     * on r^2*exp(-r^2). Compares the results with the analytic
     * expression:
     *
     *      sqrt(2)/16 * k^2 * exp(-k^2/4)
     *                                                          */
    const int sz_in = 7000;
    const int sz_out = 7000;
    assert(sz_in <= sz_max && sz_out == sz_in);

    const double dr = 0.011;
    const double dk = 0.007;
    std::for_each(grid_in, grid_in + sz_in, [&](double& x) { x = (&x - grid_in) * dr; });
    std::for_each(grid_out, grid_out + sz_out, [&](double& x) { x = (&x - grid_out) * dk; });

    std::for_each(f, f + sz_in, [&](double& x) {
        double r = (&x - f) * dr;
        x = r * r * std::exp(-r * r);
    });

    const double pref = std::sqrt(2.) / 16.;
    std::for_each(g_ref, g_ref + sz_out, [&](double& y) {
        double k = (&y - g_ref) * dk;
        y = pref * k * k * std::exp(-k * k / 4);
    });

    SphericalBesselTransformer sbt;

    sbt.direct(2, sz_in, grid_in, f, sz_out, grid_out, f);
    EXPECT_LT(max_diff(sz_out, g_ref, f), tol);
}

TEST_F(SphericalBesselTransformTest, HighOrder)
{
    /*
     * Computes the l-order spherical Bessel transforms of
     * r^l*exp(-r^2) using radrfft and direct with some high l,
     * and check the consistency between their results.
     *                                                              */
    const int l = 6;
    const int sz = 5000;
    assert(sz <= sz_max);

    const double dr = 0.01;
    const double rcut = dr * (sz - 1);
    const double dk = PI / rcut;

    std::for_each(grid_in, grid_in + sz, [&](double& x) { x = (&x - grid_in) * dr; });
    std::for_each(grid_out, grid_out + sz, [&](double& x) { x = (&x - grid_out) * dk; });

    for (int i = 0; i != sz; ++i)
    {
        double r = i * dr;
        f[i] = std::pow(r, l) * std::exp(-r * r);
    }

    SphericalBesselTransformer sbt;
    sbt.radrfft(l, sz, rcut, f, g);

    sbt.direct(l, sz, grid_in, f, sz, grid_out, g_ref);

    // NOTE: Simpson's integration gets increasingly inaccurate as k gets large
    // since the factor of (k*dr)^4 in its error becomes significant when k*dr
    // is of order 1. So we only compare the results for relatively small k.
    EXPECT_LT(max_diff(sz / 2, g_ref, g), tol);
}

TEST_F(SphericalBesselTransformTest, HeapUsage)
{
    /*
     * Tests the setter and getter of the planner flag for FFTW plan
     * creation.
     *                                                              */
    const int sz = 2000;
    assert(sz <= sz_max);

    const double dr = 0.01;
    const double rcut = dr * (sz - 1);

    SphericalBesselTransformer sbt_no_cache(false);
    SphericalBesselTransformer sbt_cached(true);

    EXPECT_EQ(sbt_no_cache.heap_usage(), 0);
    EXPECT_EQ(sbt_cached.heap_usage(), 0);

    for (int i = 0; i != sz; ++i)
    {
        f[i] = std::exp(-i * dr);
    }

    sbt_no_cache.radrfft(0, sz, rcut, f, g, 0);
    sbt_cached.radrfft(0, sz, rcut, f, g, 0);

    // zeroth-order transform does not involve tabulating jl
    // but FFT needs to allocate memory
    EXPECT_EQ(sbt_no_cache.heap_usage(), 0);
    EXPECT_EQ(sbt_cached.heap_usage(), 2*(sz-1)*sizeof(double));

    // higher-order transforms involve tabulating jl
    sbt_no_cache.radrfft(5, sz, rcut, f, g, 0);
    sbt_cached.radrfft(5, sz, rcut, f, g, 0);

    EXPECT_EQ(sbt_no_cache.heap_usage(), 0);
    EXPECT_GT(sbt_cached.heap_usage(), 2*(sz-1)*sizeof(double));

    sbt_no_cache.clear();
    sbt_cached.clear();

    EXPECT_EQ(sbt_no_cache.heap_usage(), 0);
    EXPECT_EQ(sbt_cached.heap_usage(), 0);
}


int main(int argc, char** argv)
{

#ifdef __MPI
    int nprocs, id;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
#endif

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

#ifdef __MPI
    MPI_Finalize();
#endif

    fftw_cleanup();

    return result;
}
