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
    double max_diff(const int sz, const double* const arr1, const double* const arr2);

    int sz_max = 10000;      ///< size of each buffer
    double* f = nullptr;     ///< buffer for input array
    double* g = nullptr;     ///< buffer for output array
    double* g_ref = nullptr; ///< buffer for reference array

    double* grid_in = nullptr;  ///< buffer for input grid
    double* grid_out = nullptr; ///< buffer for output grid

    double tol = 1e-9; ///< tolerance for element-wise numerical error
};

void SphericalBesselTransformTest::SetUp()
{
    f = new double[sz_max];
    g = new double[sz_max];
    g_ref = new double[sz_max];
    grid_in = new double[sz_max];
    grid_out = new double[sz_max];
}

void SphericalBesselTransformTest::TearDown()
{
    delete[] f;
    delete[] g;
    delete[] g_ref;
    delete[] grid_in;
    delete[] grid_out;
}

double SphericalBesselTransformTest::max_diff(const int sz, const double* const arr1, const double* const arr2)
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
    int sz = 10000;
    assert(sz <= sz_max);

    double dr = 0.01;
    double rcut = dr * (sz - 1);
    double dk = PI / rcut;
    double pref = std::sqrt(2. / PI) * 2.;

    std::shared_ptr<SphericalBesselTransformer> sbt = SphericalBesselTransformer::create();

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
    sbt->radrfft(0, sz, rcut, f, g, 0);
    EXPECT_LT(max_diff(sz, g_ref, g), tol);

    // first-order transform
    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * 4.0 * k / std::pow(k * k + 1, 3);
    }
    sbt->radrfft(1, sz, rcut, f, g, 0);
    EXPECT_LT(max_diff(sz, g_ref, g), tol);

    // second-order transform
    for (int i = 0; i != sz; ++i)
    {
        double k = dk * i;
        g_ref[i] = pref * 4.0 * k * k / std::pow(k * k + 1, 3);
    }
    sbt->radrfft(2, sz, rcut, f, g, 0);
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
    int sz = 5000;
    assert(sz <= sz_max);

    double dr = 0.02;
    double rcut = dr * (sz - 1);
    double dk = PI / rcut;
    double pref = std::sqrt(2. / PI) * 48.;

    std::shared_ptr<SphericalBesselTransformer> sbt = SphericalBesselTransformer::create();
    sbt->set_fftw_plan_flag(FFTW_MEASURE);

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
        sbt->radrfft(2, sz, rcut, f, g, p);
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
    double dr = 0.02;
    double pref = std::sqrt(2. / PI) * 48.;

    std::shared_ptr<SphericalBesselTransformer> sbt = SphericalBesselTransformer::create();
    sbt->set_fftw_plan_flag(FFTW_ESTIMATE);

    for (int sz = 5000; sz <= sz_max; sz += 1000)
    {

        double rcut = dr * (sz - 1);
        for (int i = 0; i != sz; ++i)
        {
            double r = i * dr;
            f[i] = r * r * std::exp(-r);
        }

        double dk = PI / rcut;
        for (int i = 0; i != sz; ++i)
        {
            double k = dk * i;
            g_ref[i] = pref * k * k / std::pow(k * k + 1, 4);
        }
        sbt->radrfft(2, sz, rcut, f, g, 0);
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
    double dr = 0.02;
    double pref = std::sqrt(2.) / 16.;

    std::shared_ptr<SphericalBesselTransformer> sbt = SphericalBesselTransformer::create();
    sbt->set_fftw_plan_flag(FFTW_ESTIMATE);

    double sz = 10000;
    double rcut = dr * (sz - 1);
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
    sbt->radrfft(2, sz, rcut, f, f, 0);
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
    int sz_in = 7000;
    int sz_out = 5000;
    assert(sz_in <= sz_max && sz_out <= sz_max);

    double dr = 0.007;
    double dk = 0.013;
    std::for_each(grid_in, grid_in + sz_in, [&](double& x) { x = (&x - grid_in) * dr; });
    std::for_each(grid_out, grid_out + sz_out, [&](double& x) { x = (&x - grid_out) * dk; });

    double pref = std::sqrt(2. / PI) * 2.;

    std::for_each(f, f + sz_in, [&](double& x) {
        double r = (&x - f) * dr;
        x = r * std::exp(-r);
    });

    // zeroth-order transform
    std::for_each(g_ref, g_ref + sz_out, [&](double& y) {
        double k = (&y - g_ref) * dk;
        y = pref * (3.0 - k * k) / std::pow(k * k + 1, 3);
    });

    std::shared_ptr<SphericalBesselTransformer> sbt = SphericalBesselTransformer::create();

    sbt->direct(0, sz_in, grid_in, f, sz_out, grid_out, g);
    EXPECT_LT(max_diff(sz_out, g_ref, g), tol);

    // first-order transform
    std::for_each(g_ref, g_ref + sz_out, [&](double& y) {
        double k = (&y - g_ref) * dk;
        y = pref * 4.0 * k / std::pow(k * k + 1, 3);
    });

    sbt->direct(1, sz_in, grid_in, f, sz_out, grid_out, g);
    EXPECT_LT(max_diff(sz_out, g_ref, g), tol);

    // second-order transform
    std::for_each(g_ref, g_ref + sz_out, [&](double& y) {
        double k = (&y - g_ref) * dk;
        y = pref * 4.0 * k * k / std::pow(k * k + 1, 3);
    });

    sbt->direct(2, sz_in, grid_in, f, sz_out, grid_out, g);
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
    int sz_in = 7000;
    int sz_out = 6000;
    assert(sz_in <= sz_max && sz_out <= sz_max);

    double dr = 0.007;
    double dk = 0.011;
    std::for_each(grid_in, grid_in + sz_in, [&](double& x) { x = (&x - grid_in) * dr; });
    std::for_each(grid_out, grid_out + sz_out, [&](double& x) { x = (&x - grid_out) * dk; });

    double pref = std::sqrt(2. / PI) * 48.;
    std::for_each(g_ref, g_ref + sz_out, [&](double& y) {
        double k = (&y - g_ref) * dk;
        y = pref * k * k / std::pow(k * k + 1, 4);
    });

    std::shared_ptr<SphericalBesselTransformer> sbt = SphericalBesselTransformer::create();

    for (int p = -2; p <= 2; ++p)
    {
        std::for_each(f, f + sz_in, [&](double& x) {
            double r = (&x - f) * dr;
            x = std::pow(r, 2 + p) * std::exp(-r);
        });

        sbt->direct(2, sz_in, grid_in, f, sz_out, grid_out, g, p);
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
    int sz_in = 7000;
    int sz_out = 7000;
    assert(sz_in <= sz_max && sz_out == sz_in);

    double dr = 0.011;
    double dk = 0.007;
    std::for_each(grid_in, grid_in + sz_in, [&](double& x) { x = (&x - grid_in) * dr; });
    std::for_each(grid_out, grid_out + sz_out, [&](double& x) { x = (&x - grid_out) * dk; });

    std::for_each(f, f + sz_in, [&](double& x) {
        double r = (&x - f) * dr;
        x = r * r * std::exp(-r * r);
    });

    double pref = std::sqrt(2.) / 16.;
    std::for_each(g_ref, g_ref + sz_out, [&](double& y) {
        double k = (&y - g_ref) * dk;
        y = pref * k * k * std::exp(-k * k / 4);
    });

    std::shared_ptr<SphericalBesselTransformer> sbt = SphericalBesselTransformer::create();

    sbt->direct(2, sz_in, grid_in, f, sz_out, grid_out, f);
    EXPECT_LT(max_diff(sz_out, g_ref, f), tol);
}

TEST_F(SphericalBesselTransformTest, HighOrder)
{
    /*
     * Computes the l-order spherical Bessel transforms of
     * r^l*exp(-r^2) using radrfft and direct with some high l,
     * and check the consistency between their results.
     *                                                              */
    int l = 6;
    int sz = 5000;
    assert(sz <= sz_max);

    double dr = 0.01;
    double rcut = dr * (sz - 1);
    double dk = PI / rcut;

    std::for_each(grid_in, grid_in + sz, [&](double& x) { x = (&x - grid_in) * dr; });
    std::for_each(grid_out, grid_out + sz, [&](double& x) { x = (&x - grid_out) * dk; });

    for (int i = 0; i != sz; ++i)
    {
        double r = i * dr;
        f[i] = std::pow(r, l) * std::exp(-r * r);
    }

    std::shared_ptr<SphericalBesselTransformer> sbt = SphericalBesselTransformer::create();
    sbt->radrfft(l, sz, rcut, f, g);

    sbt->direct(l, sz, grid_in, f, sz, grid_out, g_ref);

    // NOTE: Simpson's integration gets increasingly inaccurate as k gets large
    // since the factor of (k*dr)^4 in its error becomes significant when k*dr
    // is of order 1. So we only compare the results for relatively small k.
    EXPECT_LT(max_diff(sz / 2, g_ref, g), tol);
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
