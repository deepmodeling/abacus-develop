#include "source_base/constants.h"
#include "source_base/math_integral.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace
{
constexpr double doublethreshold = 1e-12;

TEST(GaussLegendre, PolynomialMoments)
{
    for (int n : {1, 2, 3, 5, 16, 32, 64, 128, 140, 256})
    {
        std::vector<double> x(n);
        std::vector<double> w(n);
        ModuleBase::Integral::Gauss_Legendre_grid_and_weight(n, x.data(), w.data());
        for (int k = 0; k < 2 * n; ++k)
        {
            long double sum = 0;
            for (int i = 0; i < n; ++i)
            {
                const long double node = x[i];
                sum += static_cast<long double>(w[i]) * std::pow(node, k);
            }
            const long double exact = k % 2 ? 0 : 2.L / (k + 1);
            EXPECT_NEAR(sum, exact, 8e-15L) << "n=" << n << " power=" << k;
        }
    }
}

TEST(GaussLegendre, MappedIntervalMoments)
{
    for (const std::array<double, 2> &interval : {std::array<double, 2>{0, 4}, {2, 6}, {4, 0}})
    {
        const double lo = interval[0];
        const double hi = interval[1];
        const long double half = (hi - lo) / 2;
        const long double center = (hi + lo) / 2;
        for (int n : {1, 3, 5, 16, 32, 64, 128, 140})
        {
            std::vector<double> x(n);
            std::vector<double> w(n);
            ModuleBase::Integral::Gauss_Legendre_grid_and_weight(lo, hi, n, x.data(), w.data());
            for (int k = 0; k < 2 * n; ++k)
            {
                long double sum = 0;
                for (int i = 0; i < n; ++i)
                {
                    const long double t = (static_cast<long double>(x[i]) - center) / half;
                    sum += static_cast<long double>(w[i]) * std::pow(t, k);
                }
                const long double exact = k % 2 ? 0 : 2 * half / (k + 1);
                EXPECT_NEAR(sum, exact, 2e-14L) << "n=" << n << " power=" << k;
            }
        }
    }
}

TEST(GaussLegendre, WeightUsesReturnedNodeDerivative)
{
    for (int n : {3, 5, 8, 16, 32})
    {
        std::vector<double> x(n);
        std::vector<double> w(n);
        ModuleBase::Integral::Gauss_Legendre_grid_and_weight(n, x.data(), w.data());
        for (int i = 0; i < n; ++i)
        {
            const long double z = x[i];
            long double p = 1;
            long double previous = 0;
            for (int j = 1; j <= n; ++j)
            {
                const long double older = previous;
                previous = p;
                p = ((2 * j - 1) * z * previous - (j - 1) * older) / j;
            }
            const long double derivative = n * (previous - z * p) / (1 - z * z);
            const long double exact = 2 / ((1 - z * z) * derivative * derivative);
            const long double tolerance = 8e-14L * exact;
            EXPECT_NEAR(w[i], exact, tolerance) << "n=" << n << " node=" << i;
        }
    }
}

/************************************************
 *  unit test of class Integral
 ***********************************************/

/**
 * Tested functions:
 *  - function Simpson_Integral with rab as input
 *  - function Simpson_Integral with dr as input
 *  - function Simpson_Integral_0toall
 *  - function Simpson_Integral_alltoinf
 *  - function simpson with f(x) = sin(x)
 *  - function simpson with f(x) = 1 / (1 + x^2)
 *  - function simpson with f(x) = exp(x)
 */

// generate irregular grid with sinx
void sinspace(
    double start,
    double end,
    const int nums,
    double *xv,
    double *h)
{
    double astart = std::asin(start);
    double aend = std::asin(end);
    double step = (aend - astart) / (nums - 1);

    for (int i = 0; i < nums; ++i)
    {
        const double angle = astart + i * step;
        h[i] = std::sin(angle);
    }
    // calculate the difference
    xv[0] = start;
    for (int i = 0; i < nums - 1; ++i)
    {
        h[i] = h[i + 1] - h[i];
        xv[i + 1] = xv[i] + h[i];
    }
}

TEST(LegacySimpson, SinglePointUniformGrid)
{
    for (double value : {0.0, 7.0, -3.0})
    {
        double integral = 123.0;
        ModuleBase::Integral::Simpson_Integral(1, &value, 0.2, integral);
        EXPECT_DOUBLE_EQ(integral, 0.0);
    }
}

TEST(LegacySimpson, SinglePointRadialJacobian)
{
    for (double value : {0.0, 7.0, -3.0})
    {
        const double jacobian = 0.2;
        double integral = 123.0;
        ModuleBase::Integral::Simpson_Integral(1, &value, &jacobian, integral);
        EXPECT_DOUBLE_EQ(integral, 0.0);
    }
}

class SimpsonIntegralSinx : public testing::Test
{
    /**
     * test the integral of sinx between [0,PI],
     * divide into mesh-1 intervals
     */

  protected:
    std::vector<double> func;
    std::vector<double> rab;
    int mesh = 10001;
    double dr = ModuleBase::PI / (mesh - 1);
    double asum;
    std::vector<double> asumlist;
    double expectvalue = 2.0;

    void SetUp()
    {
        func.resize(mesh);
        rab.resize(mesh);
        asumlist.resize(mesh);

        for (int i = 0; i < mesh; i++)
        {
            const double angle = ModuleBase::PI * i / (mesh - 1);
            func[i] = std::sin(angle);
            rab[i] = dr;
        }
    }
};

TEST_F(SimpsonIntegralSinx, Constructor) { EXPECT_NO_THROW(ModuleBase::Integral intgral); }

TEST_F(SimpsonIntegralSinx, SimpsonIntegralRab)
{
    ModuleBase::Integral::Simpson_Integral(mesh, func.data(), rab.data(), asum);
    EXPECT_NEAR(asum, expectvalue, doublethreshold);
    const int even_mesh = mesh - 1;
    EXPECT_DEATH(ModuleBase::Integral::Simpson_Integral(even_mesh, func.data(), rab.data(), asum), "");
}

TEST_F(SimpsonIntegralSinx, SimpsonIntegralDr)
{
    ModuleBase::Integral::Simpson_Integral(mesh, func.data(), dr, asum);
    EXPECT_NEAR(asum, expectvalue, doublethreshold);
    const int even_mesh = mesh - 1;
    EXPECT_DEATH(ModuleBase::Integral::Simpson_Integral(even_mesh, func.data(), dr, asum), "");
}

TEST_F(SimpsonIntegralSinx, SimpsonIntegral0toall)
{
    const int halfmesh = mesh / 2;
    const double expected_half = expectvalue / 2.0;
    ModuleBase::Integral::Simpson_Integral_0toall(mesh, func.data(), rab.data(), asumlist.data());
    EXPECT_NEAR(asumlist[mesh - 1], expectvalue, doublethreshold);
    EXPECT_NEAR(asumlist[mesh - 2], 1.9999999506519738901, doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh], expected_half, doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh - 1], 0.99968584073722277505, doublethreshold);
}

TEST_F(SimpsonIntegralSinx, SimpsonIntegralAlltoinf)
{
    const int halfmesh = mesh / 2;
    const double expected_half = expectvalue / 2.0;
    ModuleBase::Integral::Simpson_Integral_alltoinf(mesh, func.data(), rab.data(), asumlist.data());
    EXPECT_NEAR(asumlist[0], expectvalue, doublethreshold);
    EXPECT_NEAR(asumlist[1], 1.9999999506519754444, doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh], expected_half, doublethreshold);
    EXPECT_NEAR(asumlist[halfmesh - 1], 1.0003141592627740053, doublethreshold);
}

TEST(CumulativeSimpson, LinearFunctionOnOddAndEvenMeshes)
{
    for (const int n : {1, 2, 3, 4, 5, 6, 9, 10})
    {
        const double h = 0.125;
        std::vector<double> f(n);
        std::vector<double> rab(n, h);
        std::vector<double> prefix(n);
        std::vector<double> tail(n);
        for (int i = 0; i < n; ++i)
        {
            f[i] = 1.0 + 2.0 * i * h;
        }
        ModuleBase::Integral::Simpson_Integral_0toall(n, f.data(), rab.data(), prefix.data());
        ModuleBase::Integral::Simpson_Integral_alltoinf(n, f.data(), rab.data(), tail.data());
        const double end = (n - 1) * h;
        const double total = end + end * end;
        for (int i = 0; i < n; ++i)
        {
            const double x = i * h;
            const double expected_prefix = x + x * x;
            const double expected_tail = total - x - x * x;
            EXPECT_NEAR(prefix[i], expected_prefix, 1e-14);
            EXPECT_NEAR(tail[i], expected_tail, 1e-14);
        }
    }
}

TEST(CumulativeSimpson, RadialJacobianIsAppliedOnce)
{
    // r(i)=i+i^2/4, so a constant radial integrand has a linear index-space integrand.
    for (const int n : {2, 4, 7, 10})
    {
        std::vector<double> f(n, 1.0);
        std::vector<double> rab(n);
        std::vector<double> prefix(n);
        std::vector<double> tail(n);
        for (int i = 0; i < n; ++i)
        {
            rab[i] = 1.0 + 0.5 * i;
        }
        ModuleBase::Integral::Simpson_Integral_0toall(n, f.data(), rab.data(), prefix.data());
        ModuleBase::Integral::Simpson_Integral_alltoinf(n, f.data(), rab.data(), tail.data());
        const double end = (n - 1) + 0.25 * (n - 1) * (n - 1);
        for (int i = 0; i < n; ++i)
        {
            const double r = i + 0.25 * i * i;
            EXPECT_DOUBLE_EQ(prefix[i], r);
            const double expected_tail = end - r;
            EXPECT_DOUBLE_EQ(tail[i], expected_tail);
        }
    }
}

TEST(CumulativeSimpson, LastOddEndpointKeepsTrapezoidalRule)
{
    const double f[] = {0.0, 1.0, 4.0, 9.0};
    const double rab[] = {1.0, 1.0, 1.0, 1.0};
    double prefix[4];
    ModuleBase::Integral::Simpson_Integral_0toall(4, f, rab, prefix);
    EXPECT_DOUBLE_EQ(prefix[0], 0.0);
    EXPECT_DOUBLE_EQ(prefix[1], 0.5);
    const double expected_simpson = 8.0 / 3.0;
    const double expected_final = expected_simpson + 6.5;
    EXPECT_DOUBLE_EQ(prefix[2], expected_simpson);
    EXPECT_DOUBLE_EQ(prefix[3], expected_final);
}

TEST_F(SimpsonIntegralSinx, UniformGridOdd)
{
    double start = 0.0;
    int ngrid_max = std::pow(3, 7);
    std::vector<double> f(ngrid_max);
    for (int ngrid = 3; ngrid <= ngrid_max; ngrid *= 3)
    {
        for (double end = 0.1; end < 10 * ModuleBase::PI; end += 0.5)
        {
            double dx = (end - start) / (ngrid - 1);
            for (int i = 0; i != ngrid; ++i)
            {
                const double x = start + i * dx;
                f[i] = std::sin(x);
            }
            // error bound of composite Simpson's rule (max|sin^(4)(x)| = 1)
            double tol = (end - start) * std::pow(dx, 4) / 180;
            const double expected = std::cos(start) - std::cos(end);
            const double actual = ModuleBase::Integral::simpson(ngrid, f.data(), dx);
            const double tolerance = std::max(tol, doublethreshold);
            EXPECT_NEAR(expected, actual, tolerance);
        }
    }
}

TEST_F(SimpsonIntegralSinx, UniformGridEven)
{
    double start = 0.0;
    int ngrid_max = std::pow(2, 12);
    std::vector<double> f(ngrid_max);
    for (int ngrid = 4; ngrid <= ngrid_max; ngrid *= 2)
    {
        for (double end = 0.1; end < 10 * ModuleBase::PI; end += 0.5)
        {
            double dx = (end - start) / (ngrid - 1);
            for (int i = 0; i != ngrid; ++i)
            {
                const double x = start + i * dx;
                f[i] = std::sin(x);
            }
            // error bound of composite Simpson's rule (max|sin^(4)(x)| = 1)
            double tol =
                ((ngrid - 4) * dx - start) * std::pow(dx, 4) / 180 + (end - (ngrid - 4) * dx) * std::pow(dx, 4) / 80;
            const double expected = std::cos(start) - std::cos(end);
            const double actual = ModuleBase::Integral::simpson(ngrid, f.data(), dx);
            const double tolerance = std::max(tol, doublethreshold);
            EXPECT_NEAR(expected, actual, tolerance);
        }
    }
}

TEST_F(SimpsonIntegralSinx, LogGridOdd)
{
    double start = 0.0;
    int ngrid_max = std::pow(3, 7);
    std::vector<double> f(ngrid_max);
    std::vector<double> h(ngrid_max);
    for (int ngrid = 3; ngrid <= ngrid_max; ngrid *= 3)
    {
        for (double w = 0.1; w < 10 * ModuleBase::PI; w += 0.5)
        {
            const double base = 1.0 + w;
            const double exponent = 1.0 / (ngrid - 1);
            double b = std::pow(base, exponent);
            for (int i = 0; i != ngrid; ++i)
            {
                double x = std::pow(b, i) - 1.0 + start;
                f[i] = std::sin(x);
                const int next_index = i + 1;
                h[i] = std::pow(b, next_index) - std::pow(b, i); // the last one is not used
            }
            // FIXME error bound for irregularly-spaced Simpson's rule should be derived
            // the one below is a crude estimate
            double tol = w * std::pow(h[ngrid - 2], 4);
            const double end = start + w;
            const double expected = std::cos(start) - std::cos(end);
            const double actual = ModuleBase::Integral::simpson(ngrid, f.data(), h.data());
            const double tolerance = std::max(tol, doublethreshold);
            EXPECT_NEAR(expected, actual, tolerance);
        }
    }
}

TEST_F(SimpsonIntegralSinx, LogGridEven)
{
    double start = 0.0;
    int ngrid_max = std::pow(2, 12);
    std::vector<double> f(ngrid_max);
    std::vector<double> h(ngrid_max);
    for (int ngrid = 4; ngrid <= ngrid_max; ngrid *= 2)
    {
        for (double w = 0.1; w < 10 * ModuleBase::PI; w += 0.5)
        {
            const double base = 1.0 + w;
            const double exponent = 1.0 / (ngrid - 1);
            double b = std::pow(base, exponent);
            for (int i = 0; i != ngrid; ++i)
            {
                double x = std::pow(b, i) - 1.0 + start;
                f[i] = std::sin(x);
                const int next_index = i + 1;
                h[i] = std::pow(b, next_index) - std::pow(b, i); // the last one is not used
            }
            // FIXME error bound for irregularly-spaced Simpson's rule should be derived
            // the one below is a crude estimate
            double tol = w * std::pow(h[ngrid - 2], 4);
            const double end = start + w;
            const double expected = std::cos(start) - std::cos(end);
            const double actual = ModuleBase::Integral::simpson(ngrid, f.data(), h.data());
            const double tolerance = std::max(tol, doublethreshold);
            EXPECT_NEAR(expected, actual, tolerance);
        }
    }
}

TEST_F(SimpsonIntegralSinx, FourPoints)
{
    const double constant[4] = {1.0, 1.0, 1.0, 1.0};
    const double intervals[3] = {1.0, 2.0, 3.0};
    const double constant_integral = ModuleBase::Integral::simpson(4, constant, intervals);
    EXPECT_NEAR(6.0, constant_integral, 1e-12);
    /*!
     * This test checks whether "simpson" yields the exact integral for a cubic polynomial
     * sampled by four points. This case is handled by Simpson's 3/8 rule and should be exact.
     *                                                                                              */
    double A = 0.1;
    double B = 0.25;
    double C = 0.40;
    double D = 0.77;
    double x0 = 0.52;

    auto f = [&](double x)
    {
        const double offset = x - x0;
        return A * std::pow(offset, 3) + B * std::pow(offset, 2) + C * (offset) + D;
    };
    auto I = [&](double x)
    {
        const double offset = x - x0;
        return A * std::pow(offset, 4) / 4 + B * std::pow(offset, 3) / 3 + C * std::pow(offset, 2) / 2 + D * (offset);
    };

    double x[4];
    double y[4];

    // unevenly-spaced 4 points
    x[0] = 0.11;
    x[1] = 0.22;
    x[2] = 0.44;
    x[3] = 0.88;

    for (int i = 0; i != 4; ++i)
    {
        y[i] = f(x[i]);
    }

    double ref = I(x[3]) - I(x[0]);
    double h[3] = {x[1] - x[0], x[2] - x[1], x[3] - x[2]};
    const double irregular_integral = ModuleBase::Integral::simpson(4, y, h);
    EXPECT_DOUBLE_EQ(ref, irregular_integral);

    // evenly-spaced 4 points
    x[0] = 0.1;
    x[1] = 0.2;
    x[2] = 0.3;
    x[3] = 0.4;
    for (int i = 0; i != 4; ++i)
    {
        y[i] = f(x[i]);
    }

    ref = I(x[3]) - I(x[0]);
    const double dx = x[1] - x[0];
    const double uniform_integral = ModuleBase::Integral::simpson(4, y, dx);
    EXPECT_DOUBLE_EQ(ref, uniform_integral);
}

TEST_F(SimpsonIntegralSinx, ThreePoints)
{
    /*!
     * This test checks whether "simpson" yields the exact integral for a quadratic polynomial
     * sampled by three points.
     *                                                                                              */
    double A = 0.17;
    double B = 0.29;
    double C = 0.35;
    double x0 = 0.5;

    auto f = [&](double x)
    {
        const double offset = x - x0;
        return A * std::pow(offset, 2) + B * (offset) + C;
    };
    auto I = [&](double x)
    {
        const double offset = x - x0;
        return A * std::pow(offset, 3) / 3 + B * std::pow(offset, 2) / 2 + C * (offset);
    };

    double x[3];

    // unevenly-spaced 3 points
    x[0] = 0.1;
    x[1] = 0.2;
    x[2] = 0.4;

    double y[3];
    y[0] = f(x[0]);
    y[1] = f(x[1]);
    y[2] = f(x[2]);

    double h[2] = {x[1] - x[0], x[2] - x[1]};

    double ref = I(x[2]) - I(x[0]);
    const double irregular_integral = ModuleBase::Integral::simpson(3, y, h);
    EXPECT_DOUBLE_EQ(ref, irregular_integral);

    // evenly-spaced 4 points
    x[0] = 0.1;
    x[1] = 0.2;
    x[2] = 0.3;
    y[0] = f(x[0]);
    y[1] = f(x[1]);
    y[2] = f(x[2]);

    ref = I(x[2]) - I(x[0]);
    const double dx = x[1] - x[0];
    const double uniform_integral = ModuleBase::Integral::simpson(3, y, dx);
    EXPECT_DOUBLE_EQ(ref, uniform_integral);
}

class SimpsonIntegralITF : public testing::Test
{
};

TEST_F(SimpsonIntegralITF, UniformGridOdd)
{
    double start = -1.0;
    double end = 1.0;
    const int ngrid_max = 10000;
    std::vector<double> f(ngrid_max);
    for (int ngrid = 5; ngrid <= ngrid_max; ngrid += 2)
    {
        const double dx = (end - start) / (ngrid - 1);
        for (int i = 0; i < ngrid; ++i)
        {
            double x = start + i * dx;
            f[i] = 1.0 / (1.0 + x * x);
        }
        double tol = (end - start) * std::pow(dx, 4) * 24 / 180;
        const double expected = std::atan(end) - std::atan(start);
        const double actual = ModuleBase::Integral::simpson(ngrid, f.data(), dx);
        const double tolerance = std::max(tol, doublethreshold);
        EXPECT_NEAR(expected, actual, tolerance);
    }
}

TEST_F(SimpsonIntegralITF, SinGridOdd)
{
    double start = -1.0;
    double end = 1.0;
    const int ngrid_max = 10000;
    std::vector<double> xv(ngrid_max);
    std::vector<double> h(ngrid_max);
    std::vector<double> f(ngrid_max);
    for (int ngrid = 5; ngrid <= ngrid_max; ngrid += 2)
    {
        sinspace(start, end, ngrid, xv.data(), h.data());
        for (int i = 0; i < ngrid; ++i)
        {
            f[i] = 1.0 / (1.0 + xv[i] * xv[i]);
        }

        // crude estimate for irregular-grid error bound
        double dx = h[ngrid / 2];
        double tol = (end - start) * std::pow(dx, 4);
        const double expected = std::atan(end) - std::atan(start);
        const double actual = ModuleBase::Integral::simpson(ngrid, f.data(), h.data());
        const double tolerance = std::max(tol, doublethreshold);
        EXPECT_NEAR(expected, actual, tolerance);
    }
}

class SimpsonIntegralExp : public testing::Test
{
};

TEST_F(SimpsonIntegralExp, UniformGridOdd)
{
    double start = 0.0;
    double end = 1.0;
    const int ngrid_max = 10000;
    std::vector<double> f(ngrid_max);
    for (int ngrid = 5; ngrid <= ngrid_max; ngrid += 2)
    {
        const double dx = (end - start) / (ngrid - 1);
        for (int i = 0; i < ngrid; ++i)
        {
            double x = start + i * dx;
            f[i] = std::exp(x);
        }
        double tol = (end - start) * std::exp(1) * std::pow(dx, 4) / 180;
        const double expected = std::exp(end) - std::exp(start);
        const double actual = ModuleBase::Integral::simpson(ngrid, f.data(), dx);
        const double tolerance = std::max(tol, doublethreshold);
        EXPECT_NEAR(expected, actual, tolerance);
    }
}

TEST_F(SimpsonIntegralExp, SinGridOdd)
{
    double start = 0.0;
    double end = 1.0;
    const int ngrid_max = 10000;
    std::vector<double> xv(ngrid_max);
    std::vector<double> h(ngrid_max);
    std::vector<double> f(ngrid_max);
    // skip ngrid = 3 since the errors exceeds the threshold
    for (int ngrid = 5; ngrid <= ngrid_max; ngrid += 2)
    {
        sinspace(start, end, ngrid, xv.data(), h.data());
        for (int i = 0; i < ngrid; ++i)
        {
            f[i] = std::exp(xv[i]);
        }

        double dx = h[ngrid / 2];
        double tol = (end - start) * std::pow(dx, 4);
        const double expected = std::exp(end) - std::exp(start);
        const double actual = ModuleBase::Integral::simpson(ngrid, f.data(), h.data());
        const double tolerance = std::max(tol, doublethreshold);
        EXPECT_NEAR(expected, actual, tolerance);
    }
}

} // namespace
