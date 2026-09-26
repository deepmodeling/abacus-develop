#include "../xc_rvv10.h"

#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>

namespace
{
using namespace Rvv10;

struct SampledChannels
{
    std::array<double, channel_count> theta;
    std::array<double, channel_count> dtheta_dn;
    std::array<double, channel_count> dtheta_dsigma;
};

SampledChannels sample_channels(const SplineBasis& basis, const LocalValues& local)
{
    SampledChannels result;
    for (int a = 0; a < channel_count; ++a)
    {
        const Channel value = basis.channel(local, a);
        result.theta[a] = value.theta;
        result.dtheta_dn[a] = value.dtheta_dn;
        result.dtheta_dsigma[a] = value.dtheta_dsigma;
    }
    return result;
}

void sample_basis(const SplineBasis& basis,
                  double q,
                  std::array<double, channel_count>& values,
                  std::array<double, channel_count>& derivatives)
{
    // Unit weight and dq/dn expose p(q) and p'(q) through the production API.
    const SampledChannels result = sample_channels(basis, {q, 1.0, 0.0, 1.0, 0.0});
    values = result.theta;
    derivatives = result.dtheta_dn;
}

TEST(Rvv10Local, IndependentHartreeConvertedReference)
{
    // 80-digit reference from the published Hartree formulas, converted to Ry.
    // q = omega/kappa is unchanged; D_Ry = D_Ha/sqrt(8).
    const LocalModel model(6.3, 0.0093);
    const LocalValues v = model.evaluate(0.01, 1.e-6);
    EXPECT_NEAR(v.q0, 0.02592385049709885768591636917636525, 2.e-16);
    EXPECT_NEAR(v.dq0_dn, 0.8639844620232583766468373760040083, 2.e-14);
    EXPECT_NEAR(v.dq0_dsigma, 0.5755515201463192674655772561435343, 2.e-14);
    EXPECT_NEAR(v.weight, 0.0001593788101414985286692409912164584, 2.e-18);
    EXPECT_NEAR(v.dweight_dn, 0.01195341076061238965019307434123438, 2.e-16);
    EXPECT_NEAR(model.beta(), 0.009009695931388873409570041483960262, 2.e-16);
}

TEST(Rvv10Local, DensityAndSigmaDerivatives)
{
    const LocalModel model(6.3, 0.0093);
    for (double n: {0.001, 0.01, 0.1})
    {
        for (double sigma: {1.e-8, 1.e-6, 1.e-4})
        {
            const LocalValues v = model.evaluate(n, sigma);
            const double hn = n * 1.e-5;
            const double hs = sigma * 1.e-4;
            const LocalValues np = model.evaluate(n + hn, sigma);
            const LocalValues nm = model.evaluate(n - hn, sigma);
            EXPECT_NEAR(v.dq0_dn, (np.q0 - nm.q0) / (2 * hn), 2.e-7 * (1 + std::abs(v.dq0_dn)));
            EXPECT_NEAR(v.dweight_dn, (np.weight - nm.weight) / (2 * hn), 1.e-10);
            EXPECT_NEAR(v.dq0_dsigma,
                        (model.evaluate(n, sigma + hs).q0 - model.evaluate(n, sigma - hs).q0) / (2 * hs),
                        2.e-6 * (1 + std::abs(v.dq0_dsigma)));
        }
    }
}

TEST(Rvv10Local, ZeroAndSmallGradientAreContinuous)
{
    const LocalModel model(6.3, 0.0093);
    const LocalValues zero = model.evaluate(0.001, 0.0);
    const LocalValues tiny = model.evaluate(0.001, 1.e-14);
    EXPECT_DOUBLE_EQ(zero.dq0_dsigma, 0.0);
    EXPECT_GT(tiny.dq0_dsigma, 0.0); // Do not copy QE's sigma <= 1.e-12 derivative mask.
    EXPECT_NEAR(tiny.q0, zero.q0, 1.e-15);
    const LocalValues twice = model.evaluate(0.001, 2.e-14);
    EXPECT_NEAR(twice.dq0_dsigma / tiny.dq0_dsigma, 2.0, 1.e-10);
}

TEST(Rvv10Local, DensityCutoffAndLowerClip)
{
    const LocalModel model(6.3, 0.0093);
    for (double n: {-0.1, 0.0, density_cutoff})
    {
        const LocalValues v = model.evaluate(n, 1.e-3);
        EXPECT_DOUBLE_EQ(v.q0, q_cut);
        EXPECT_DOUBLE_EQ(v.weight, 0.0);
        EXPECT_DOUBLE_EQ(v.dweight_dn, 0.0);
        EXPECT_DOUBLE_EQ(v.dq0_dn, 0.0);
        EXPECT_DOUBLE_EQ(v.dq0_dsigma, 0.0);
    }
    const LocalValues clipped = model.evaluate(1.e-10, 0.0);
    EXPECT_DOUBLE_EQ(clipped.q0, q_min);
    // QE clips q0 for interpolation but retains the analytic dq0/dn used by
    // the potential; only the explicit density cutoff zeros derivatives.
    EXPECT_GT(clipped.dq0_dn, 0.0);
    EXPECT_DOUBLE_EQ(clipped.dq0_dsigma, 0.0);
    EXPECT_GT(clipped.weight, 0.0);
    // beta is a separate local term, not masked together with theta.
    EXPECT_GT(model.beta(), 0.0);
}

TEST(Rvv10Local, SaturationDoesNotOverflow)
{
    const LocalModel model(6.3, 0.0093);
    for (double sigma: {1.0, 1.e100, 1.e300})
    {
        const LocalValues v = model.evaluate(0.01, sigma);
        EXPECT_DOUBLE_EQ(v.q0, q_cut);
        EXPECT_DOUBLE_EQ(v.dq0_dn, 0.0);
        EXPECT_DOUBLE_EQ(v.dq0_dsigma, 0.0);
        EXPECT_TRUE(std::isfinite(v.weight));
    }
}

TEST(Rvv10Local, ExplicitParametersAndInvalidInput)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(LocalModel(0.0, 0.0093), std::invalid_argument);
    EXPECT_THROW(LocalModel(-1.0, 0.0093), std::invalid_argument);
    EXPECT_THROW(LocalModel(6.3, -0.1), std::invalid_argument);
    EXPECT_THROW(LocalModel(nan, 0.0093), std::invalid_argument);
    EXPECT_THROW(LocalModel(6.3, inf), std::invalid_argument);
    const LocalModel model(6.3, 0.0093);
    EXPECT_THROW(model.evaluate(nan, 0.0), std::invalid_argument);
    EXPECT_THROW(model.evaluate(0.01, inf), std::invalid_argument);
    EXPECT_THROW(model.evaluate(0.01, -1.e-6), std::invalid_argument);
    const LocalModel zero_c(6.3, 0.0);
    EXPECT_DOUBLE_EQ(zero_c.evaluate(0.01, 1.e300).q0, zero_c.evaluate(0.01, 0.0).q0);
    const LocalModel other_b(12.6, 0.0093);
    EXPECT_NEAR(model.beta() / other_b.beta(), std::sqrt(8.0), 1.e-14);
    EXPECT_NEAR(model.evaluate(0.01, 0.0).weight / other_b.evaluate(0.01, 0.0).weight, std::sqrt(8.0), 1.e-14);
}

TEST(Rvv10Kernel, AnalyticLimitsAndSymmetry)
{
    EXPECT_DOUBLE_EQ(real_kernel(0.01, 0.2, 0.0), -12.0);
    EXPECT_DOUBLE_EQ(real_kernel(0.01, 0.2, 3.0), real_kernel(0.2, 0.01, 3.0));
    const double q = 0.02;
    const double r = 4.0;
    EXPECT_NEAR(real_kernel(q, q, r), -12.0 / std::pow(1 + q * r * r, 3), 1.e-14);
    const double far = 1.e5;
    EXPECT_NEAR(real_kernel(0.01, 0.2, far) * std::pow(far, 6), -24.0 / (0.01 * 0.2 * (0.01 + 0.2)), 0.002);
    EXPECT_TRUE(std::isfinite(real_kernel(0.01, 0.2, 1.e300)));
    EXPECT_THROW(real_kernel(0.0, 0.2, 1.0), std::invalid_argument);
    EXPECT_THROW(real_kernel(0.01, 0.2, -1.0), std::invalid_argument);
}

TEST(Rvv10Kernel, AvoidPrematureDistanceSquareOverflow)
{
    // q*r^2 = 1.e10 is finite even though r^2 alone overflows double.
    EXPECT_NEAR(real_kernel(1.e-300, 1.e-300, 1.e155), -1.19999999964000000007e-29, 1.e-42);
}

TEST(Rvv10Spline, CardinalValuesAtEveryKnot)
{
    const SplineBasis basis;
    std::array<double, channel_count> p;
    std::array<double, channel_count> dp;
    for (int i = 0; i < channel_count; ++i)
    {
        sample_basis(basis, q_mesh()[i], p, dp);
        for (int j = 0; j < channel_count; ++j)
        {
            EXPECT_NEAR(p[j], i == j ? 1.0 : 0.0, 2.e-14);
        }
    }
}

TEST(Rvv10Spline, PartitionOfUnityAndDerivative)
{
    const SplineBasis basis;
    std::array<double, channel_count> p;
    std::array<double, channel_count> dp;
    std::array<double, channel_count> pp;
    std::array<double, channel_count> dpp;
    std::array<double, channel_count> pm;
    std::array<double, channel_count> dpm;
    for (int i = 0; i < channel_count - 1; ++i)
    {
        const double q = 0.5 * (q_mesh()[i] + q_mesh()[i + 1]);
        const double h = (q_mesh()[i + 1] - q_mesh()[i]) * 1.e-5;
        sample_basis(basis, q, p, dp);
        sample_basis(basis, q + h, pp, dpp);
        sample_basis(basis, q - h, pm, dpm);
        double sum_p = 0;
        double sum_dp = 0;
        for (int j = 0; j < channel_count; ++j)
        {
            sum_p += p[j];
            sum_dp += dp[j];
            EXPECT_NEAR(dp[j], (pp[j] - pm[j]) / (2 * h), 2.e-6 * (1 + std::abs(dp[j])));
        }
        EXPECT_NEAR(sum_p, 1.0, 1.e-13);
        EXPECT_NEAR(sum_dp, 0.0, 1.e-10);
    }
    EXPECT_THROW(sample_basis(basis, q_min * 0.5, p, dp), std::out_of_range);
    EXPECT_THROW(sample_basis(basis, q_cut + 0.001, p, dp), std::out_of_range);
}

TEST(Rvv10Spline, NaturalRatherThanDefaultNotAKnotBoundaries)
{
    const SplineBasis basis;
    std::array<double, channel_count> p;
    std::array<double, channel_count> d0;
    std::array<double, channel_count> d1;
    std::array<double, channel_count> d2;
    for (int endpoint: {0, channel_count - 1})
    {
        const double q = q_mesh()[endpoint];
        const double width = endpoint == 0 ? q_mesh()[1] - q : q - q_mesh()[endpoint - 1];
        const double h = (endpoint == 0 ? 1 : -1) * width * 1.e-3;
        sample_basis(basis, q, p, d0);
        sample_basis(basis, q + h, p, d1);
        sample_basis(basis, q + 2 * h, p, d2);
        for (int j = 0; j < channel_count; ++j)
        {
            // One-sided derivative of the quadratic p' is exact to roundoff.
            EXPECT_NEAR((-3 * d0[j] + 4 * d1[j] - d2[j]) / (2 * h), 0.0, 5.e-5);
        }
    }
}

TEST(Rvv10Spline, ThetaAndItsLocalChainRule)
{
    const LocalModel model(6.3, 0.0093);
    const SplineBasis basis;
    const double n = 0.01;
    const double sigma = 1.e-6;
    const double hn = 1.e-7;
    const double hs = 1.e-9;
    const SampledChannels t = sample_channels(basis, model.evaluate(n, sigma));
    const SampledChannels np = sample_channels(basis, model.evaluate(n + hn, sigma));
    const SampledChannels nm = sample_channels(basis, model.evaluate(n - hn, sigma));
    const SampledChannels sp = sample_channels(basis, model.evaluate(n, sigma + hs));
    const SampledChannels sm = sample_channels(basis, model.evaluate(n, sigma - hs));
    double sum_t = 0;
    double sum_n = 0;
    double sum_s = 0;
    for (int a = 0; a < channel_count; ++a)
    {
        sum_t += t.theta[a];
        sum_n += t.dtheta_dn[a];
        sum_s += t.dtheta_dsigma[a];
        EXPECT_NEAR(t.dtheta_dn[a], (np.theta[a] - nm.theta[a]) / (2 * hn), 1.e-9);
        EXPECT_NEAR(t.dtheta_dsigma[a], (sp.theta[a] - sm.theta[a]) / (2 * hs), 1.e-9);
    }
    const LocalValues v = model.evaluate(n, sigma);
    EXPECT_NEAR(sum_t, v.weight, 1.e-17);
    EXPECT_NEAR(sum_n, v.dweight_dn, 1.e-14);
    EXPECT_NEAR(sum_s, 0.0, 1.e-15);
    const SampledChannels vacuum = sample_channels(basis, model.evaluate(0.0, 0.0));
    for (int a = 0; a < channel_count; ++a)
    {
        EXPECT_DOUBLE_EQ(vacuum.theta[a], 0.0);
        EXPECT_DOUBLE_EQ(vacuum.dtheta_dn[a], 0.0);
        EXPECT_DOUBLE_EQ(vacuum.dtheta_dsigma[a], 0.0);
    }
}
TEST(Rvv10Spline, SaturatedLargeWeightKeepsZeroGradientDerivativeFinite)
{
    const LocalModel model(1.e-52, 0.0093);
    const SplineBasis basis;
    const LocalValues v = model.evaluate(1.e308, 0.0);
    ASSERT_TRUE(std::isfinite(v.weight));
    ASSERT_DOUBLE_EQ(v.q0, q_cut);
    const SampledChannels t = sample_channels(basis, v);
    for (int a = 0; a < channel_count; ++a)
    {
        EXPECT_TRUE(std::isfinite(t.theta[a]));
        EXPECT_TRUE(std::isfinite(t.dtheta_dn[a]));
        EXPECT_DOUBLE_EQ(t.dtheta_dsigma[a], 0.0);
    }
}

TEST(Rvv10Spline, InvalidChannelAndNonfiniteArgument)
{
    const LocalModel model(6.3, 0.0093);
    const SplineBasis basis;
    LocalValues local = model.evaluate(0.01, 1.e-6);
    EXPECT_THROW(basis.channel(local, -1), std::out_of_range);
    EXPECT_THROW(basis.channel(local, channel_count), std::out_of_range);
    local.q0 = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(basis.channel(local, 0), std::out_of_range);
    local.q0 = std::numeric_limits<double>::infinity();
    EXPECT_THROW(basis.channel(local, 0), std::out_of_range);
}

TEST(Rvv10Table, FiniteRadiusZeroModeAndFourierDiagonal)
{
    const KernelTable table;
    const KernelValues zero = table.evaluate(0.0);
    const double pi = std::acos(-1.0);
    for (int a: {0, 10, 19})
    {
        const double q = q_mesh()[a];
        const double x = std::sqrt(q) * 100.0;
        // Closed-form integral to R=100, independently integrated analytically.
        const double exact = -6 * pi * std::pow(q, -1.5) * (std::atan(x) - x * (1 - x * x) / std::pow(1 + x * x, 2));
        EXPECT_NEAR(zero.value[a * channel_count + a], exact, 3.e-7 * std::abs(exact));
    }
    // At the largest q the R=100 tail is small; compare two nonzero mesh modes
    // with the analytic infinite-radius transform -3*pi^2/q^(3/2)*(1+t)*exp(-t).
    for (int mode: {1, 5})
    {
        const double g = mode * 2 * pi / 100;
        const double q = 0.5;
        const double t = g / std::sqrt(q);
        const double exact = -3 * pi * pi * std::pow(q, -1.5) * (1 + t) * std::exp(-t);
        EXPECT_NEAR(table.evaluate(g).value[399], exact, 0.0001);
    }
}

TEST(Rvv10Table, SymmetryDerivativeAndRange)
{
    const KernelTable table;
    const double g = 0.173;
    const double h = 1.e-6;
    const KernelValues v = table.evaluate(g);
    const KernelValues plus = table.evaluate(g + h);
    const KernelValues minus = table.evaluate(g - h);
    for (int a = 0; a < channel_count; ++a)
    {
        for (int b = 0; b < channel_count; ++b)
        {
            const int ab = a * channel_count + b;
            const int ba = b * channel_count + a;
            EXPECT_DOUBLE_EQ(v.value[ab], v.value[ba]);
            EXPECT_DOUBLE_EQ(v.derivative[ab], v.derivative[ba]);
            EXPECT_NEAR(v.derivative[ab],
                        (plus.value[ab] - minus.value[ab]) / (2 * h),
                        1.e-7 * (1 + std::abs(v.derivative[ab])));
        }
    }
    EXPECT_THROW(table.evaluate(-0.01), std::out_of_range);
    EXPECT_THROW(table.evaluate(table.maximum_g()), std::out_of_range);
    EXPECT_THROW(table.evaluate(std::numeric_limits<double>::quiet_NaN()), std::out_of_range);
}

TEST(Rvv10Table, ValueOnlyPathPreservesAnalyticLimitsAndRangeChecks)
{
    const KernelTable table;
    const auto zero = table.values(0.0);
    const double pi = std::acos(-1.0);
    for (int a: {0, 10, 19})
    {
        const double q = q_mesh()[a];
        const double x = std::sqrt(q) * 100.0;
        // Same independent finite-radius integral as the full-output path.
        const double exact = -6 * pi * std::pow(q, -1.5) * (std::atan(x) - x * (1 - x * x) / std::pow(1 + x * x, 2));
        EXPECT_NEAR(zero[a * channel_count + a], exact, 3.e-7 * std::abs(exact));
    }
    for (int mode: {1, 5})
    {
        const double g = mode * 2 * pi / 100;
        const double t = g / std::sqrt(0.5);
        const double exact = -3 * pi * pi * std::pow(0.5, -1.5) * (1 + t) * std::exp(-t);
        EXPECT_NEAR(table.values(g)[399], exact, 0.0001);
    }
    EXPECT_THROW(table.values(-0.01), std::out_of_range);
    EXPECT_THROW(table.values(table.maximum_g()), std::out_of_range);
    EXPECT_THROW(table.values(std::numeric_limits<double>::quiet_NaN()), std::out_of_range);
    EXPECT_THROW(table.values(std::numeric_limits<double>::infinity()), std::out_of_range);
}

// Test-only five-point quadrature, not a substitute for a PW_Basis FFT test.
// Its periodic centered gradient has the adjoint -divergence exactly.
double discrete_energy(const std::array<double, 5>& density, std::array<double, 5>& potential)
{
    const LocalModel model(6.3, 0.0093);
    const SplineBasis basis;
    const double dx = 0.8;
    const double dv = dx * dx * dx;
    std::array<double, 5> grad;
    std::array<double, 5> vector_field;
    std::array<SampledChannels, 5> fields;
    for (int i = 0; i < 5; ++i)
    {
        grad[i] = (density[(i + 1) % 5] - density[(i + 4) % 5]) / (2 * dx);
        fields[i] = sample_channels(basis, model.evaluate(density[i], grad[i] * grad[i]));
    }
    double energy = 0.0;
    for (int i = 0; i < 5; ++i)
    {
        std::array<double, channel_count> u = {};
        for (int j = 0; j < 5; ++j)
        {
            const int steps = std::abs(i - j);
            const double r = dx * (steps > 2 ? 5 - steps : steps);
            for (int a = 0; a < channel_count; ++a)
            {
                for (int b = 0; b < channel_count; ++b)
                {
                    u[a] += dv * real_kernel(q_mesh()[a], q_mesh()[b], r) * fields[j].theta[b];
                }
            }
        }
        energy += dv * model.beta() * density[i];
        potential[i] = model.beta();
        vector_field[i] = 0.0;
        for (int a = 0; a < channel_count; ++a)
        {
            energy += 0.5 * dv * fields[i].theta[a] * u[a];
            potential[i] += u[a] * fields[i].dtheta_dn[a];
            vector_field[i] += 2 * grad[i] * u[a] * fields[i].dtheta_dsigma[a];
        }
    }
    for (int i = 0; i < 5; ++i)
    {
        potential[i] -= (vector_field[(i + 1) % 5] - vector_field[(i + 4) % 5]) / (2 * dx);
    }
    return energy;
}

TEST(Rvv10Variation, DiscreteEnergyDerivativeIncludesNegativeDivergence)
{
    std::array<double, 5> density = {{0.01, 0.012, 0.008, 0.015, 0.006}};
    std::array<double, 5> potential;
    std::array<double, 5> unused;
    const double energy = discrete_energy(density, potential);
    EXPECT_TRUE(std::isfinite(energy));
    const double h = 1.e-7;
    const double dv = 0.8 * 0.8 * 0.8;
    for (int i = 0; i < 5; ++i)
    {
        std::array<double, 5> plus = density;
        std::array<double, 5> minus = density;
        plus[i] += h;
        minus[i] -= h;
        const double fd = (discrete_energy(plus, unused) - discrete_energy(minus, unused)) / (2 * h * dv);
        EXPECT_NEAR(potential[i], fd, 2.e-10);
    }
}
} // namespace
