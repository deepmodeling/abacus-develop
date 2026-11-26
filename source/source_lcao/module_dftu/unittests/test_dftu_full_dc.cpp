#include "../dftu_full_dc.h"
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
using Complex = std::complex<double>;
using DFTU_LCAO::evaluate_full_dc;

std::vector<Complex> example()
{
    return {0.8, Complex(0.1, 0.2), Complex(0.1, -0.2), 0.3};
}

double contraction(const std::vector<Complex>& v, const std::vector<Complex>& n, const int size)
{
    Complex sum = 0.0;
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            sum += v[i * size + j] * n[j * size + i];
        }
    }
    EXPECT_NEAR(sum.imag(), 0.0, 1e-12);
    return sum.real();
}

// Independent two-orbital oracle in the spin-summed variable D = 2*n.
double spin_sum_energy(const std::vector<Complex>& density, const double u, const double alpha)
{
    const double d0 = density[0].real();
    const double d1 = density[3].real();
    const double trace = d0 + d1;
    const double square_trace = d0 * d0 + d1 * d1 + 2.0 * std::norm(density[1]);
    return 0.5 * u * (alpha * trace + (1.0 - alpha) * trace * trace / 4.0 - square_trace / 2.0);
}
} // namespace

TEST(DftuFullDc, AnalyticComplexShell)
{
    const std::vector<Complex> n = example();
    std::vector<Complex> v;
    double fll_correction = 0.0;
    const double fll = evaluate_full_dc(n, 2, 2.0, 1.0, v, fll_correction);
    EXPECT_NEAR(fll, 0.27, 1e-14);
    EXPECT_NEAR(fll_correction, 0.83, 1e-14);
    EXPECT_NEAR(std::abs(v[0] + 0.6), 0.0, 1e-14);
    EXPECT_NEAR(std::abs(v[3] - 0.4), 0.0, 1e-14);
    EXPECT_NEAR(std::abs(v[1] - Complex(-0.2, -0.4)), 0.0, 1e-14);
    double amf_correction = 0.0;
    const double amf = evaluate_full_dc(n, 2, 2.0, 0.0, v, amf_correction);
    EXPECT_NEAR(amf, -0.225, 1e-14);
    EXPECT_NEAR(amf_correction, 0.225, 1e-14);
    EXPECT_NEAR(std::abs(v[0] + 0.5), 0.0, 1e-14);
    EXPECT_NEAR(std::abs(v[3] - 0.5), 0.0, 1e-14);
    double mixed_correction = 0.0;
    const double mixed = evaluate_full_dc(n, 2, 2.0, 0.5, v, mixed_correction);
    EXPECT_NEAR(mixed, 0.0225, 1e-14);
    EXPECT_NEAR(mixed_correction, 0.5275, 1e-14);
}

TEST(DftuFullDc, EnergyLedgerUsesTraceContraction)
{
    const std::vector<Complex> n = example();
    for (double alpha : {0.0, 0.2, 0.5, 1.0})
    {
        std::vector<Complex> v;
        double energy_correction = 0.0;
        const double energy = evaluate_full_dc(n, 2, 3.7, alpha, v, energy_correction);
        const double band = contraction(v, n, 2);
        EXPECT_NEAR(energy - band, energy_correction, 1e-13);
    }
}

TEST(DftuFullDc, RealAndComplexAgree)
{
    const std::vector<double> n = {0.8, 0.1, 0.1, 0.3};
    const std::vector<Complex> nc(n.begin(), n.end());
    for (double alpha : {0.0, 0.5, 1.0})
    {
        std::vector<double> vr;
        std::vector<Complex> vc;
        double er_correction = 0.0;
        const double er = evaluate_full_dc(n, 2, 2.0, alpha, vr, er_correction);
        double ec_correction = 0.0;
        const double ec = evaluate_full_dc(nc, 2, 2.0, alpha, vc, ec_correction);
        EXPECT_DOUBLE_EQ(er, ec);
        EXPECT_DOUBLE_EQ(er_correction, ec_correction);
        for (int i = 0; i < 4; ++i)
        {
            EXPECT_DOUBLE_EQ(vr[i], vc[i].real());
            EXPECT_DOUBLE_EQ(vc[i].imag(), 0.0);
        }
    }
}

TEST(DftuFullDc, UniformOccupationHasZeroAmf)
{
    for (int m : {1, 3, 5, 7})
    {
        std::vector<double> n(m * m, 0.0);
        for (int i = 0; i < m; ++i)
        {
            n[i * m + i] = 0.37;
        }
        std::vector<double> v;
        double e_correction = 0.0;
        const double e = evaluate_full_dc(n, m, 4.0, 0.0, v, e_correction);
        EXPECT_NEAR(e, 0.0, 1e-14);
        EXPECT_NEAR(e_correction, 0.0, 1e-14);
        for (double value : v)
        {
            EXPECT_NEAR(value, 0.0, 1e-14);
        }
    }
}

TEST(DftuFullDc, EndpointsAndContinuityWithoutBranches)
{
    const std::vector<Complex> n = example();
    std::vector<Complex> vf;
    std::vector<Complex> va;
    double ef_correction = 0.0;
    const double ef = evaluate_full_dc(n, 2, 2.0, 1.0, vf, ef_correction);
    double ea_correction = 0.0;
    const double ea = evaluate_full_dc(n, 2, 2.0, 0.0, va, ea_correction);
    for (double alpha : {0.0, 1e-9, 0.25, 0.5, 0.75, 1.0 - 1e-9, 1.0})
    {
        std::vector<Complex> v;
        double e_correction = 0.0;
        const double e = evaluate_full_dc(n, 2, 2.0, alpha, v, e_correction);
        EXPECT_NEAR(e, alpha * ef + (1.0 - alpha) * ea, 1e-14);
        EXPECT_NEAR(e_correction, alpha * ef_correction + (1.0 - alpha) * ea_correction, 1e-14);
        for (int i = 0; i < 4; ++i)
        {
            EXPECT_NEAR(std::abs(v[i] - alpha * vf[i] - (1.0 - alpha) * va[i]), 0.0, 1e-14);
        }
    }
}

TEST(DftuFullDc, FiniteDifferencesIncludingImaginaryOffDiagonal)
{
    const std::vector<Complex> n = example();
    const std::vector<std::vector<Complex>> directions = {
        {1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 1.0, 0.0},
        {0.0, Complex(0.0, 1.0), Complex(0.0, -1.0), 0.0}};
    const double step = 1e-6;
    for (double alpha : {0.0, 0.5, 1.0})
    {
        std::vector<Complex> v;
        double band_correction = 0.0;
        evaluate_full_dc(n, 2, 2.0, alpha, v, band_correction);
        for (const std::vector<Complex>& d : directions)
        {
            std::vector<Complex> plus = n;
            std::vector<Complex> minus = n;
            for (int i = 0; i < 4; ++i)
            {
                plus[i] += step * d[i];
                minus[i] -= step * d[i];
            }
            std::vector<Complex> unused;
            double ep_correction = 0.0;
            const double ep = evaluate_full_dc(plus, 2, 2.0, alpha, unused, ep_correction);
            double em_correction = 0.0;
            const double em = evaluate_full_dc(minus, 2, 2.0, alpha, unused, em_correction);
            EXPECT_NEAR((ep - em) / (2.0 * step), contraction(v, d, 2), 2e-10);
        }
    }
}

TEST(DftuFullDc, OrbitalPhaseCovariance)
{
    const std::vector<Complex> n = example();
    const Complex phase = std::polar(1.0, 0.73);
    std::vector<Complex> rotated = n;
    rotated[1] *= phase;
    rotated[2] *= std::conj(phase);
    for (double alpha : {0.0, 0.5, 1.0})
    {
        std::vector<Complex> v;
        std::vector<Complex> vr;
        double e_correction = 0.0;
        const double e = evaluate_full_dc(n, 2, 2.0, alpha, v, e_correction);
        double er_correction = 0.0;
        const double er = evaluate_full_dc(rotated, 2, 2.0, alpha, vr, er_correction);
        EXPECT_NEAR(e, er, 1e-14);
        EXPECT_NEAR(e_correction, er_correction, 1e-14);
        EXPECT_NEAR(std::abs(vr[1] - phase * v[1]), 0.0, 1e-14);
    }
}

TEST(DftuFullDc, NonmagneticSpinSumHasAnalyticEnergyAndUndoubledPotential)
{
    const std::vector<Complex> spin_sum = {1.6, Complex(0.2, 0.4), Complex(0.2, -0.4), 0.6};
    std::vector<Complex> per_spin = spin_sum;
    for (Complex& value : per_spin)
    {
        value *= 0.5;
    }
    for (double alpha : {0.0, 0.35, 1.0})
    {
        SCOPED_TRACE(alpha);
        std::vector<Complex> v;
        double e_correction = 0.0;
        const double e = evaluate_full_dc(per_spin, 2, 2.0, alpha, v, e_correction);
        // Closed-form values for N_sigma=1.1 and Tr(n_sigma^2)=0.83.
        EXPECT_NEAR(2.0 * e, -0.45 + 0.99 * alpha, 1e-14);
        EXPECT_NEAR(2.0 * e_correction, 0.45 + 1.21 * alpha, 1e-14);
        EXPECT_NEAR(v[0].real(), -0.5 - 0.1 * alpha, 1e-14);
        EXPECT_NEAR(v[3].real(), 0.5 - 0.1 * alpha, 1e-14);
        EXPECT_NEAR(std::abs(v[1] - Complex(-0.2, -0.4)), 0.0, 1e-14);
        const double band = contraction(v, spin_sum, 2);
        EXPECT_NEAR(2.0 * e - band, 2.0 * e_correction, 1e-14);
    }
}

TEST(DftuFullDc, SpinSumChainRuleAgainstIndependentFunctional)
{
    const std::vector<Complex> spin_sum = {1.6, Complex(0.2, 0.4), Complex(0.2, -0.4), 0.6};
    const std::vector<std::vector<Complex>> directions = {
        {1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 1.0, 0.0},
        {0.0, Complex(0.0, 1.0), Complex(0.0, -1.0), 0.0}};
    std::vector<Complex> per_spin = spin_sum;
    for (Complex& value : per_spin)
    {
        value *= 0.5;
    }
    const double step = 1e-6;
    for (double alpha : {0.0, 0.35, 1.0})
    {
        SCOPED_TRACE(alpha);
        std::vector<Complex> v;
        double band_correction = 0.0;
        evaluate_full_dc(per_spin, 2, 2.0, alpha, v, band_correction);
        for (const std::vector<Complex>& direction : directions)
        {
            std::vector<Complex> plus = spin_sum;
            std::vector<Complex> minus = spin_sum;
            for (int i = 0; i < 4; ++i)
            {
                plus[i] += step * direction[i];
                minus[i] -= step * direction[i];
            }
            const double ep = spin_sum_energy(plus, 2.0, alpha);
            const double em = spin_sum_energy(minus, 2.0, alpha);
            const double derivative = (ep - em) / (2.0 * step);
            const double predicted = contraction(v, direction, 2);
            EXPECT_NEAR(derivative, predicted, 3e-10);
        }
    }
}

TEST(DftuFullDc, PolarizedShellsUseSeparateSpinAverages)
{
    const std::vector<Complex> up = example();
    const std::vector<Complex> down = {0.2, Complex(0.03, -0.04), Complex(0.03, 0.04), 0.1};
    for (double alpha : {0.0, 0.35, 1.0})
    {
        SCOPED_TRACE(alpha);
        std::vector<Complex> vu;
        std::vector<Complex> vd;
        double eu_correction = 0.0;
        const double eu = evaluate_full_dc(up, 2, 2.0, alpha, vu, eu_correction);
        double ed_correction = 0.0;
        const double ed = evaluate_full_dc(down, 2, 2.0, alpha, vd, ed_correction);
        // N_up/M=0.55, N_down/M=0.15, not their common spin average 0.35.
        EXPECT_NEAR(eu + ed, -0.235 + 0.75 * alpha, 1e-14);
        EXPECT_NEAR(eu_correction + ed_correction, 0.235 + 0.65 * alpha, 1e-14);
        EXPECT_NEAR(vu[0].real(), -0.5 - 0.1 * alpha, 1e-14);
        EXPECT_NEAR(vd[0].real(), -0.1 + 0.7 * alpha, 1e-14);
        EXPECT_NEAR(vd[3].real(), 0.1 + 0.7 * alpha, 1e-14);
        EXPECT_NEAR(std::abs(vd[1] - Complex(-0.06, 0.08)), 0.0, 1e-14);
    }
}

TEST(DftuFullDc, ZeroUAndNoPersistentState)
{
    const std::vector<Complex> n = example();
    std::vector<Complex> v;
    double original_correction = 0.0;
    const double original = evaluate_full_dc(n, 2, 2.0, 0.5, v, original_correction);
    const std::vector<Complex> original_v = v;
    double zero_correction = 0.0;
    const double zero = evaluate_full_dc(n, 2, 0.0, 0.5, v, zero_correction);
    EXPECT_DOUBLE_EQ(zero, 0.0);
    EXPECT_DOUBLE_EQ(zero_correction, 0.0);
    for (Complex value : v)
    {
        EXPECT_DOUBLE_EQ(std::abs(value), 0.0);
    }
    double again_correction = 0.0;
    const double again = evaluate_full_dc(n, 2, 2.0, 0.5, v, again_correction);
    EXPECT_EQ(original_v, v);
    EXPECT_DOUBLE_EQ(original, again);
    EXPECT_DOUBLE_EQ(original_correction, again_correction);
}

TEST(DftuFullDc, InvalidInputsDoNotOverwriteOutputs)
{
    std::vector<Complex> n = example();
    std::vector<Complex> v = {42.0};
    double band_correction = 17.0;
    EXPECT_THROW(evaluate_full_dc(n, 0, 2.0, 0.5, v, band_correction), std::invalid_argument);
    EXPECT_THROW(evaluate_full_dc(n, 3, 2.0, 0.5, v, band_correction), std::invalid_argument);
    for (double alpha : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()})
    {
        EXPECT_THROW(evaluate_full_dc(n, 2, 2.0, alpha, v, band_correction), std::invalid_argument);
    }
    n[1] = Complex(0.1, 0.3);
    EXPECT_THROW(evaluate_full_dc(n, 2, 2.0, 0.5, v, band_correction), std::invalid_argument);
    n = example();
    n[0] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(evaluate_full_dc(n, 2, 2.0, 0.5, v, band_correction), std::invalid_argument);
    ASSERT_EQ(v.size(), 1U);
    EXPECT_EQ(v[0], Complex(42.0, 0.0));
    EXPECT_DOUBLE_EQ(band_correction, 17.0);
}
