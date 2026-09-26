#include "../xc_functional.h"
#include "../xc_rvv10.h"
#include "../xc_rvv10_pw.h"
#include "source_basis/module_pw/pw_basis.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
const double pi = std::acos(-1.0);

class Rvv10PW : public testing::Test
{
  protected:
    ModulePW::PW_Basis pw{"cpu", "double"};

    void SetUp() override
    {
        const ModuleBase::Matrix3 lattice(1, 0, 0, 0, 1.25, 0, 0, 0, 1.5);
        pw.initgrids(10.0, lattice, 6, 8, 10);
        pw.initparameters(false, 100.0, 1, true);
        pw.setuptransform();
        pw.collect_local_pw();
    }

    std::vector<double> wave() const
    {
        std::vector<double> v(pw.nrxx);
        for (int x = 0; x < pw.nx; ++x)
            for (int y = 0; y < pw.ny; ++y)
                for (int z = 0; z < pw.nz; ++z)
                    v[(x * pw.ny + y) * pw.nz + z] = std::cos(2 * pi * x / pw.nx) + 0.7 * std::sin(2 * pi * y / pw.ny)
                                                     + 0.4 * std::cos(2 * pi * z / pw.nz);
        return v;
    }
};

TEST_F(Rvv10PW, RealFFTNormalizationAndSpectralDerivatives)
{
    const std::vector<double> n = wave();
    std::vector<double> back(pw.nrxx);
    std::vector<double> laplacian(pw.nrxx);
    std::vector<std::complex<double>> ng(pw.npw);
    std::vector<ModuleBase::Vector3<double>> gradient(pw.nrxx);
    pw.real2recip(n.data(), ng.data());
    pw.recip2real(ng.data(), back.data());
    XC_Functional::grad_rho(ng.data(), gradient.data(), &pw, pw.tpiba);
    XC_Functional::grad_dot(gradient.data(), laplacian.data(), &pw, pw.tpiba);
    double real_norm = 0;
    double reciprocal_norm = 0;
    for (int x = 0; x < pw.nx; ++x)
        for (int y = 0; y < pw.ny; ++y)
            for (int z = 0; z < pw.nz; ++z)
            {
                const int i = (x * pw.ny + y) * pw.nz + z;
                EXPECT_NEAR(back[i], n[i], 2.e-14);
                EXPECT_NEAR(gradient[i].x, -2 * pi / 10 * std::sin(2 * pi * x / pw.nx), 2.e-14);
                EXPECT_NEAR(gradient[i].y, 0.7 * 2 * pi / 12.5 * std::cos(2 * pi * y / pw.ny), 2.e-14);
                EXPECT_NEAR(gradient[i].z, -0.4 * 2 * pi / 15 * std::sin(2 * pi * z / pw.nz), 2.e-14);
                const double exact = -std::pow(2 * pi / 10, 2) * std::cos(2 * pi * x / pw.nx)
                                     - 0.7 * std::pow(2 * pi / 12.5, 2) * std::sin(2 * pi * y / pw.ny)
                                     - 0.4 * std::pow(2 * pi / 15, 2) * std::cos(2 * pi * z / pw.nz);
                EXPECT_NEAR(laplacian[i], exact, 3.e-14);
                real_norm += n[i] * n[i] / pw.nxyz;
            }
    for (const auto& v: ng)
        reciprocal_norm += std::norm(v);
    EXPECT_NEAR(real_norm, reciprocal_norm, 2.e-14);
}
TEST_F(Rvv10PW, ConstantDensityEnergyAndValenceContraction)
{
    const Rvv10::SerialEvaluator evaluator(6.3, 0.0093);
    std::vector<double> total(pw.nrxx, 0.01);
    std::vector<double> valence(pw.nrxx, 0.006);
    const auto value = evaluator.evaluate(pw, total, valence);
    const Rvv10::LocalModel local(6.3, 0.0093);
    const Rvv10::SplineBasis basis;
    const Rvv10::KernelTable kernel;
    std::array<Rvv10::Channel, Rvv10::channel_count> t;
    for (int a = 0; a < Rvv10::channel_count; ++a)
        t[a] = basis.channel(local.evaluate(0.01, 0.0), a);
    const auto k = kernel.evaluate(0.0);
    double e = local.beta() * 0.01;
    double v = local.beta();
    for (int a = 0; a < 20; ++a)
        for (int b = 0; b < 20; ++b)
        {
            e += 0.5 * t[a].theta * k.value[a * 20 + b] * t[b].theta;
            v += t[a].dtheta_dn * k.value[a * 20 + b] * t[b].theta;
        }
    EXPECT_NEAR(value.energy, pw.omega * e, 1.e-12);
    for (double vr: value.potential)
        EXPECT_NEAR(vr, v, 1.e-12);
    EXPECT_NEAR(value.vtxc, pw.omega * 0.006 * v, 1.e-12);
    const auto all_valence = evaluator.evaluate(pw, total, total);
    EXPECT_DOUBLE_EQ(value.energy, all_valence.energy);
    EXPECT_NEAR(all_valence.vtxc / value.vtxc, 0.01 / 0.006, 1.e-12);
}

TEST_F(Rvv10PW, NonuniformEnergyVariationWithFixedCore)
{
    const Rvv10::SerialEvaluator evaluator(6.3, 0.0093);
    const auto shape = wave();
    std::vector<double> total(pw.nrxx);
    std::vector<double> valence(pw.nrxx);
    std::vector<double> direction(pw.nrxx);
    for (int i = 0; i < pw.nrxx; ++i)
    {
        valence[i] = 0.008 + 0.002 * shape[i];
        total[i] = valence[i] + 0.003;
        direction[i] = 0.3 + shape[i];
    }
    const auto value = evaluator.evaluate(pw, total, valence);
    const double dv = pw.omega / pw.nxyz;
    const double h = 1.e-7;
    double derivative = 0;
    double contraction = 0;
    for (int i = 0; i < pw.nrxx; ++i)
    {
        derivative += value.potential[i] * direction[i] * dv;
        contraction += value.potential[i] * valence[i] * dv;
    }
    auto plus = total;
    auto minus = total;
    for (int i = 0; i < pw.nrxx; ++i)
    {
        plus[i] += h * direction[i];
        minus[i] -= h * direction[i];
    }
    const double fd
        = (evaluator.evaluate(pw, plus, valence).energy - evaluator.evaluate(pw, minus, valence).energy) / (2 * h);
    EXPECT_NEAR(derivative, fd, 2.e-7);
    EXPECT_NEAR(value.vtxc, contraction, 1.e-13);
    // Local perturbation also probes modes removed by the density PW cutoff.
    plus = total;
    minus = total;
    plus[13] += h;
    minus[13] -= h;
    const double point_fd
        = (evaluator.evaluate(pw, plus, valence).energy - evaluator.evaluate(pw, minus, valence).energy) / (2 * h * dv);
    EXPECT_NEAR(value.potential[13], point_fd, 2.e-9);
}

TEST_F(Rvv10PW, TranslationCovarianceAndRepeatability)
{
    const Rvv10::SerialEvaluator evaluator(6.3, 0.0093);
    auto n = wave();
    auto shifted = n;
    for (auto& v: n)
        v = 0.01 + 0.002 * v;
    for (int x = 0; x < pw.nx; ++x)
        for (int y = 0; y < pw.ny; ++y)
            for (int z = 0; z < pw.nz; ++z)
                shifted[(x * pw.ny + y) * pw.nz + z] = n[(((x + 1) % pw.nx) * pw.ny + y) * pw.nz + z];
    const auto a = evaluator.evaluate(pw, n, n);
    const auto b = evaluator.evaluate(pw, shifted, shifted);
    EXPECT_NEAR(a.energy, b.energy, 1.e-12);
    for (int x = 0; x < pw.nx; ++x)
        for (int y = 0; y < pw.ny; ++y)
            for (int z = 0; z < pw.nz; ++z)
                EXPECT_NEAR(b.potential[(x * pw.ny + y) * pw.nz + z],
                            a.potential[(((x + 1) % pw.nx) * pw.ny + y) * pw.nz + z],
                            1.e-12);
    const auto repeated = evaluator.evaluate(pw, n, n);
    EXPECT_DOUBLE_EQ(a.energy, repeated.energy);
}

TEST_F(Rvv10PW, RejectUnsupportedLayoutsAndMalformedDensity)
{
    const Rvv10::SerialEvaluator evaluator(6.3, 0.0093);
    const std::vector<double> n(pw.nrxx, 0.01);
    const std::vector<double> short_n(1, 0.01);
    EXPECT_THROW(evaluator.evaluate(pw, n, short_n), std::invalid_argument);
    pw.gamma_only = true;
    EXPECT_THROW(evaluator.evaluate(pw, n, n), std::invalid_argument);
    pw.gamma_only = false;
    pw.poolnproc = 2;
    EXPECT_THROW(evaluator.evaluate(pw, n, n), std::invalid_argument);
    pw.poolnproc = 1;
    auto bad = n;
    bad[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(evaluator.evaluate(pw, bad, n), std::invalid_argument);
    EXPECT_THROW(evaluator.evaluate(pw, n, bad), std::invalid_argument);
    ModulePW::PW_Basis empty("cpu", "double");
    EXPECT_THROW(evaluator.evaluate(empty, {}, {}), std::invalid_argument);
}

class Rvv10VacuumPW : public testing::Test
{
  protected:
    ModulePW::PW_Basis pw{"cpu", "double"};

    void SetUp() override
    {
        // A nonorthogonal cell and a density cutoff smaller than the FFT box
        // exercise the reciprocal metric and projection in grad_rho/grad_dot.
        const ModuleBase::Matrix3 lattice(1, 0, 0, -0.5, std::sqrt(0.75), 0, 0, 0, 2);
        pw.initgrids(10.0, lattice, 6, 8, 12);
        pw.initparameters(false, 4.0, 1, true);
        pw.setuptransform();
        pw.collect_local_pw();
    }
};

TEST_F(Rvv10VacuumPW, ProjectedGradientAndDivergenceAreDiscreteAdjoints)
{
    // Independent scalar/vector fields include modes outside the density basis.
    // The nonorthogonal metric and the same G projection must enter both sides.
    ASSERT_LT(pw.npw, pw.nxyz);
    std::vector<double> scalar(pw.nrxx);
    std::vector<double> projected(pw.nrxx);
    std::vector<double> divergence(pw.nrxx);
    std::vector<std::complex<double>> reciprocal(pw.npw);
    std::vector<ModuleBase::Vector3<double>> gradient(pw.nrxx);
    std::vector<ModuleBase::Vector3<double>> field(pw.nrxx);
    for (int x = 0; x < pw.nx; ++x)
        for (int y = 0; y < pw.ny; ++y)
            for (int z = 0; z < pw.nz; ++z)
            {
                const int i = (x * pw.ny + y) * pw.nz + z;
                const double ax = 2 * pi * x / pw.nx;
                const double ay = 2 * pi * y / pw.ny;
                const double az = 2 * pi * z / pw.nz;
                const double high_a = std::cos(3 * ax + 3 * ay);
                const double high_b = std::sin(4 * ay + 5 * az);
                scalar[i] = std::cos(ax) + 0.3 * std::sin(ay) + 0.2 * std::cos(az)
                            + 0.11 * high_a + 0.09 * high_b;
                field[i].x = std::sin(ax) + 0.2 * std::cos(ay + az) + 0.07 * high_a;
                field[i].y = 0.8 * std::cos(ay) + 0.3 * std::sin(ax + az) + 0.09 * high_b;
                field[i].z = -0.4 * std::sin(az) + 0.4 * std::cos(ax - ay) + 0.05 * high_a;
            }
    pw.real2recip(scalar.data(), reciprocal.data());
    pw.recip2real(reciprocal.data(), projected.data());
    XC_Functional::grad_rho(reciprocal.data(), gradient.data(), &pw, pw.tpiba);
    XC_Functional::grad_dot(field.data(), divergence.data(), &pw, pw.tpiba);
    const double dv = pw.omega / pw.nxyz;
    double forward = 0;
    double reverse = 0;
    double removed_norm = 0;
    double scale = 0;
    for (int i = 0; i < pw.nrxx; ++i)
    {
        const double left = dv * (gradient[i].x * field[i].x + gradient[i].y * field[i].y
                                  + gradient[i].z * field[i].z);
        const double right = -dv * scalar[i] * divergence[i];
        forward += left;
        reverse += right;
        scale += std::abs(left) + std::abs(right);
        removed_norm += dv * std::pow(scalar[i] - projected[i], 2);
    }
    EXPECT_GT(removed_norm, 1.e-3); // Not accidentally a full-grid identity.
    EXPECT_GT(std::abs(forward), 1.0); // A zero contraction cannot test the sign.
    EXPECT_NEAR(forward, reverse, 2.e-13 * (1 + scale));
}

TEST_F(Rvv10VacuumPW, SmallGradientPotentialIsDiscreteEnergyDerivative)
{
    // Regression for suppressing q0 derivatives when |grad n|^2 <= 1e-12:
    // small absolute gradients do not imply a small |grad n|/n in vacuum.
    // The two scales put the same smooth field below and across that threshold.
    const Rvv10::SerialEvaluator evaluator(6.3, 0.0093);
    const double scales[] = {1.e-10, 1.e-5};
    const double relative_steps[] = {1.e-3, 3.e-4, 1.e-4};
    const double dv = pw.omega / pw.nxyz;
    ASSERT_LT(pw.npw, pw.nxyz);
    for (double scale: scales)
    {
        SCOPED_TRACE(testing::Message() << "density scale = " << scale);
        std::vector<double> total(pw.nrxx);
        std::vector<double> valence(pw.nrxx);
        double max_sigma = 0;
        for (int x = 0; x < pw.nx; ++x)
            for (int y = 0; y < pw.ny; ++y)
                for (int z = 0; z < pw.nz; ++z)
                {
                    const int i = (x * pw.ny + y) * pw.nz + z;
                    const double ax = 2 * pi * x / pw.nx;
                    const double ay = 2 * pi * y / pw.ny;
                    const double az = 2 * pi * z / pw.nz;
                    total[i] = scale * (1 + 0.2 * std::cos(ax) + 0.1 * std::sin(ay) + 0.15 * std::cos(az));
                    valence[i] = total[i] - 0.25 * scale; // Core held fixed below.
                    // Cartesian reciprocal vectors of the stated hexagonal cell.
                    const double gx = -0.2 * scale * (2 * pi / 10) * std::sin(ax);
                    const double gy
                        = gx / std::sqrt(3.0) + 0.1 * scale * (4 * pi / (10 * std::sqrt(3.0))) * std::cos(ay);
                    const double gz = -0.15 * scale * (2 * pi / 20) * std::sin(az);
                    const double sigma = gx * gx + gy * gy + gz * gz;
                    if (sigma > max_sigma)
                        max_sigma = sigma;
                    ASSERT_GT(valence[i], 10 * Rvv10::density_cutoff);
                }
        if (scale == scales[0])
            EXPECT_LT(max_sigma, 1.e-12);
        else
            EXPECT_GT(max_sigma, 1.e-12);

        const auto value = evaluator.evaluate(pw, total, valence);
        ASSERT_TRUE(std::isfinite(value.energy));
        double valence_contraction = 0;
        for (int i = 0; i < pw.nrxx; ++i)
        {
            ASSERT_TRUE(std::isfinite(value.potential[i]));
            valence_contraction += dv * valence[i] * value.potential[i];
        }
        EXPECT_NEAR(value.vtxc, valence_contraction, 1.e-14);

        // At scale 1e-5 these points lie below/above the old gradient threshold,
        // respectively. Both remain far above the density cutoff at every h.
        const int points[] = {(0 * pw.ny + 2) * pw.nz + 3, (1 * pw.ny + 2) * pw.nz + 3};
        for (int point: points)
        {
            SCOPED_TRACE(testing::Message() << "grid index = " << point);
            for (double relative_step: relative_steps)
            {
                SCOPED_TRACE(testing::Message() << "relative step = " << relative_step);
                const double h = relative_step * total[point];
                auto total_plus = total;
                auto total_minus = total;
                auto valence_plus = valence;
                auto valence_minus = valence;
                total_plus[point] += h;
                total_minus[point] -= h;
                valence_plus[point] += h;
                valence_minus[point] -= h;
                const double fd = (evaluator.evaluate(pw, total_plus, valence_plus).energy
                                   - evaluator.evaluate(pw, total_minus, valence_minus).energy)
                                  / (2 * h * dv);
                // The energy finite difference is independent of the potential
                // assembly, including its projected divergence and beta term.
                EXPECT_NEAR(value.potential[point], fd, 2.e-6);
            }
        }
    }
}

TEST_F(Rvv10VacuumPW, InactiveDensityRetainsFiniteBetaEnergyAndPotential)
{
    // The nonlocal channels vanish below the density cutoff, but beta*n and
    // its constant derivative must remain, even for small negative FFT noise.
    const Rvv10::SerialEvaluator evaluator(6.3, 0.0093);
    const double density[] = {-1.e-13, 0, 5.e-13, 1.e-12};
    const double beta = std::pow(3.0, 0.75) * std::pow(6.3, -1.5) / 16;
    const double dv = pw.omega / pw.nxyz;
    std::vector<double> total(pw.nrxx);
    std::vector<double> valence(pw.nrxx);
    double expected_energy = 0;
    double expected_vtxc = 0;
    for (int i = 0; i < pw.nrxx; ++i)
    {
        total[i] = density[i % 4];
        valence[i] = 0.5 * total[i];
        expected_energy += dv * beta * total[i];
        expected_vtxc += dv * beta * valence[i];
    }
    const auto value = evaluator.evaluate(pw, total, valence);
    EXPECT_NEAR(value.energy, expected_energy, 1.e-23);
    EXPECT_NEAR(value.vtxc, expected_vtxc, 1.e-23);
    for (double potential: value.potential)
        EXPECT_NEAR(potential, beta, 1.e-15);
}

TEST(Rvv10Vacuum, MaterialAndVacuumShareOnePeriodicDensity)
{
    // A smooth periodic slab spans ordinary density, small active density and
    // inactive vacuum in the same cell. Missing cross-region coupling or a
    // wrong divergence changes the independent energy finite difference.
    const Rvv10::SerialEvaluator evaluator(6.3, 0.0093);
    for (int nz: {24, 32})
    {
        SCOPED_TRACE(testing::Message() << "nz=" << nz);
        ModulePW::PW_Basis basis("cpu", "double");
        basis.initgrids(10.0, ModuleBase::Matrix3(1, 0, 0, -0.5, std::sqrt(0.75), 0, 0, 0, 3), 6, 8, nz);
        basis.initparameters(false, 4.0, 1, true);
        basis.setuptransform();
        basis.collect_local_pw();
        const double dv = basis.omega / basis.nxyz;
        std::vector<double> total(basis.nrxx);
        std::vector<double> valence(basis.nrxx);
        for (int x = 0; x < basis.nx; ++x)
            for (int y = 0; y < basis.ny; ++y)
                for (int z = 0; z < nz; ++z)
                {
                    const int i = (x * basis.ny + y) * nz + z;
                    const double profile = std::exp(-16 * (1 - std::cos(2 * pi * z / nz)));
                    total[i] = 5.e-14 + 0.02 * profile * (1 + 0.1 * std::cos(2 * pi * x / basis.nx));
                    valence[i] = total[i] - 0.002 * profile;
                }
        ASSERT_GT(total[0], 1.e-2);
        ASSERT_LT(total[nz / 2], Rvv10::density_cutoff);
        const auto value = evaluator.evaluate(basis, total, valence);
        ASSERT_TRUE(std::isfinite(value.energy));
        double contraction = 0;
        for (int i = 0; i < basis.nrxx; ++i)
        {
            ASSERT_TRUE(std::isfinite(value.potential[i]));
            contraction += dv * valence[i] * value.potential[i];
        }
        EXPECT_NEAR(value.vtxc, contraction, 1.e-12);
        for (int point: {0, nz / 4})
        {
            SCOPED_TRACE(testing::Message() << "point=" << point);
            ASSERT_GT(total[point], 100 * Rvv10::density_cutoff);
            // In the dilute tail n~2.5e-9, h~1e-13 loses the energy
            // difference to roundoff. Use three larger, still local steps;
            // do not relax the potential tolerance or cross the density cut.
            const double relative_step = point == 0 ? 1.e-3 : 3.e-2;
            for (double factor: {1.0, 0.3, 0.1})
            {
                const double h = factor * relative_step * total[point];
                auto plus = total;
                auto minus = total;
                plus[point] += h;
                minus[point] -= h;
                // Energy only: valence is not an input to the energy functional.
                const double fd = (evaluator.evaluate(basis, plus, valence).energy
                                   - evaluator.evaluate(basis, minus, valence).energy)
                                  / (2 * h * dv);
                EXPECT_NEAR(value.potential[point], fd, 2.e-5) << "h=" << h;
            }
        }
    }
}
} // namespace

#ifdef __LIBXC
// The built-in PBE control still obtains its printable names from Libxc.
// Repeated calls expose missing ownership handling under LeakSanitizer.
TEST(XCOutputInfo, RepeatedBuiltinFunctionalNames)
{
    XC_Functional::set_xc_type("PBE");
    for (int repeat = 0; repeat < 32; ++repeat)
        EXPECT_EQ(XC_Functional::output_info(), " XC:\tgga_x_pbe\tgga_c_pbe\t");
}

TEST(Rvv10SCF, FunctionalSelectionRemainsIndependentOfNonlocalCorrection)
{
    XC_Functional::set_xc_type("GGA_X_RPW86+GGA_C_PBE");
    const auto base_ids = XC_Functional::get_func_id();
    EXPECT_EQ(XC_Functional::get_func_type(), 2);
    XC_Functional::set_xc_type("PBE");
    EXPECT_NE(base_ids, XC_Functional::get_func_id());
    XC_Functional::set_xc_type("GGA_X_RPW86+GGA_C_PBE");
    EXPECT_EQ(base_ids, XC_Functional::get_func_id());
    XC_Functional::set_xc_type("PBE");
}
#endif
