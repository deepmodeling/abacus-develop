#include "../xc_functional.h"
#include "../libxc_abacus.h"
#include "source_base/constants.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_global.h"
#include "source_base/parallel_reduce.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "gtest/gtest.h"
#include <mpi.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <tuple>
#include <vector>

// The focused target links PW/XC objects; the fixture owns charge storage.
Charge::Charge() {}
Charge::~Charge() {}
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
SepPot::SepPot() {}
SepPot::~SepPot() {}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

namespace
{
int test_rank = 0;
int test_size = 1;
const int noncollinear_spin = 4;
const int variational_gga = 2;
const double hybrid_alpha = 0.0;
const double hse_omega = 0.11;
typedef std::tuple<double, double, ModuleBase::matrix> VxcResult;
enum Branch { Smooth, Negative, Saturated, Radial };

// Each branch runs through both the built-in and LibXC production dispatch.
class RealPwNcgga : public testing::TestWithParam<std::tuple<bool, int>>
{
  protected:
    ModulePW::PW_Basis pw;
    UnitCell cell;
    Charge charge;
    std::array<std::vector<double>, 4> density;
    std::array<std::vector<double>, 4> direction;
    std::array<double*, 4> pointers;
    std::vector<double> core;
    std::vector<std::complex<double>> core_g;

    double sum(double value)
    {
        Parallel_Reduce::reduce_pool(value);
        return value;
    }

    void SetUp() override
    {
        pw.initmpi(test_size, test_rank, MPI_COMM_WORLD);
        // Unequal MPI2 slabs and a finite PW cutoff exercise the projection.
        pw.initgrids(7.0, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 24, 10, 9);
        pw.initparameters(false, 80.0, 2, false);
        pw.setuptransform();
        pw.collect_local_pw();
        cell.omega = pw.omega;
        cell.tpiba = pw.tpiba;
        cell.magnet.lsign_ = false;
        for (int c = 0; c < noncollinear_spin; ++c)
        {
            density[c].resize(pw.nrxx);
            direction[c].resize(pw.nrxx);
            pointers[c] = density[c].data();
        }
        core.resize(pw.nrxx);
        core_g.resize(pw.npw);
        charge.rhopw = &pw;
        charge.nrxx = pw.nrxx;
        charge.nxyz = pw.nxyz;
        charge.ngmc = pw.npw;
        charge.nspin = noncollinear_spin;
        charge.rho = pointers.data();
        charge.rho_core = core.data();
        charge.rhog_core = core_g.data();
        const int branch = std::get<1>(GetParam());
        const double totals[] = {2.2, -1.55, 0.43, 0.030};
        const double magnitudes[] = {0.62, 0.45, 0.95, 4.0e-4};
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            const double x = ModuleBase::TWO_PI * (ir / (pw.ny * pw.nplane)) / pw.nx;
            const double y = ModuleBase::TWO_PI * ((ir / pw.nplane) % pw.ny) / pw.ny;
            const double z = ModuleBase::TWO_PI * (ir % pw.nplane + pw.startz_current) / pw.nz;
            density[0][ir] = totals[branch] * (1.0 + 0.08 * std::sin(x) + 0.03 * std::cos(y + z));
            const double magnitude = magnitudes[branch] * (1.0 + 0.1 * std::cos(5 * x) + 0.03 * std::sin(y));
            const double theta = 0.7 + 0.32 * std::sin(3 * x) + 0.18 * std::cos(7 * x) + 0.1 * std::cos(y);
            const double phi = 0.4 + 0.27 * std::cos(4 * x) - 0.16 * std::sin(6 * x) + 0.09 * std::sin(z);
            density[1][ir] = magnitude * std::sin(theta) * std::cos(phi);
            density[2][ir] = magnitude * std::sin(theta) * std::sin(phi);
            density[3][ir] = magnitude * std::cos(theta);
            core[ir] = 0.05 * totals[branch] * (1.0 + 0.1 * std::cos(3 * x + z));
            if (branch == Negative) ASSERT_LT(density[0][ir] + core[ir], 0.0);
            if (branch == Saturated) ASSERT_GT(magnitude, density[0][ir] + core[ir]);
            if (branch == Radial) ASSERT_LT(magnitude, 1e-3);
            for (int c = 0; c < noncollinear_spin; ++c)
                direction[c][ir] = 0.1 + 0.2 * std::cos((c + 2) * x + y) + 0.1 * std::sin(z + c);
        }
        pw.real2recip(core.data(), core_g.data());
        XC_Functional::set_xc_type(std::get<0>(GetParam()) ? "GGA_X_PBE+GGA_C_PBE" : "PBE");
    }

    VxcResult evaluate()
    {
        return XC_Functional::v_xc(pw.nrxx, &charge, &cell, noncollinear_spin,
                                   true, false, variational_gga, hybrid_alpha, hse_omega);
    }

    double inner(const ModuleBase::matrix& potential, const bool perturbation)
    {
        double value = 0.0;
        for (int c = 0; c < noncollinear_spin; ++c)
            for (int ir = 0; ir < pw.nrxx; ++ir)
                value += potential(c, ir) * (perturbation ? direction[c][ir] : density[c][ir]);
        return sum(value) * pw.omega / pw.nxyz;
    }
};

TEST_P(RealPwNcgga, GradientAndDivergenceAreAdjoints)
{
    std::vector<std::complex<double>> reciprocal(pw.npw);
    std::vector<ModuleBase::Vector3<double>> field(pw.nrxx);
    std::vector<ModuleBase::Vector3<double>> gradient(pw.nrxx);
    std::vector<double> divergence(pw.nrxx);
    for (int ir = 0; ir < pw.nrxx; ++ir)
        field[ir] = ModuleBase::Vector3<double>(direction[1][ir], direction[2][ir], direction[3][ir]);
    pw.real2recip(direction[0].data(), reciprocal.data());
    XC_Functional::grad_rho(reciprocal.data(), gradient.data(), &pw, pw.tpiba);
    XC_Functional::grad_dot(field.data(), divergence.data(), &pw, pw.tpiba);
    double identity = 0.0;
    double norm = 0.0;
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        const double left = field[ir] * gradient[ir];
        const double right = divergence[ir] * direction[0][ir];
        identity += left + right;
        norm += std::abs(left) + std::abs(right);
    }
    const double global_norm = sum(norm);
    EXPECT_GT(global_norm, 1e-4);
    EXPECT_LE(std::abs(sum(identity)), 5e-11 * global_norm);
}

TEST_P(RealPwNcgga, DensityAndCoreDerivatives)
{
    const VxcResult reference = evaluate();
    const ModuleBase::matrix& potential = std::get<2>(reference);
    EXPECT_NEAR(std::get<1>(reference), inner(potential, false), 2e-12 * std::max(1.0, std::abs(std::get<1>(reference))));
    // Four independent spin channels plus a core-density perturbation.
    const bool radial = std::get<1>(GetParam()) == Radial;
    const double step = radial ? 2e-4 : 2e-3;
    const int refinements = radial ? 10 : 3;
    for (int channel = 0; channel <= noncollinear_spin; ++channel)
    {
        SCOPED_TRACE(channel);
        const bool is_core = channel == noncollinear_spin;
        const int c = is_core ? 0 : channel;
        std::vector<double>& values = is_core ? core : density[c];
        const std::vector<double> original = values;
        double analytic = 0.0;
        for (int ir = 0; ir < pw.nrxx; ++ir)
            analytic += potential(c, ir) * direction[c][ir];
        analytic = sum(analytic) * pw.omega / pw.nxyz;
        double previous = 0.0;
        double best = 1.0;
        for (int level = 0; level < refinements; ++level)
        {
            const double epsilon = step / (1 << level);
            double energies[2];
            for (int sign = 0; sign < 2; ++sign)
            {
                for (int ir = 0; ir < pw.nrxx; ++ir)
                    values[ir] = original[ir] + (sign == 0 ? epsilon : -epsilon) * direction[c][ir];
                if (is_core) pw.real2recip(core.data(), core_g.data());
                energies[sign] = std::get<0>(evaluate());
            }
            std::copy(original.begin(), original.end(), values.begin());
            if (is_core) pw.real2recip(core.data(), core_g.data());
            const double error = std::abs((energies[0] - energies[1]) / (2 * epsilon) - analytic);
            const double scale = std::max(1.0, std::abs(analytic));
            if (level > 0 && level <= 2) EXPECT_LE(error, 0.4 * previous + 5e-9 * scale);
            best = std::min(best, error);
            if (level == refinements - 1) EXPECT_LE(best, 3e-8 * scale);
            previous = error;
        }
    }
}

TEST_P(RealPwNcgga, SixComponentStressDerivative)
{
    double energy = 0.0;
    double vtxc = 0.0;
    ModuleBase::matrix potential;
    std::vector<double> stress;
    XC_Functional::gradcorr(energy, vtxc, potential, &charge, &pw, &cell, stress,
                            true, noncollinear_spin, true, false, variational_gga, hybrid_alpha, hse_omega);
    ASSERT_EQ(stress.size(), 9U);
    const std::vector<ModuleBase::Vector3<double>> original(pw.gcar, pw.gcar + pw.npw);
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column <= row; ++column)
        {
            SCOPED_TRACE(row * 3 + column);
            const double analytic = sum(stress[row * 3 + column]) / pw.nxyz;
            double previous = 0.0;
            double best = 1.0;
            for (int level = 0; level < 4; ++level)
            {
                const double epsilon = 2e-3 / (1 << level);
                double energies[2];
                for (int sign = 0; sign < 2; ++sign)
                {
                    // Inverse of I + epsilon e_row e_column, holding density fixed.
                    const double strain = sign == 0 ? epsilon : -epsilon;
                    const double factor = row == column ? strain / (1 + strain) : strain;
                    for (int ig = 0; ig < pw.npw; ++ig)
                    {
                        pw.gcar[ig] = original[ig];
                        pw.gcar[ig][column] -= factor * original[ig][row];
                    }
                    energies[sign] = std::get<0>(evaluate());
                }
                std::copy(original.begin(), original.end(), pw.gcar);
                const double error = std::abs(-(energies[0] - energies[1]) / (2 * epsilon * pw.omega) - analytic);
                const double tolerance = 3e-4 * std::max(1e-10, std::abs(analytic)) + 1e-11;
                if (level > 0 && previous > 1e-10 && error > 1e-10)
                    EXPECT_LE(error, 0.4 * previous + tolerance);
                best = std::min(best, error);
                if (level == 3) EXPECT_LE(best, tolerance);
                previous = error;
            }
        }
}

TEST_P(RealPwNcgga, SpinRotationAndInversion)
{
    const VxcResult reference = evaluate();
    // Cyclic permutation is a proper global rotation; inversion is checked separately.
    for (int inversion = 0; inversion < 2; ++inversion)
    {
        const std::array<std::vector<double>, 4> original = density;
        for (int c = 1; c < noncollinear_spin; ++c)
            for (int ir = 0; ir < pw.nrxx; ++ir)
                density[c][ir] = (inversion == 0 ? 1 : -1) * original[c % 3 + 1][ir];
        const VxcResult rotated = evaluate();
        EXPECT_NEAR(std::get<0>(reference), std::get<0>(rotated), 1e-10);
        for (int c = 0; c < noncollinear_spin; ++c)
            for (int ir = 0; ir < pw.nrxx; ++ir)
                EXPECT_NEAR(std::get<2>(rotated)(c, ir),
                            (c == 0 ? 1 : (inversion == 0 ? 1 : -1))
                                * std::get<2>(reference)(c == 0 ? 0 : c % 3 + 1, ir), 1e-11);
        for (int c = 0; c < noncollinear_spin; ++c)
            std::copy(original[c].begin(), original[c].end(), density[c].begin());
    }
}

TEST_P(RealPwNcgga, ZeroMagnetizationIsFinite)
{
    for (int c = 1; c < noncollinear_spin; ++c)
        std::fill(density[c].begin(), density[c].end(), 0.0);
    const VxcResult result = evaluate();
    EXPECT_TRUE(std::isfinite(std::get<0>(result)));
    for (int c = 1; c < noncollinear_spin; ++c)
        for (int ir = 0; ir < pw.nrxx; ++ir)
            EXPECT_NEAR(std::get<2>(result)(c, ir), 0.0, 1e-13);
}

INSTANTIATE_TEST_SUITE_P(ProductionBranches, RealPwNcgga,
                        testing::Combine(testing::Bool(), testing::Values(Smooth, Negative, Saturated, Radial)));
} // namespace

int main(int argc, char** argv)
{
    int threads = 1;
    Parallel_Global::read_pal_param(argc, argv, test_size, threads, test_rank);
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    Parallel_Global::finalize_mpi();
    return result;
}
