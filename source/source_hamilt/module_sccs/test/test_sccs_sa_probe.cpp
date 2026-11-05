#include "sccs_sa_test.h"
#include "source_base/parallel_reduce.h"

#include <complex>

using SccsSolventAware = SccsTest::SolventAwareTest;
using SccsTest::origin;
using SccsTest::pool_max;
using SccsTest::probe_parameters;
using SccsTest::setup_basis;

namespace
{
// Normalized transform of the isolated probe, 4 pi int r^2 p(r) j0(G r) dr
// over the same cut at width + 5 spread, divided by its G = 0 value.
double isolated_probe_transform(double g, double width, double spread)
{
    const double cutoff = width + 5.0 * spread;
    const int intervals = 40000;
    const double step = cutoff / intervals;
    double value = 0.0;
    double norm = 0.0;
    for (int i = 0; i <= intervals; ++i)
    {
        const double r = i * step;
        const double weight = (i == 0 || i == intervals) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
        const double argument = (r - width) / spread;
        const double radial = r * r * std::erfc(argument);
        const double phase = g * r;
        const double bessel = phase > 1.0e-12 ? std::sin(phase) / phase : 1.0;
        value += weight * radial * bessel;
        norm += weight * radial;
    }
    return value / norm;
}

// Kernel entry of G = (fx, fy, fz) in units of 2 pi/a, from whichever rank holds it.
double kernel_entry(const ModulePW::PW_Basis& basis, const std::vector<double>& kernel, double fx)
{
    double values[2] = {0.0, 0.0};
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const ModuleBase::Vector3<double>& g = basis.gdirect[ig];
        if (std::abs(g.x - fx) < 1.0e-12 && std::abs(g.y) < 1.0e-12 && std::abs(g.z) < 1.0e-12)
        {
            values[0] += kernel[ig];
            values[1] += 1.0;
        }
    }
    Parallel_Reduce::reduce_pool(values, 2);
    EXPECT_EQ(values[1], 1.0);
    return values[0];
}

// Largest deviation of the kernel from the isolated-probe transform over the
// G vectors with |G| <= maximum_g (Bohr^-1).
double oblique_probe_error(const ModuleBase::Matrix3& lattice,
                           double scale,
                           const ModuleSccs::SolventAwareParameters& parameters,
                           double maximum_g,
                           double& compared)
{
    ModulePW::PW_Basis basis("cpu", "double");
    setup_basis(basis, scale, lattice, 120.0);
    const std::vector<double> kernel = ModuleSccs::sa_probe_kernel(basis, lattice, scale, parameters);
    const double width = parameters.solvent_radius * parameters.radial_scale;
    const double tpiba = ModuleBase::TWO_PI / scale;
    double error = 0.0;
    compared = 0.0;
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const double g = tpiba * std::sqrt(basis.gg[ig]);
        if (g > maximum_g)
        {
            continue;
        }
        const double expected = isolated_probe_transform(g, width, parameters.radial_spread);
        const double deviation = std::abs(kernel[ig] - expected);
        error = std::max(error, deviation);
        compared += 1.0;
    }
    Parallel_Reduce::reduce_pool(compared);
    return pool_max(basis, error);
}
} // namespace

// In a cell wider than the probe the kernel is the transform of Environ's
// minimum-image probe, normalized on the grid, and it keeps constants.
TEST_F(SccsSolventAware, WideCellProbeMatchesEnvironMinimumImageProbe)
{
    // At 120 Ry the grid error of the sampled erfc edge is below 1e-10.
    const ModuleSccs::SolventAwareParameters probe = probe_parameters(3.0);
    set_up_cube(24.0, 120.0, probe);
    const double width = parameters.solvent_radius * parameters.radial_scale;
    const double spread = parameters.radial_spread;
    // Environ function_erfc cuts the probe at width + 5 spread.
    const double environ_cutoff = width + 5.0 * spread;
    std::vector<double> environ_probe(basis.nrxx, 0.0);
    double integral = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double distance = displacement(ir).norm();
        if (distance <= environ_cutoff)
        {
            const double argument = (distance - width) / spread;
            environ_probe[ir] = std::erfc(argument);
        }
        integral += environ_probe[ir];
    }
    Parallel_Reduce::reduce_pool(integral);
    const double volume = length * length * length;
    integral *= volume / basis.nxyz;
    std::vector<std::complex<double>> environ_g(basis.npw);
    basis.real2recip(environ_probe.data(), environ_g.data());
    double difference = 0.0;
    for (int ig = 0; ig < basis.npw; ++ig)
    {
        const std::complex<double> expected = environ_g[ig] * (volume / integral);
        const double deviation = std::abs(kernel[ig] - expected);
        difference = std::max(difference, deviation);
    }
    EXPECT_LT(pool_max(basis, difference), 1.0e-13);
    const double kernel_at_zero = kernel_entry(basis, kernel, 0.0);
    EXPECT_NEAR(kernel_at_zero, 1.0, 1.0e-12);
    const double first_shell = ModuleBase::TWO_PI / length;
    const double kernel_at_shell = kernel_entry(basis, kernel, 1.0);
    const double isolated_at_shell = isolated_probe_transform(first_shell, width, spread);
    EXPECT_NEAR(kernel_at_shell, isolated_at_shell, 1.0e-8);

    // Any constant: the unit-integral probe keeps it.
    const double constant_value = 0.37;
    const std::vector<double> constant(basis.nrxx, constant_value);
    const std::vector<double> convolved = ModuleSccs::convolve_probe(kernel, basis, constant);
    double constant_error = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double deviation = std::abs(convolved[ir] - constant_value);
        constant_error = std::max(constant_error, deviation);
    }
    EXPECT_LT(pool_max(basis, constant_error), 1.0e-12);
}

// Image counts must follow the cell heights, not the lattice-vector lengths:
// a hexagonal cell narrower than the probe, and a cell sheared to 10 degrees
// whose 2.1 Bohr height needs images n = +-4 along a 12 Bohr vector. Counts
// from the lengths miss probe weight (kernel error 1.9e-2, against 8e-15).
TEST_F(SccsSolventAware, ObliqueCellProbeSumsPeriodicImages)
{
    const double sqrt3_half = 0.5 * std::sqrt(3.0);
    const ModuleBase::Matrix3 hexagonal(1.0, 0.0, 0.0,
                                        -0.5, sqrt3_half, 0.0,
                                        0.0, 0.0, 1.2);
    const double shear_angle = 10.0 * ModuleBase::PI / 180.0;
    const double shear_cos = std::cos(shear_angle);
    const double shear_sin = std::sin(shear_angle);
    const ModuleBase::Matrix3 sheared(1.0, 0.0, 0.0,
                                      shear_cos, shear_sin, 0.0,
                                      0.0, 0.0, 1.0);
    const ModuleSccs::SolventAwareParameters probe = probe_parameters(3.0);
    double hexagonal_count = 0.0;
    const double hexagonal_error = oblique_probe_error(hexagonal, 9.0, probe, 1.5, hexagonal_count);
    double sheared_count = 0.0;
    const double sheared_error = oblique_probe_error(sheared, 12.0, probe, 1.5, sheared_count);
    EXPECT_GT(hexagonal_count, 20.0);
    EXPECT_GT(sheared_count, 20.0);
    EXPECT_LT(hexagonal_error, 1.0e-6);
    EXPECT_LT(sheared_error, 1.0e-6);
}

// A thick shell around a small void: the void cannot hold the probe and is
// filled, while the open solvent outside the shell stays solvent.
TEST_F(SccsSolventAware, FillsAnEnclosedVoidAndKeepsTheOpenSolvent)
{
    const ModuleSccs::SolventAwareParameters probe = probe_parameters(3.0);
    set_up_cube(30.0, 20.0, probe);
    const std::vector<double> outer = sphere(origin, 8.0, 0.5);
    const std::vector<double> inner = sphere(origin, 1.5, 0.5);
    std::vector<double> shell(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        shell[ir] = outer[ir] - inner[ir];
    }
    const ModuleSccs::SolventAwareBoundary filled = fill(shell);
    // center local value, center filled value, filled volume
    double sums[3] = {0.0, 0.0, 0.0};
    double outside_filled = 0.0;
    const double rounding = 1.0e-15;
    const double upper_bound = 1.0 + rounding;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double lower_bound = shell[ir] - rounding;
        EXPECT_GE(filled.boundary[ir], lower_bound);
        EXPECT_LE(filled.boundary[ir], upper_bound);
        const double distance = displacement(ir).norm();
        if (distance < 1.0e-12)
        {
            sums[0] += shell[ir];
            sums[1] += filled.boundary[ir];
        }
        if (distance > 10.5)
        {
            outside_filled = std::max(outside_filled, filled.boundary[ir]);
        }
        sums[2] += (filled.boundary[ir] - shell[ir]) * length * length * length / basis.nxyz;
    }
    Parallel_Reduce::reduce_pool(sums, 3);
    EXPECT_LT(sums[0], 1.0e-4);
    EXPECT_GT(sums[1], 0.999);
    EXPECT_LT(pool_max(basis, outside_filled), 1.0e-6);
    // The void holds about 4 pi 1.5^3 / 3 = 14 Bohr^3.
    EXPECT_GT(sums[2], 5.0);
    EXPECT_LT(sums[2], 40.0);
}
