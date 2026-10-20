#include "sccs_test.h"
#include "../sccs_cav_deriv.h"
#include "../sccs_cavity.h"
#include "../sccs_response.h"
#include "../sccs_sa.h"

#include "source_hamilt/module_xc/xc_functional.h"

#include "source_base/parallel_reduce.h"

using SccsCavityDerivativesTest = SccsTest::PwTest;

namespace
{
std::vector<ModuleBase::Vector3<double>> fft_gradient(const std::vector<double>& values,
                                                      const std::vector<double>& filter,
                                                      const ModulePW::PW_Basis& basis,
                                                      double tpiba)
{
    std::vector<std::complex<double>> values_g(basis.npw);
    basis.real2recip(values.data(), values_g.data());
    for (int ig = 0; ig < basis.npw && !filter.empty(); ++ig)
    {
        values_g[ig] *= filter[ig];
    }
    std::vector<ModuleBase::Vector3<double>> gradient(values.size());
    XC_Functional::grad_rho(values_g.data(), gradient.data(), &basis, tpiba);
    return gradient;
}

// grad s_sa = (1 - f) grad s + (1 - s) f' p * grad s for a given grad s.
std::vector<ModuleBase::Vector3<double>> filled_gradient(const ModuleSccs::SccsResponse& response,
                                                         const std::vector<double>& kernel,
                                                         const ModulePW::PW_Basis& basis,
                                                         const std::vector<ModuleBase::Vector3<double>>& local)
{
    const std::vector<ModuleBase::Vector3<double>> fraction = ModuleSccs::convolve_grad(kernel, basis,
                                                                                         local);
    std::vector<ModuleBase::Vector3<double>> gradient(local.size());
    for (std::size_t i = 0; i < local.size(); ++i)
    {
        const double empty = 1.0 - response.solvent_aware.filling.filling[i];
        const double solvent = 1.0 - response.solvent_aware.local[i];
        gradient[i] = local[i] * empty + fraction[i] * (solvent * response.solvent_aware.filling.dfilling[i]);
    }
    return gradient;
}
} // namespace

TEST_F(SccsCavityDerivativesTest, OpenBoundaryCoefficientUsesSwitchingFunctionFft)
{
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    cavity.epsilon_bulk = 5.0;
    const std::vector<double> mode = cosine_mode(0);
    std::vector<double> density(basis.nrxx);
    // Manufacture a single-shell boundary rather than deriving the expected
    // coefficient with the same FFT helpers as production.
    const double amplitude = 0.45;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double solute = 0.5 + amplitude * mode[ir];
        double lower = cavity.density_min;
        double upper = cavity.density_max;
        for (int iteration = 0; iteration < 64; ++iteration)
        {
            const double middle = 0.5 * (lower + upper);
            const ModuleSccs::CavityPoint point = ModuleSccs::evaluate_cavity(middle, cavity);
            if (point.solute < solute)
            {
                lower = middle;
            }
            else
            {
                upper = middle;
            }
        }
        density[ir] = 0.5 * (lower + upper);
    }
    const bool open_boundary = true;
    ModuleSccs::SccsResponse response;
    ModuleSccs::CavityDerivatives derivatives;
    const std::vector<double>& coefficient = derivatives.coefficient;
    ModuleSccs::make_cav_derivs(density, cavity, SccsTest::no_probe, basis, tpiba, open_boundary,
                                response, derivatives);
    const double log_bulk = std::log(cavity.epsilon_bulk);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const double angle = ModuleBase::TWO_PI * ix / basis.nx;
        const double solute = 0.5 + amplitude * mode[ir];
        const double log_epsilon = log_bulk * (1.0 - solute);
        const double epsilon = std::exp(log_epsilon);
        const double gradient = -amplitude * tpiba * std::sin(angle);
        const double laplacian = -amplitude * tpiba * tpiba * mode[ir];
        const double expected = epsilon * (-0.5 * log_bulk * laplacian
                                           + 0.25 * log_bulk * log_bulk * gradient * gradient)
                                / ModuleBase::FOUR_PI;
        EXPECT_NEAR(response.solute[ir], solute, 1e-14);
        EXPECT_NEAR(coefficient[ir], expected, 1e-12);
        const double expected_gradient = -log_bulk * gradient;
        EXPECT_NEAR(response.grad_log_epsilon[ir].x, expected_gradient, 1e-12);
        EXPECT_NEAR(response.grad_log_epsilon[ir].y, 0.0, 1e-12);
        EXPECT_NEAR(response.grad_log_epsilon[ir].z, 0.0, 1e-12);
    }
    const bool periodic_boundary = false;
    ModuleSccs::SccsResponse periodic;
    ModuleSccs::CavityDerivatives periodic_derivatives;
    const std::vector<double>& periodic_coefficient = periodic_derivatives.coefficient;
    ModuleSccs::make_cav_derivs(density, cavity, SccsTest::no_probe, basis, tpiba, periodic_boundary,
                                periodic, periodic_derivatives);
    double difference = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double coefficient_difference = coefficient[ir] - periodic_coefficient[ir];
        const double local_difference = std::abs(coefficient_difference);
        if (local_difference > difference)
        {
            difference = local_difference;
        }
    }
    Parallel_Reduce::reduce_max_pool(basis.poolnproc, difference);
    EXPECT_GT(difference, 1e-8);

}

namespace
{
// Density whose switching function is 0.5 + amplitude * mode, found by
// bisection of the cavity function at every grid point.
std::vector<double> shell_density(const ModuleSccs::CavityParameters& cavity,
                                  const std::vector<double>& mode,
                                  double amplitude)
{
    std::vector<double> density(mode.size());
    for (std::size_t ir = 0; ir < mode.size(); ++ir)
    {
        const double solute = 0.5 + amplitude * mode[ir];
        double lower = cavity.density_min;
        double upper = cavity.density_max;
        for (int iteration = 0; iteration < 64; ++iteration)
        {
            const double middle = 0.5 * (lower + upper);
            const ModuleSccs::CavityPoint point = ModuleSccs::evaluate_cavity(middle, cavity);
            if (point.solute < solute)
            {
                lower = middle;
            }
            else
            {
                upper = middle;
            }
        }
        density[ir] = 0.5 * (lower + upper);
    }
    return density;
}
} // namespace

// The lowpass multiplies the single cosine shell by f(G) = erfc(G^2 - 0.5)/2
// = erfc(0.5)/2 before the FFT derivatives of the switching function.
TEST_F(SccsCavityDerivativesTest, OpenBoundaryLowpassFiltersSwitchingFunctionDerivatives)
{
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    cavity.epsilon_bulk = 5.0;
    cavity.lowpass_p1 = basis.ggecut;
    cavity.lowpass_p2 = 0.5;
    const std::vector<double> mode = cosine_mode(0);
    const double amplitude = 0.45;
    const std::vector<double> density = shell_density(cavity, mode, amplitude);
    const bool open_boundary = true;
    ModuleSccs::SccsResponse response;
    ModuleSccs::CavityDerivatives derivatives;
    ModuleSccs::make_cav_derivs(density, cavity, SccsTest::no_probe, basis, tpiba, open_boundary,
                                response, derivatives);
    const double log_bulk = std::log(cavity.epsilon_bulk);
    const double weight = 0.5 * std::erfc(0.5);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        const double angle = ModuleBase::TWO_PI * ix / basis.nx;
        const double gradient = -amplitude * tpiba * std::sin(angle) * weight;
        const double laplacian = -amplitude * tpiba * tpiba * mode[ir] * weight;
        const double expected = response.epsilon[ir] * (-0.5 * log_bulk * laplacian
                                      + 0.25 * log_bulk * log_bulk * gradient * gradient)
                                / ModuleBase::FOUR_PI;
        EXPECT_NEAR(derivatives.coefficient[ir], expected, 1e-12);
        EXPECT_NEAR(derivatives.gradient[ir].x, gradient, 1e-12);
    }
}

// The filled dielectric: periodic cells chain through the density; open
// boundaries differentiate s_sa on the filtered grid with lowpass, otherwise
// fill the FFT gradient of the local s analytically.
TEST_F(SccsCavityDerivativesTest, FilledDielectricPathsFollowBoundaryAndLowpass)
{
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    cavity.epsilon_bulk = 5.0;
    cavity.solvent_aware.solvent_radius = 1.5;
    cavity.solvent_aware.filling_threshold = 0.67;
    cavity.solvent_aware.filling_spread = 0.1;
    const ModuleBase::Matrix3 lattice;
    const std::vector<double> kernel = ModuleSccs::sa_probe_kernel(basis, lattice, length, cavity.solvent_aware);
    const std::vector<double> mode = cosine_mode(0);
    std::vector<double> density(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        density[ir] = 1e-3 * (1.0 + 0.3 * mode[ir]);
    }
    const double log_bulk = std::log(cavity.epsilon_bulk);
    const std::vector<double> no_filter;
    for (int path = 0; path < 3; ++path)
    {
        const bool open_boundary = path > 0;
        const bool lowpass = path == 2;
        cavity.lowpass_p1 = lowpass ? basis.ggecut : -1.0;
        cavity.lowpass_p2 = lowpass ? 0.5 : -1.0;
        ModuleSccs::SccsResponse response;
        ModuleSccs::CavityDerivatives derivatives;
        ModuleSccs::make_cav_derivs(density, cavity, kernel, basis, tpiba, open_boundary, response,
                                    derivatives);
        double filled_volume = response.solvent_aware.volume;
        EXPECT_GT(filled_volume, 50.0);
        std::vector<ModuleBase::Vector3<double>> expected;
        if (lowpass)
        {
            expected = fft_gradient(response.solute, derivatives.filter, basis, tpiba);
        }
        else if (open_boundary)
        {
            const std::vector<ModuleBase::Vector3<double>> local = fft_gradient(response.solvent_aware.local, no_filter,
                                                                                basis, tpiba);
            expected = filled_gradient(response, kernel, basis, local);
        }
        else
        {
            std::vector<ModuleBase::Vector3<double>> local = fft_gradient(density, no_filter, basis, tpiba);
            for (int ir = 0; ir < basis.nrxx; ++ir)
            {
                local[ir] *= response.dsolute_drho[ir];
            }
            expected = filled_gradient(response, kernel, basis, local);
        }
        double error = 0.0;
        double scale = 0.0;
        for (int ir = 0; ir < basis.nrxx; ++ir)
        {
            const ModuleBase::Vector3<double> target = expected[ir] * (-log_bulk);
            const ModuleBase::Vector3<double> difference = response.grad_log_epsilon[ir] - target;
            const double target_norm = target.norm();
            error = std::max(error, difference.norm());
            scale = std::max(scale, target_norm);
            const double epsilon_exponent = log_bulk * (1.0 - response.solute[ir]);
            EXPECT_NEAR(response.epsilon[ir], std::exp(epsilon_exponent), 1e-12);
        }
        Parallel_Reduce::reduce_max_pool(basis.poolnproc, error);
        Parallel_Reduce::reduce_max_pool(basis.poolnproc, scale);
        EXPECT_GT(scale, 1e-2);
        const double tolerance = 1e-12 * scale;
        EXPECT_LT(error, tolerance);
        EXPECT_EQ(derivatives.gradient.empty(), !lowpass);
    }
}

// The solvent radius alone switches the filling on, the same on every pool
// rank: the other solvent-aware parameters without a radius give the plain
// cavity on every path.
TEST_F(SccsCavityDerivativesTest, SolventAwareParametersWithoutARadiusLeaveTheCavityUnfilled)
{
    ModuleSccs::CavityParameters plain;
    plain.density_min = 1e-4;
    plain.density_max = 5e-3;
    plain.epsilon_bulk = 5.0;
    ModuleSccs::CavityParameters parameters_only = plain;
    parameters_only.solvent_aware.radial_scale = 1.5;
    parameters_only.solvent_aware.radial_spread = 0.3;
    parameters_only.solvent_aware.filling_threshold = 0.6;
    parameters_only.solvent_aware.filling_spread = 0.05;
    const std::vector<double> mode = cosine_mode(0);
    std::vector<double> density(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        density[ir] = 1e-3 * (1.0 + 0.3 * mode[ir]);
    }
    for (int path = 0; path < 3; ++path)
    {
        const bool open_boundary = path > 0;
        const bool lowpass = path == 2;
        const double p1 = lowpass ? basis.ggecut : -1.0;
        const double p2 = lowpass ? 0.5 : -1.0;
        plain.lowpass_p1 = p1;
        plain.lowpass_p2 = p2;
        parameters_only.lowpass_p1 = p1;
        parameters_only.lowpass_p2 = p2;
        ModuleSccs::SccsResponse expected;
        ModuleSccs::CavityDerivatives expected_derivatives;
        ModuleSccs::make_cav_derivs(density, plain, SccsTest::no_probe, basis, tpiba, open_boundary,
                                    expected, expected_derivatives);
        ModuleSccs::SccsResponse response;
        ModuleSccs::CavityDerivatives derivatives;
        ModuleSccs::make_cav_derivs(density, parameters_only, SccsTest::no_probe, basis, tpiba,
                                    open_boundary, response, derivatives);
        EXPECT_FALSE(response.solvent_aware.enabled());
        EXPECT_EQ(response.solute, expected.solute);
        EXPECT_EQ(derivatives.coefficient, expected_derivatives.coefficient);
    }
}

// A filled cavity without its probe kernel stops the run instead of reading
// past the kernel; a forked death test cannot share the pool collectives of a
// multi-rank run.
TEST_F(SccsCavityDerivativesTest, SolventAwareCavityWithoutAProbeKernelStops)
{
    if (SccsTest::pool_size != 1)
    {
        GTEST_SKIP();
    }
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    cavity.epsilon_bulk = 5.0;
    cavity.solvent_aware.solvent_radius = 1.0;
    const std::vector<double> density(basis.nrxx, 1e-3);
    const bool open_boundary = false;
    ModuleSccs::SccsResponse response;
    ModuleSccs::CavityDerivatives derivatives;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(ModuleSccs::make_cav_derivs(density, cavity, SccsTest::no_probe, basis, tpiba, open_boundary,
                                            response, derivatives),
                ::testing::ExitedWithCode(1), "");
    testing::internal::GetCapturedStdout();
}
