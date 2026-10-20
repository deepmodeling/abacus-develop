#include "sccs_sa_test.h"
#include "source_base/parallel_reduce.h"

using SccsSolventAware = SccsTest::SolventAwareTest;
using SccsTest::DensitySurface;
using SccsTest::origin;
using SccsTest::pool_dot;
using SccsTest::pool_max;
using SccsTest::probe_parameters;
using SccsTest::shifted;
using SccsTest::spectral_gradient;
using SccsTest::spectral_laplacian;

// E(s) = int g s_sa(s): the adjoint applied to g is dE/ds, including the
// nonlocal probe term that the pointwise (1 - f) g misses (by 27 percent).
TEST_F(SccsSolventAware, AdjointIsTheDerivativeOfTheFilledBoundary)
{
    ModuleSccs::SolventAwareParameters probe = probe_parameters(1.5);
    probe.filling_threshold = 0.45;
    probe.filling_spread = 0.1;
    set_up_cube(16.0, 20.0, probe);
    const std::vector<double> local = sphere(origin, 3.5, 0.6);
    const ModuleBase::Vector3<double> center(2.5, 1.0, 0.0);
    const std::vector<double> direction = sphere(center, 2.0, 0.8);
    std::vector<double> weight(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const ModuleBase::Vector3<double> r = displacement(ir);
        const double angle_x = ModuleBase::TWO_PI * r.x / length;
        const double angle_y = ModuleBase::TWO_PI * r.y / length;
        weight[ir] = 1.0 + 0.5 * std::cos(angle_x) + 0.3 * std::sin(angle_y);
    }
    ModuleSccs::FilledCavity cavity;
    cavity.active = true;
    cavity.local = local;
    cavity.filling = fill(local);
    const std::vector<double> adjoint = ModuleSccs::sa_adjoint(cavity, kernel, basis, weight);
    std::vector<double> pointwise(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        pointwise[ir] = (1.0 - cavity.filling.filling[ir]) * weight[ir];
    }
    const double step = 1.0e-5;
    const double backward = -step;
    const std::vector<double> plus = shifted(local, direction, step);
    const std::vector<double> minus = shifted(local, direction, backward);
    const ModuleSccs::SolventAwareBoundary filled_plus = fill(plus);
    const ModuleSccs::SolventAwareBoundary filled_minus = fill(minus);
    const double energy_plus = pool_dot(weight, filled_plus.boundary);
    const double energy_minus = pool_dot(weight, filled_minus.boundary);
    const double finite_difference = (energy_plus - energy_minus) / (2.0 * step);
    const double analytic = pool_dot(adjoint, direction);
    const double local_only = pool_dot(pointwise, direction);
    const double scale = std::abs(finite_difference);
    // Measured 9e-10 relative.
    const double tolerance = 1.0e-8 * scale;
    EXPECT_NEAR(analytic, finite_difference, tolerance);
    const double local_miss = std::abs(local_only - finite_difference);
    const double local_floor = 1.0e-2 * scale;
    EXPECT_GT(local_miss, local_floor);
}

// For a well-resolved boundary the Environ chain formulas agree with spectral
// derivatives of s_sa up to the grid error of the composition. The filling
// switches over filling_spread / |grad c|, which must span several grid
// points; with Environ's filling_spread 0.02 it usually does not.
TEST_F(SccsSolventAware, ChainDerivativesMatchSpectralDerivativesOfTheFilledBoundary)
{
    ModuleSccs::SolventAwareParameters probe = probe_parameters(1.5);
    probe.filling_threshold = 0.45;
    probe.filling_spread = 0.3;
    set_up_cube(16.0, 80.0, probe);
    const std::vector<double> local = sphere(origin, 3.5, 1.5);
    const ModuleSccs::SolventAwareBoundary filled = fill(local);
    std::vector<ModuleBase::Vector3<double>> gradient = spectral_gradient(local, basis, tpiba);
    std::vector<double> laplacian = spectral_laplacian(local, basis, tpiba);
    const std::vector<ModuleBase::Vector3<double>> fraction_gradient
        = ModuleSccs::convolve_grad(kernel, basis, gradient);
    const std::vector<double> fraction_laplacian = ModuleSccs::convolve_probe(kernel, basis, laplacian);
    ModuleSccs::sa_chain_derivs(local, filled, fraction_gradient, fraction_laplacian, gradient,
                                laplacian);

    const std::vector<ModuleBase::Vector3<double>> expected_gradient = spectral_gradient(filled.boundary, basis,
                                                                                         tpiba);
    const std::vector<double> expected_laplacian = spectral_laplacian(filled.boundary, basis, tpiba);
    double gradient_error = 0.0;
    double gradient_scale = 0.0;
    double laplacian_error = 0.0;
    double laplacian_scale = 0.0;
    double filling_change = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const ModuleBase::Vector3<double> gradient_difference = gradient[ir] - expected_gradient[ir];
        const double gradient_norm = expected_gradient[ir].norm();
        const double laplacian_difference = std::abs(laplacian[ir] - expected_laplacian[ir]);
        const double laplacian_norm = std::abs(expected_laplacian[ir]);
        const double change = filled.boundary[ir] - local[ir];
        gradient_error = std::max(gradient_error, gradient_difference.norm());
        gradient_scale = std::max(gradient_scale, gradient_norm);
        laplacian_error = std::max(laplacian_error, laplacian_difference);
        laplacian_scale = std::max(laplacian_scale, laplacian_norm);
        filling_change = std::max(filling_change, change);
    }
    EXPECT_GT(pool_max(basis, filling_change), 0.05);
    // Measured 8e-6 and 1e-4 of the maxima.
    const double gradient_tolerance = 5.0e-5 * pool_max(basis, gradient_scale);
    const double laplacian_tolerance = 5.0e-4 * pool_max(basis, laplacian_scale);
    EXPECT_LT(pool_max(basis, gradient_error), gradient_tolerance);
    EXPECT_LT(pool_max(basis, laplacian_error), laplacian_tolerance);
}

// A density saddle on a grid node inside the switching range: the continuum
// curvature gives -tr(H)/r there (about 2e8 at r = 1e-8); the discrete
// derivative stays of the size of its neighbours and remains exact.
TEST_F(SccsSolventAware, ChainSurfaceDerivativeStaysBoundedAtADensitySaddle)
{
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1e-4;
    cavity.density_max = 5e-3;
    cavity.epsilon_bulk = 1.0;
    cavity.solvent_aware = probe_parameters(0.5);
    cavity.solvent_aware.filling_threshold = 0.99;
    set_up_cube(12.0, 60.0, cavity.solvent_aware);
    const std::vector<double> density = gaussian_pair(1.5);
    const double regularization = 1.0e-8;
    const DensitySurface surface = surface_of_density(density, cavity, regularization);
    double saddle = 0.0;
    double largest = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double value = std::abs(surface.chain.density_potential[ir]);
        if (displacement(ir).norm() < 1.0e-12)
        {
            saddle = value;
        }
        else
        {
            largest = std::max(largest, value);
        }
    }
    Parallel_Reduce::reduce_pool(saddle);
    largest = pool_max(basis, largest);
    EXPECT_LT(saddle, largest);
    const ModuleBase::Vector3<double> center(0.3, 0.2, 0.1);
    const std::vector<double> direction = sphere(center, 1.0, 0.5);
    const double error = surface_derivative_error(density, direction, cavity, regularization, 1.0e-7);
    EXPECT_LT(error, 5.0e-6);
}
