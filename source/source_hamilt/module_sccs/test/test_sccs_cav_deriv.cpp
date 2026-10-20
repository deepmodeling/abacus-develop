#include "sccs_test.h"
#include "../sccs_cav_deriv.h"
#include "../sccs_cavity.h"
#include "../sccs_response.h"

#include "source_base/parallel_reduce.h"

using SccsCavityDerivativesTest = SccsTest::PwTest;

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
    ModuleSccs::make_cav_derivs(density, cavity, basis, tpiba, open_boundary,
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
    ModuleSccs::make_cav_derivs(density, cavity, basis, tpiba, periodic_boundary,
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

    cavity.lowpass_p1 = basis.ggecut;
    cavity.lowpass_p2 = 0.5;
    ModuleSccs::make_cav_derivs(density, cavity, basis, tpiba, open_boundary,
                                response, derivatives);
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
        EXPECT_NEAR(coefficient[ir], expected, 1e-12);
        EXPECT_NEAR(derivatives.gradient[ir].x, gradient, 1e-12);
    }
}
