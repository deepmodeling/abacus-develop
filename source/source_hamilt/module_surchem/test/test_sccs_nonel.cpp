#include "../sccs/sccs_cavity.h"
#include "../sccs/sccs_nonel.h"

#include "gtest/gtest.h"

#include <cmath>
#include <vector>

namespace
{

ModuleSccs::NonElectrostaticResult evaluate_from_density(const std::vector<double>& density,
                                                         const ModuleSccs::UniformGrid& grid,
                                                         const ModuleSccs::CavityParameters& cavity,
                                                         const ModuleSccs::NonElectrostaticParameters& parameters)
{
    std::vector<double> solute(density.size());
    std::vector<double> derivative(density.size());
    for (std::size_t index = 0; index < density.size(); ++index)
    {
        const ModuleSccs::CavityPoint point = ModuleSccs::evaluate_cavity(density[index], cavity);
        solute[index] = point.solute;
        derivative[index] = point.dsolute_drho;
    }
    return ModuleSccs::evaluate_non_electrostatic(grid, parameters, solute, derivative);
}

TEST(SccsNonElectrostatic, UniformSoluteHasNoRegularizedSurface)
{
    ModuleSccs::UniformGrid grid;
    grid.nx = 4;
    grid.ny = 3;
    grid.nz = 2;
    grid.spacing_x = 0.5;
    grid.spacing_y = 0.7;
    grid.spacing_z = 0.9;
    ModuleSccs::NonElectrostaticParameters parameters;
    parameters.surface_tension = 0.8;
    parameters.pressure = -0.2;
    parameters.surface_regularization = 1.0e-6;
    const std::vector<double> solute(24, 0.25);
    const std::vector<double> derivative(24, 1.0);

    const ModuleSccs::NonElectrostaticResult result
        = ModuleSccs::evaluate_non_electrostatic(grid, parameters, solute, derivative);
    const double cell_volume = 4.0 * 3.0 * 2.0 * 0.5 * 0.7 * 0.9;
    EXPECT_NEAR(result.surface, 0.0, 1.0e-15);
    EXPECT_NEAR(result.volume, 0.25 * cell_volume, 1.0e-14);
    for (std::size_t index = 0; index < result.density_potential.size(); ++index)
    {
        EXPECT_NEAR(result.density_potential[index], parameters.pressure, 1.0e-15);
    }
}

TEST(SccsNonElectrostatic, DensityPotentialIsTheDiscreteEnergyDerivative)
{
    ModuleSccs::UniformGrid grid;
    grid.nx = 4;
    grid.ny = 3;
    grid.nz = 2;
    grid.spacing_x = 0.6;
    grid.spacing_y = 0.7;
    grid.spacing_z = 0.8;
    ModuleSccs::CavityParameters cavity;
    cavity.density_min = 1.0e-4;
    cavity.density_max = 5.0e-3;
    cavity.epsilon_bulk = 78.3;
    ModuleSccs::NonElectrostaticParameters parameters;
    parameters.surface_tension = 0.003;
    parameters.pressure = -0.0002;
    parameters.surface_regularization = 1.0e-5;
    std::vector<double> density(24);
    std::vector<double> direction(24);
    for (std::size_t index = 0; index < density.size(); ++index)
    {
        density[index] = 2.0e-4 + 4.0e-3 * (0.2 + 0.6 * static_cast<double>(index) / density.size());
        direction[index] = std::sin(0.7 * static_cast<double>(index));
    }

    const ModuleSccs::NonElectrostaticResult center
        = evaluate_from_density(density, grid, cavity, parameters);
    const double step = 1.0e-8;
    std::vector<double> plus = density;
    std::vector<double> minus = density;
    for (std::size_t index = 0; index < density.size(); ++index)
    {
        plus[index] += step * direction[index];
        minus[index] -= step * direction[index];
    }
    const ModuleSccs::NonElectrostaticResult result_plus
        = evaluate_from_density(plus, grid, cavity, parameters);
    const ModuleSccs::NonElectrostaticResult result_minus
        = evaluate_from_density(minus, grid, cavity, parameters);
    const double finite = ((result_plus.surface_energy + result_plus.volume_energy)
                           - (result_minus.surface_energy + result_minus.volume_energy))
                          / (2.0 * step);
    const double volume_element = grid.spacing_x * grid.spacing_y * grid.spacing_z;
    double analytic = 0.0;
    for (std::size_t index = 0; index < density.size(); ++index)
    {
        analytic += center.density_potential[index] * direction[index] * volume_element;
    }
    EXPECT_NEAR(analytic, finite, std::abs(finite) * 2.0e-8);
}

TEST(SccsNonElectrostatic, RejectsMismatchedArrays)
{
    ModuleSccs::UniformGrid grid;
    grid.nx = 2;
    grid.ny = 2;
    grid.nz = 2;
    grid.spacing_x = 1.0;
    grid.spacing_y = 1.0;
    grid.spacing_z = 1.0;
    ModuleSccs::NonElectrostaticParameters parameters;
    parameters.surface_regularization = 1.0e-6;
    EXPECT_THROW(ModuleSccs::evaluate_non_electrostatic(grid,
                                                        parameters,
                                                        std::vector<double>(8),
                                                        std::vector<double>(7)),
                 std::invalid_argument);
}

} // namespace
