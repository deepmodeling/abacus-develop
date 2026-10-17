#include "sccs_test.h"
#include "../sccs_functional.h"
#include "../sccs_parameters.h"
#include "../sccs_response.h"
#include "source_base/parallel_reduce.h"

using SccsFunctionalTest = SccsTest::PwTest;

TEST_F(SccsFunctionalTest, UniformDielectricEnergyAndElectronDerivative)
{
    const std::vector<double> mode = cosine_mode(0);
    std::vector<double> charge = mode;
    for (double& value : charge)
    {
        value *= 1e-3;
    }
    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 1e-4;
    config.cavity.density_max = 5e-3;
    config.cavity.epsilon_bulk = 5.0;
    config.surface_regularization = 1e-8;
    ModuleSccs::PolarizationSolverParameters solver;
    const std::vector<double> density(basis.nrxx, 0.0);
    const std::vector<double> cold;
    ModuleSccs::SccsResponse response;
    ModuleSccs::solve_sccs_response(density, charge, config.cavity, solver, cold,
                                    basis, tpiba, response);
    ModuleSccs::FunctionalResult result;
    ModuleSccs::evaluate_functional(charge, response, config, basis, tpiba, result);
    const double kernel = ModuleBase::FOUR_PI / (tpiba * tpiba);
    const double expected = 0.25 * basis.omega * 1e-6 * kernel * (0.2 - 1.0);
    EXPECT_NEAR(result.reaction_energy, expected, 1e-13);
    EXPECT_DOUBLE_EQ(result.surface_energy, 0.0);
    EXPECT_DOUBLE_EQ(result.volume_energy, 0.0);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double potential = 0.8 * kernel * charge[ir];
        EXPECT_NEAR(result.electron_potential[ir], potential, 1e-13);
    }
}

TEST_F(SccsFunctionalTest, NonElectrostaticDiscreteEnergyDerivative)
{
    ModuleSccs::SccsConfig config = ModuleSccs::make_sccs_config(ModuleSccs::Preset::Vacuum);
    config.surface_tension = 0.01;
    config.pressure = -0.001;
    config.surface_regularization = 0.02;
    const std::vector<double> mode = cosine_mode(0);
    ModuleSccs::SccsResponse response;
    response.solute.resize(basis.nrxx);
    response.dsolute_drho.assign(basis.nrxx, 1.0);
    response.polarization.potential.assign(basis.nrxx, 0.0);
    response.cavity_potential.assign(basis.nrxx, 0.0);
    const std::vector<double> charge(basis.nrxx, 0.0);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        response.solute[ir] = 0.5 + 0.2 * mode[ir];
    }
    ModuleSccs::FunctionalResult center;
    ModuleSccs::evaluate_functional(charge, response, config, basis, tpiba, center);
    // Avoid cancellation of grid-integrated energies at excessively small steps.
    const double step = 1e-5;
    double predicted = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        predicted += center.electron_potential[ir] * mode[ir];
        response.solute[ir] += step * mode[ir];
    }
    Parallel_Reduce::reduce_pool(predicted);
    predicted *= basis.omega / basis.nxyz;
    ModuleSccs::FunctionalResult plus;
    ModuleSccs::FunctionalResult minus;
    ModuleSccs::evaluate_functional(charge, response, config, basis, tpiba, plus);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        response.solute[ir] -= 2.0 * step * mode[ir];
    }
    ModuleSccs::evaluate_functional(charge, response, config, basis, tpiba, minus);
    const double fd = (plus.surface_energy + plus.volume_energy - minus.surface_energy - minus.volume_energy)
                      / (2.0 * step);
    EXPECT_NEAR(fd, predicted, 1e-8);
    const double expected_volume = 0.5 * basis.omega;
    EXPECT_NEAR(center.volume, expected_volume, 1e-10);
}
