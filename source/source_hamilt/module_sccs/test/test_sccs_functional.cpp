#include "sccs_test.h"
#include "../sccs_functional.h"
#include "../sccs_parameters.h"
#include "../sccs_response.h"
#include "../sccs_pw_coulomb.h"
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
    ModuleSccs::solve_response(density, charge, config.cavity, solver, cold,
                               basis, tpiba, response);
    ModuleSccs::FunctionalResult result;
    ModuleSccs::eval_functional(charge, response, config, basis, tpiba, result);
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
    ModuleSccs::SccsConfig config = ModuleSccs::preset_config(ModuleSccs::Preset::Vacuum);
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
    ModuleSccs::eval_functional(charge, response, config, basis, tpiba, center);
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
    ModuleSccs::eval_functional(charge, response, config, basis, tpiba, plus);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        response.solute[ir] -= 2.0 * step * mode[ir];
    }
    ModuleSccs::eval_functional(charge, response, config, basis, tpiba, minus);
    const double fd = (plus.surface_energy + plus.volume_energy - minus.surface_energy - minus.volume_energy)
                      / (2.0 * step);
    EXPECT_NEAR(fd, predicted, 1e-8);
    const double expected_volume = 0.5 * basis.omega;
    EXPECT_NEAR(center.volume, expected_volume, 1e-10);
}

TEST_F(SccsFunctionalTest, ExplicitPeriodicOperatorPreservesResponseAndFunctional)
{
    ModuleSccs::SccsConfig config = ModuleSccs::preset_config(ModuleSccs::Preset::WaterNeutral);
    const std::vector<double> mode = cosine_mode(0);
    std::vector<double> charge(basis.nrxx);
    std::vector<double> density(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        charge[ir] = 1e-3 * mode[ir];
        density[ir] = 1e-3 + 1e-5 * mode[ir];
    }
    ModuleSccs::PolarizationSolverParameters solver;
    const std::vector<double> cold;
    ModuleSccs::SccsResponse legacy;
    ModuleSccs::SccsResponse explicit_response;
    ModuleSccs::solve_response(density, charge, config.cavity, solver, cold,
                               basis, tpiba, legacy);
    ModuleSccs::PeriodicCoulombOperator coulomb(basis, tpiba);
    EXPECT_FALSE(coulomb.has_correction());
    ModuleSccs::solve_response(density, charge, config.cavity, solver, cold,
                               basis, tpiba, coulomb, explicit_response);
    EXPECT_EQ(legacy.polarization.iterations, explicit_response.polarization.iterations);
    EXPECT_EQ(legacy.polarization.potential, explicit_response.polarization.potential);
    EXPECT_EQ(legacy.cavity_potential, explicit_response.cavity_potential);
    EXPECT_EQ(legacy.restart_potential, explicit_response.restart_potential);
    ModuleSccs::FunctionalResult legacy_result;
    ModuleSccs::FunctionalResult explicit_result;
    ModuleSccs::eval_functional(charge, legacy, config, basis, tpiba,
                                legacy_result);
    ModuleSccs::eval_functional(charge, explicit_response, config, basis, tpiba,
                                coulomb, explicit_result);
    EXPECT_DOUBLE_EQ(legacy_result.reaction_energy, explicit_result.reaction_energy);
    EXPECT_DOUBLE_EQ(legacy_result.surface_energy, explicit_result.surface_energy);
    EXPECT_DOUBLE_EQ(legacy_result.volume_energy, explicit_result.volume_energy);
    EXPECT_EQ(legacy_result.reaction_potential, explicit_result.reaction_potential);
    EXPECT_EQ(legacy_result.electron_potential, explicit_result.electron_potential);
}

// The cavity-density derivative is produced only for the core-Gaussian forces
// of the full cavity; the electron potential does not depend on it.
TEST_F(SccsFunctionalTest, CavityDensityDerivativeOnlyForTheFullCavity)
{
    const std::vector<double> mode = cosine_mode(0);
    std::vector<double> density(basis.nrxx);
    std::vector<double> charge(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        density[ir] = 1e-3 * (1.0 + 0.3 * mode[ir]);
        charge[ir] = -density[ir];
    }
    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 1e-4;
    config.cavity.density_max = 5e-3;
    config.cavity.epsilon_bulk = 5.0;
    config.surface_regularization = 1e-8;
    config.surface_tension = 1e-5;
    config.pressure = 1e-6;
    ModuleSccs::PolarizationSolverParameters solver;
    const std::vector<double> cold;
    ModuleSccs::SccsResponse response;
    ModuleSccs::solve_response(density, charge, config.cavity, solver, cold,
                               basis, tpiba, response);
    ModuleSccs::FunctionalResult electronic;
    ModuleSccs::eval_functional(charge, response, config, basis, tpiba, electronic);
    EXPECT_TRUE(electronic.cavity_potential.empty());
    config.core_electrons = true;
    ModuleSccs::FunctionalResult full;
    ModuleSccs::eval_functional(charge, response, config, basis, tpiba, full);
    ASSERT_EQ(full.cavity_potential.size(), static_cast<std::size_t>(basis.nrxx));
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_EQ(full.electron_potential[ir], electronic.electron_potential[ir]);
        const double expected = full.electron_potential[ir] + full.reaction_potential[ir];
        EXPECT_NEAR(full.cavity_potential[ir], expected, 1e-14);
    }
}
