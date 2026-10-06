#include "sccs_test.h"
#include "../sccs_ionic_force.h"
#include "../sccs_ionic_charge.h"
#include "../sccs_response.h"
#include "../sccs_parameters.h"
#include "../sccs_functional.h"

#include "source_cell/cell_tools.h"

using SccsIonicForceTest = SccsTest::PwTest;

TEST_F(SccsIonicForceTest, AnalyticFourierForceAndGaugeInvariance)
{
    std::vector<unitcell::AtomData> atoms(1);
    atoms[0].valence_charge = 2.0;
    atoms[0].position.x = length / 4.0;
    std::vector<double> potential = cosine_mode(0);
    const double spread = ModuleSccs::gaussian_ion_spread;
    std::vector<ModuleBase::Vector3<double>> forces;
    ModuleSccs::gaussian_ionic_force(atoms, potential, basis, tpiba, spread, forces);
    const double exponent = -0.25 * spread * spread * tpiba * tpiba;
    const double expected = 2.0 * tpiba * std::exp(exponent);
    EXPECT_NEAR(forces[0].x, expected, 1e-12);
    EXPECT_NEAR(forces[0].y, 0.0, 1e-12);
    EXPECT_NEAR(forces[0].z, 0.0, 1e-12);
    atoms[0].position.x += length;
    for (double& value : potential) { value += 4.0; }
    ModuleSccs::gaussian_ionic_force(atoms, potential, basis, tpiba, spread, forces);
    EXPECT_NEAR(forces[0].x, expected, 1e-12);
    potential.assign(basis.nrxx, 3.0);
    ModuleSccs::gaussian_ionic_force(atoms, potential, basis, tpiba, spread, forces);
    EXPECT_NEAR(forces[0].x, 0.0, 1e-12);
    EXPECT_NEAR(forces[0].y, 0.0, 1e-12);
    EXPECT_NEAR(forces[0].z, 0.0, 1e-12);
}

TEST_F(SccsIonicForceTest, ReactionEnergyFiniteDifferenceAtFixedNonuniformCavity)
{
    std::vector<unitcell::AtomData> atoms(2);
    atoms[0].valence_charge = 1.0;
    atoms[0].position = ModuleBase::Vector3<double>(2.1, 3.2, 4.3);
    atoms[1].valence_charge = 1.0;
    atoms[1].position = ModuleBase::Vector3<double>(6.4, 5.1, 2.7);
    std::vector<double> density = cosine_mode(0);
    for (double& value : density) { value = (2.0 / basis.omega) * (1.0 + 0.1 * value); }
    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 1e-3;
    config.cavity.density_max = 3e-3;
    config.cavity.epsilon_bulk = 5.0;
    config.surface_regularization = 1e-8;
    ModuleSccs::PolarizationSolverParameters solver;
    solver.tolerance_rms = 1e-13;
    solver.tolerance_max = 1e-12;
    const double spread = ModuleSccs::gaussian_ion_spread;
    const std::vector<double> cold_start;
    auto evaluate = [&](const std::vector<unitcell::AtomData>& displaced,
                        ModuleSccs::FunctionalResult& functional) {
        std::vector<double> charge;
        ModuleSccs::gaussian_ionic_density(displaced, basis, tpiba, spread, charge);
        for (int ir = 0; ir < basis.nrxx; ++ir) { charge[ir] -= density[ir]; }
        ModuleSccs::SccsResponse response;
        ModuleSccs::solve_sccs_response(density, charge, config.cavity, solver, cold_start, basis, tpiba, response);
        ModuleSccs::evaluate_functional(charge, response, config, basis, tpiba, functional);
    };
    ModuleSccs::FunctionalResult baseline;
    evaluate(atoms, baseline);
    std::vector<ModuleBase::Vector3<double>> forces;
    ModuleSccs::gaussian_ionic_force(atoms, baseline.reaction_potential, basis,
                                     tpiba, spread, forces);
    const double step = 1e-4;
    for (std::size_t ia = 0; ia < atoms.size(); ++ia)
    {
        for (int axis = 0; axis < 3; ++axis)
        {
            std::vector<unitcell::AtomData> plus = atoms;
            std::vector<unitcell::AtomData> minus = atoms;
            plus[ia].position[axis] += step;
            minus[ia].position[axis] -= step;
            ModuleSccs::FunctionalResult positive;
            ModuleSccs::FunctionalResult negative;
            evaluate(plus, positive);
            evaluate(minus, negative);
            const double finite_difference = -(positive.reaction_energy - negative.reaction_energy) / (2.0 * step);
            EXPECT_NEAR(forces[ia][axis], finite_difference, 1e-8);
        }
    }
    config.cavity.epsilon_bulk = 1.0;
    evaluate(atoms, baseline);
    ModuleSccs::gaussian_ionic_force(atoms, baseline.reaction_potential, basis,
                                     tpiba, spread, forces);
    for (const ModuleBase::Vector3<double>& force : forces)
    {
        EXPECT_NEAR(force.x, 0.0, 1e-12);
        EXPECT_NEAR(force.y, 0.0, 1e-12);
        EXPECT_NEAR(force.z, 0.0, 1e-12);
    }
}
