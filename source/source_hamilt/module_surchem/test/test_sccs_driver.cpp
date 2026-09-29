#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_driver.h"
#include "../sccs/sccs_pw_charge.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

ModuleSccs::SccsResult evaluate_uniform_charge(const double net_charge,
                                               ModuleSccs::SccsState& state,
                                               ModulePW::PW_Basis& basis,
                                               const ModuleBase::Matrix3& lattice,
                                               const double length)
{
    const double volume = length * length * length;
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const double ionic_charge = 1.0;
    const double electron_count = ionic_charge - net_charge;
    const std::vector<double> electron_density(basis.nrxx, electron_count / volume);
    const std::vector<double> ionic_density(basis.nrxx, ionic_charge / volume);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(basis, lattice, length);
    const ModuleSccs::PccGeometry pcc
        = ModuleSccs::pcc_geometry(lattice, length, 1.0e-10);
    const ModuleBase::Vector3<double> origin = pcc.origin;

    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 1.0e-2;
    config.cavity.density_max = 2.0e-2;
    config.cavity.epsilon_bulk = 5.0;
    config.surface_regularization = 1.0e-6;
    config.boundary = ModuleSccs::Boundary::Pcc0d;
    config.max_iterations = 100;
    config.mixing = 0.7;
    config.tolerance_rms = 1.0e-14;
    config.tolerance_max = 1.0e-14;
    const ModuleSccs::Pcc2dGeometry pcc_2d_geometry;
    const ModuleSccs::SerialChargeReduction charge_reduction;
    const ModuleSccs::SerialPolarizationReduction polarization_reduction;
    return ModuleSccs::evaluate_pw_sccs(electron_density,
                                        ionic_density,
                                        electron_count,
                                        ionic_charge,
                                        1.0e-10,
                                        positions,
                                        origin,
                                        config,
                                        pcc,
                                        pcc_2d_geometry,
                                        basis,
                                        ModuleBase::TWO_PI / length,
                                        volume_element,
                                        charge_reduction,
                                        polarization_reduction,
                                        state);
}

TEST(SccsDriver, EvaluatesNeutralAndFixedChargePcc2dSources)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(0.8, 0.0, 0.0,
                                      0.0, 1.2, 0.0,
                                      0.0, 0.0, 1.0);
    const double scale = 10.0;
    basis.initgrids(scale, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    const double volume = 0.8 * 1.2 * scale * scale * scale;
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(basis, lattice, scale);
    const ModuleBase::Vector3<double> origin = ModuleSccs::cell_center(lattice, scale);
    const ModuleSccs::Pcc2dGeometry geometry
        = ModuleSccs::pcc_2d_geometry(lattice, scale, 1.0e-10);
    std::vector<double> electron_density(basis.nrxx, 1.0 / volume);
    std::vector<double> ionic_density(basis.nrxx, 0.0);
    for (int index = 0; index < basis.nrxx; ++index)
    {
        ionic_density[index]
            = (1.0 + 0.4 * std::sin(ModuleBase::TWO_PI * positions[index].y
                                    / geometry.parameters.cell_length_y))
              / volume;
    }

    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 1.0e-2;
    config.cavity.density_max = 2.0e-2;
    config.cavity.epsilon_bulk = 5.0;
    config.surface_regularization = 1.0e-6;
    config.boundary = ModuleSccs::Boundary::Pcc2d;
    config.max_iterations = 100;
    config.mixing = 0.7;
    config.tolerance_rms = 1.0e-13;
    config.tolerance_max = 1.0e-13;
    const ModuleSccs::PccGeometry pcc;
    const ModuleSccs::SerialChargeReduction charge_reduction;
    const ModuleSccs::SerialPolarizationReduction polarization_reduction;
    ModuleSccs::SccsState state;
    const ModuleSccs::SccsResult neutral
        = ModuleSccs::evaluate_pw_sccs(electron_density,
                                       ionic_density,
                                       1.0,
                                       1.0,
                                       1.0e-10,
                                       positions,
                                       origin,
                                       config,
                                       pcc,
                                       geometry,
                                       basis,
                                       ModuleBase::TWO_PI / scale,
                                       volume_element,
                                       charge_reduction,
                                       polarization_reduction,
                                       state);
    EXPECT_NEAR(neutral.charge.net_charge, 0.0, 1.0e-12);
    EXPECT_NEAR(neutral.solute_moments_2d.charge, 0.0, 1.0e-12);
    EXPECT_NEAR(neutral.polarization_moments_2d.charge, 0.0, 1.0e-12);
    EXPECT_GT(std::abs(neutral.solute_moments_2d.dipole_y), 1.0e-3);
    EXPECT_GT(neutral.smooth_vacuum_pcc_energy, 0.0);
    EXPECT_TRUE(state.valid);

    state.reset();
    ionic_density.assign(basis.nrxx, 1.1 / volume);
    const ModuleSccs::SccsResult cation
        = ModuleSccs::evaluate_pw_sccs(electron_density,
                                       ionic_density,
                                       1.0,
                                       1.1,
                                       1.0e-10,
                                       positions,
                                       origin,
                                       config,
                                       pcc,
                                       geometry,
                                       basis,
                                       ModuleBase::TWO_PI / scale,
                                       volume_element,
                                       charge_reduction,
                                       polarization_reduction,
                                       state);
    EXPECT_NEAR(cation.charge.net_charge, 0.1, 1.0e-12);
    EXPECT_NEAR(cation.solute_moments_2d.charge, 0.1, 1.0e-12);
    EXPECT_NEAR(cation.polarization_moments_2d.charge, -0.08, 1.0e-12);
    EXPECT_NEAR(cation.screened_moments_2d.charge, 0.02, 1.0e-12);
    EXPECT_NEAR(cation.electrostatic.reaction_energy,
                -0.8 * cation.smooth_vacuum_pcc_energy,
                2.0e-12);
    EXPECT_TRUE(state.valid);
}

TEST(SccsDriver, PeriodicSqrtCgWarmStartsFromStoredPotential)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    basis.initgrids(length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    const double volume = length * length * length;
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(basis, lattice, length);
    const ModuleBase::Vector3<double> origin = ModuleSccs::cell_center(lattice, length);
    // A cavity-crossing electron mode on a uniform ionic background.
    std::vector<double> electron_density(basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const int ix = ir / (basis.ny * basis.nplane);
        electron_density[ir] = 0.009 + 0.008 * std::cos(ModuleBase::TWO_PI * ix / basis.nx);
    }
    const double electron_count = 0.009 * volume;
    const std::vector<double> ionic_density(basis.nrxx, electron_count / volume);

    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 2.4e-3;
    config.cavity.density_max = 1.55e-2;
    config.cavity.epsilon_bulk = 78.3;
    config.surface_regularization = 1.0e-6;
    config.boundary = ModuleSccs::Boundary::Periodic;
    config.max_iterations = 200;
    config.mixing = 0.5;
    config.tolerance_rms = 1.0e-11;
    config.tolerance_max = 1.0e-10;
    const ModuleSccs::PccGeometry pcc;
    const ModuleSccs::Pcc2dGeometry pcc_2d;
    const ModuleSccs::SerialChargeReduction charge_reduction;
    const ModuleSccs::SerialPolarizationReduction polarization_reduction;
    const double tpiba = ModuleBase::TWO_PI / length;
    ModuleSccs::SccsState state;
    const ModuleSccs::SccsResult cold
        = ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, electron_count,
                                       electron_count, 1.0e-10, positions, origin, config,
                                       pcc, pcc_2d, basis, tpiba, volume_element,
                                       charge_reduction, polarization_reduction, state);
    ASSERT_TRUE(state.valid);
    ASSERT_EQ(state.potential.size(), static_cast<std::size_t>(basis.nrxx));
    EXPECT_FALSE(cold.response.polarization.warm_started);
    ASSERT_GT(cold.response.polarization.iterations, 1);

    const ModuleSccs::SccsResult warm
        = ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, electron_count,
                                       electron_count, 1.0e-10, positions, origin, config,
                                       pcc, pcc_2d, basis, tpiba, volume_element,
                                       charge_reduction, polarization_reduction, state);
    EXPECT_TRUE(warm.response.polarization.warm_started);
    EXPECT_LT(warm.response.polarization.iterations, cold.response.polarization.iterations);
    EXPECT_NEAR(warm.electrostatic.reaction_energy, cold.electrostatic.reaction_energy, 1.0e-10);

    state.reset();
    EXPECT_TRUE(state.potential.empty());
}

// A charged solute inside a resolved cavity: the far field of the PCC sqrt-CG
// solution must satisfy Gauss's law, -(1 - 1/epsilon) times the solute charge.
TEST(SccsDriver, ChargedPcc2dSqrtCgPolarizationSatisfiesGaussLaw)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 2.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double scale = 10.0;
    basis.initgrids(scale, lattice, 320.0);
    basis.initparameters(false, 320.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    const double volume = 2.0 * scale * scale * scale;
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(basis, lattice, scale);
    const ModuleBase::Vector3<double> origin = ModuleSccs::cell_center(lattice, scale);
    ModuleSccs::Pcc2dGeometry geometry = ModuleSccs::pcc_2d_geometry(lattice, scale, 1.0e-10);
    geometry.origin_y = origin.y;
    const double electron_count = 0.8;
    const double ionic_charge = 1.0;
    const double electron_width = 1.5;
    const double ion_width = 0.5;
    std::vector<double> electron_density(basis.nrxx);
    std::vector<double> ionic_density(basis.nrxx);
    double electron_sum = 0.0;
    double ionic_sum = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double dx = positions[ir].x - origin.x;
        const double dy = positions[ir].y - origin.y;
        const double dz = positions[ir].z - origin.z;
        const double r2 = dx * dx + dy * dy + dz * dz;
        electron_density[ir] = std::exp(-r2 / (electron_width * electron_width));
        ionic_density[ir] = std::exp(-r2 / (ion_width * ion_width));
        electron_sum += electron_density[ir] * volume_element;
        ionic_sum += ionic_density[ir] * volume_element;
    }
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        electron_density[ir] *= electron_count / electron_sum;
        ionic_density[ir] *= ionic_charge / ionic_sum;
    }

    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 1.0e-4;
    config.cavity.density_max = 5.0e-3;
    config.cavity.epsilon_bulk = 78.3;
    config.surface_regularization = 1.0e-6;
    config.boundary = ModuleSccs::Boundary::Pcc2d;
    config.max_iterations = 300;
    config.mixing = 0.5;
    config.tolerance_rms = 1.0e-12;
    config.tolerance_max = 1.0e-10;
    const ModuleSccs::PccGeometry pcc;
    const ModuleSccs::SerialChargeReduction charge_reduction;
    const ModuleSccs::SerialPolarizationReduction polarization_reduction;
    ModuleSccs::SccsState state;
    const ModuleSccs::SccsResult result
        = ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, electron_count,
                                       ionic_charge, 1.0e-10, positions, origin, config,
                                       pcc, geometry, basis, ModuleBase::TWO_PI / scale,
                                       volume_element, charge_reduction,
                                       polarization_reduction, state);
    const double expected = -(1.0 - 1.0 / config.cavity.epsilon_bulk) * 0.2;
    const double far_field = result.response.far_field_polarization_charge;
    const double density_integral = result.polarization_moments_2d.charge;
    std::cout << "PCC2D_SQRT_CG_GAUSS far_field " << far_field << " density_integral "
              << density_integral << " expected " << expected << " iterations "
              << result.response.polarization.iterations << std::endl;
    EXPECT_NEAR(result.solute_moments_2d.charge, 0.2, 1.0e-10);
    // Both converge with the grid: far field 2.6e-5 at 160 Ry and 2e-6 here,
    // the dielectric_of_potential integral 6.5e-5 at 160 Ry and 5e-6 here.
    EXPECT_NEAR(far_field, expected, 1.0e-5);
    EXPECT_NEAR(density_integral, expected, 3.0e-5);
    EXPECT_GT(result.response.polarization.iterations, 1);
}

// With electron density above density_max everywhere, no bulk solvent reaches
// the open boundary and Gauss's law for eps_bulk cannot hold: the PCC2D check
// must stop before the warm-start state is stored.
TEST(SccsDriver, Pcc2dStopsWhenBulkSolventDoesNotReachTheOpenBoundary)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 2.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double scale = 10.0;
    basis.initgrids(scale, lattice, 40.0);
    basis.initparameters(false, 40.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    const double volume = 2.0 * scale * scale * scale;
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSccs::pw_grid_positions(basis, lattice, scale);
    const ModuleBase::Vector3<double> origin = ModuleSccs::cell_center(lattice, scale);
    ModuleSccs::Pcc2dGeometry geometry = ModuleSccs::pcc_2d_geometry(lattice, scale, 1.0e-10);
    geometry.origin_y = origin.y;
    const double uniform_density = 1.0e-2;
    const double electron_count = uniform_density * volume;
    const double ionic_charge = electron_count + 1.0;
    const double ion_width = 0.5;
    const std::vector<double> electron_density(basis.nrxx, uniform_density);
    std::vector<double> ionic_density(basis.nrxx);
    double ionic_sum = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const double dx = positions[ir].x - origin.x;
        const double dy = positions[ir].y - origin.y;
        const double dz = positions[ir].z - origin.z;
        const double r2 = dx * dx + dy * dy + dz * dz;
        ionic_density[ir] = std::exp(-r2 / (ion_width * ion_width));
        ionic_sum += ionic_density[ir] * volume_element;
    }
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        ionic_density[ir] *= ionic_charge / ionic_sum;
    }

    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 1.0e-4;
    config.cavity.density_max = 5.0e-3;
    config.cavity.epsilon_bulk = 78.3;
    config.surface_regularization = 1.0e-6;
    config.boundary = ModuleSccs::Boundary::Pcc2d;
    config.max_iterations = 50;
    config.mixing = 0.5;
    config.tolerance_rms = 1.0e-12;
    config.tolerance_max = 1.0e-10;
    const ModuleSccs::PccGeometry pcc;
    const ModuleSccs::SerialChargeReduction charge_reduction;
    const ModuleSccs::SerialPolarizationReduction polarization_reduction;
    ModuleSccs::SccsState state;
    const double tpiba = ModuleBase::TWO_PI / scale;
    std::string message;
    try
    {
        ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, electron_count,
                                     ionic_charge, 1.0e-10, positions, origin, config, pcc,
                                     geometry, basis, tpiba, volume_element, charge_reduction,
                                     polarization_reduction, state);
    }
    catch (const std::runtime_error& error)
    {
        message = error.what();
    }
    EXPECT_NE(message.find("pcc_2d far-field polarization charge"), std::string::npos)
        << message;
    EXPECT_FALSE(state.valid);
    EXPECT_TRUE(state.potential.empty());
}

TEST(SccsDriver, PreservesChargeAndCombinesPccEnergyPotentialAndState)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    basis.initgrids(length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    ModuleSccs::SccsState state;
    const ModuleSccs::SccsResult cation
        = evaluate_uniform_charge(1.0, state, basis, lattice, length);
    ASSERT_EQ(cation.response.polarization.status, ModuleSccs::PolarizationStatus::Converged);
    EXPECT_NEAR(cation.charge.net_charge, 1.0, 1.0e-12);
    EXPECT_NEAR(cation.solute_moments.charge, 1.0, 1.0e-12);
    EXPECT_NEAR(cation.polarization_moments.charge, -0.8, 1.0e-12);
    EXPECT_NEAR(cation.screened_moments.charge, 0.2, 1.0e-12);
    EXPECT_NEAR(cation.electrostatic.reaction_energy,
                -0.8 * cation.smooth_vacuum_pcc_energy,
                2.0e-12);
    EXPECT_NEAR(cation.non_electrostatic.surface_energy, 0.0, 1.0e-14);
    EXPECT_NEAR(cation.non_electrostatic.volume_energy, 0.0, 1.0e-14);
    ASSERT_TRUE(state.valid);

    // In a uniform dielectric the preconditioner is exact, so the warm-start
    // step reproduces the stored solution without any CG iteration.
    const ModuleSccs::SccsResult reused
        = evaluate_uniform_charge(1.0, state, basis, lattice, length);
    EXPECT_TRUE(reused.response.polarization.warm_started);
    EXPECT_EQ(reused.response.polarization.iterations, 0);
    for (std::size_t index = 0; index < reused.electron_potential_hartree.size(); ++index)
    {
        EXPECT_NEAR(reused.electron_potential_hartree[index],
                    reused.electrostatic.electron_potential[index],
                    1.0e-14);
    }

    // An incompatible grid signature must reject the cached potential, even if
    // it contains invalid data. The result must match a clean cold start.
    state.tpiba *= 2.0;
    const double invalid_value = std::numeric_limits<double>::quiet_NaN();
    std::fill(state.potential.begin(), state.potential.end(), invalid_value);
    const ModuleSccs::SccsResult invalidated
        = evaluate_uniform_charge(1.0, state, basis, lattice, length);
    EXPECT_EQ(invalidated.response.polarization.iterations, cation.response.polarization.iterations);
    EXPECT_DOUBLE_EQ(invalidated.electrostatic.reaction_energy, cation.electrostatic.reaction_energy);
    for (std::size_t index = 0; index < cation.electron_potential_hartree.size(); ++index)
    {
        EXPECT_DOUBLE_EQ(invalidated.electron_potential_hartree[index], cation.electron_potential_hartree[index]);
    }

    state.reset();
    const ModuleSccs::SccsResult anion
        = evaluate_uniform_charge(-1.0, state, basis, lattice, length);
    EXPECT_NEAR(anion.charge.net_charge, -1.0, 1.0e-12);
    EXPECT_NEAR(anion.solute_moments.charge, -1.0, 1.0e-12);
    EXPECT_NEAR(anion.polarization_moments.charge, 0.8, 1.0e-12);
    EXPECT_NEAR(anion.electrostatic.reaction_energy,
                cation.electrostatic.reaction_energy,
                2.0e-12);
}

} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    int process_count = 1;
    int thread_count = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, rank);
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
