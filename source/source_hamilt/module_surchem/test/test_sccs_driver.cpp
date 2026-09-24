#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs_driver.h"
#include "../sccs_pw_charge.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
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

    const ModuleSccs::SccsResult reused
        = evaluate_uniform_charge(1.0, state, basis, lattice, length);
    EXPECT_EQ(reused.response.polarization.iterations, 1);
    for (std::size_t index = 0; index < reused.electron_potential_hartree.size(); ++index)
    {
        EXPECT_NEAR(reused.electron_potential_hartree[index],
                    reused.electrostatic.electron_potential[index],
                    1.0e-14);
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
    MPI_Init(&argc, &argv);
    MPI_Comm_split(MPI_COMM_WORLD, 0, 1, &POOL_WORLD);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Comm_free(&POOL_WORLD);
    MPI_Finalize();
#endif
    return result;
}
