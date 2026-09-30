#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../sccs/sccs_driver.h"
#include "../common/pw_grid.h"
#include "../sccs/sccs_pw_coulomb.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/timer.h"
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

// Electronic solvent mode: no core-electron density in the cavity.
const std::vector<double> no_core_density;

// Center of the cell spanned by the scaled lattice rows.
ModuleBase::Vector3<double> cell_center(const ModuleBase::Matrix3& lattice, const double scale)
{
    const double x = 0.5 * scale * (lattice.e11 + lattice.e21 + lattice.e31);
    const double y = 0.5 * scale * (lattice.e12 + lattice.e22 + lattice.e32);
    const double z = 0.5 * scale * (lattice.e13 + lattice.e23 + lattice.e33);
    return ModuleBase::Vector3<double>(x, y, z);
}

ModuleSccs::SccsResult evaluate_uniform_charge(const double net_charge,
                                               ModuleSccs::SccsState& state,
                                               ModulePW::PW_Basis& basis,
                                               const ModuleBase::Matrix3& lattice,
                                               const double length)
{
    const double volume = length * length * length;
    const double ionic_charge = 1.0;
    const double electron_count = ionic_charge - net_charge;
    const std::vector<double> electron_density(basis.nrxx, electron_count / volume);
    const std::vector<double> ionic_density(basis.nrxx, ionic_charge / volume);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSurchem::pw_grid_positions(basis, lattice, length);
    const ModulePcc::PccGeometry pcc
        = ModulePcc::pcc_geometry(lattice, length, 1.0e-10);

    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 1.0e-2;
    config.cavity.density_max = 2.0e-2;
    config.cavity.epsilon_bulk = 5.0;
    config.surface_regularization = 1.0e-6;
    config.boundary = ModulePcc::Boundary::Pcc0d;
    config.max_iterations = 100;
    config.tolerance_rms = 1.0e-14;
    config.tolerance_max = 1.0e-14;
    const ModulePcc::Pcc2dGeometry pcc_2d_geometry;
    const ModuleSurchem::SerialChargeReduction charge_reduction;
    return ModuleSccs::evaluate_pw_sccs(electron_density,
                                        ionic_density,
                                        no_core_density,
                                        electron_count,
                                        ionic_charge,
                                        1.0e-10,
                                        positions,
                                        config,
                                        pcc,
                                        pcc_2d_geometry,
                                        basis,
                                        charge_reduction,
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
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSurchem::pw_grid_positions(basis, lattice, scale);
    const ModulePcc::Pcc2dGeometry geometry
        = ModulePcc::pcc_2d_geometry(lattice, scale, 1, 1.0e-10);
    std::vector<double> electron_density(basis.nrxx, 1.0 / volume);
    std::vector<double> ionic_density(basis.nrxx, 0.0);
    for (int index = 0; index < basis.nrxx; ++index)
    {
        ionic_density[index]
            = (1.0 + 0.4 * std::sin(ModuleBase::TWO_PI * positions[index].y
                                    / geometry.parameters.cell_length))
              / volume;
    }

    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 1.0e-2;
    config.cavity.density_max = 2.0e-2;
    config.cavity.epsilon_bulk = 5.0;
    config.surface_regularization = 1.0e-6;
    config.boundary = ModulePcc::Boundary::Pcc2d;
    config.max_iterations = 100;
    config.tolerance_rms = 1.0e-13;
    config.tolerance_max = 1.0e-13;
    const ModulePcc::PccGeometry pcc;
    const ModuleSurchem::SerialChargeReduction charge_reduction;
    ModuleSccs::SccsState state;
    const ModuleSccs::SccsResult neutral
        = ModuleSccs::evaluate_pw_sccs(electron_density,
                                       ionic_density,
                                       no_core_density,
                                       1.0,
                                       1.0,
                                       1.0e-10,
                                       positions,
                                       config,
                                       pcc,
                                       geometry,
                                       basis,
                                       charge_reduction,
                                       state);
    EXPECT_NEAR(neutral.charge.net_charge, 0.0, 1.0e-12);
    EXPECT_NEAR(neutral.solute_moments_2d.charge, 0.0, 1.0e-12);
    EXPECT_NEAR(neutral.polarization_moments_2d.charge, 0.0, 1.0e-12);
    EXPECT_GT(std::abs(neutral.solute_moments_2d.dipole), 1.0e-3);
    EXPECT_GT(neutral.smooth_vacuum_pcc_energy, 0.0);
    EXPECT_TRUE(state.valid);

    state = ModuleSccs::SccsState();
    ionic_density.assign(basis.nrxx, 1.1 / volume);
    const ModuleSccs::SccsResult cation
        = ModuleSccs::evaluate_pw_sccs(electron_density,
                                       ionic_density,
                                       no_core_density,
                                       1.0,
                                       1.1,
                                       1.0e-10,
                                       positions,
                                       config,
                                       pcc,
                                       geometry,
                                       basis,
                                       charge_reduction,
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
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSurchem::pw_grid_positions(basis, lattice, length);
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
    config.boundary = ModulePcc::Boundary::Periodic;
    config.max_iterations = 200;
    config.tolerance_rms = 1.0e-11;
    config.tolerance_max = 1.0e-10;
    const ModulePcc::PccGeometry pcc;
    const ModulePcc::Pcc2dGeometry pcc_2d;
    const ModuleSurchem::SerialChargeReduction charge_reduction;
    ModuleSccs::SccsState state;
    const ModuleSccs::SccsResult cold
        = ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, no_core_density, electron_count,
                                       electron_count, 1.0e-10, positions, config,
                                       pcc, pcc_2d, basis,
                                       charge_reduction, state);
    ASSERT_TRUE(state.valid);
    ASSERT_EQ(state.potential.size(), static_cast<std::size_t>(basis.nrxx));
    EXPECT_FALSE(cold.response.polarization.warm_started);
    ASSERT_GT(cold.response.polarization.iterations, 1);

    const ModuleSccs::SccsResult warm
        = ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, no_core_density, electron_count,
                                       electron_count, 1.0e-10, positions, config,
                                       pcc, pcc_2d, basis,
                                       charge_reduction, state);
    EXPECT_TRUE(warm.response.polarization.warm_started);
    EXPECT_LT(warm.response.polarization.iterations, cold.response.polarization.iterations);
    EXPECT_NEAR(warm.electrostatic.reaction_energy, cold.electrostatic.reaction_energy, 1.0e-10);

    state = ModuleSccs::SccsState();
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
        = ModuleSurchem::pw_grid_positions(basis, lattice, scale);
    const ModuleBase::Vector3<double> origin = cell_center(lattice, scale);
    ModulePcc::Pcc2dGeometry geometry = ModulePcc::pcc_2d_geometry(lattice, scale, 1, 1.0e-10);
    geometry.origin = origin.y;
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
    config.boundary = ModulePcc::Boundary::Pcc2d;
    config.max_iterations = 300;
    config.tolerance_rms = 1.0e-12;
    config.tolerance_max = 1.0e-10;
    const ModulePcc::PccGeometry pcc;
    const ModuleSurchem::SerialChargeReduction charge_reduction;
    ModuleSccs::SccsState state;
    const ModuleSccs::SccsResult result
        = ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, no_core_density, electron_count,
                                       ionic_charge, 1.0e-10, positions, config,
                                       pcc, geometry, basis,
                                       charge_reduction,
                                       state);
    const double expected = -(1.0 - 1.0 / config.cavity.epsilon_bulk) * 0.2;
    const double far_field = result.response.far_field_polarization_charge;
    const double density_integral = result.polarization_moments_2d.charge;
    std::cout << "PCC2D_SQRT_CG_GAUSS far_field " << far_field << " density_integral "
              << density_integral << " expected " << expected << " iterations "
              << result.response.polarization.iterations << std::endl;
    EXPECT_NEAR(result.solute_moments_2d.charge, 0.2, 1.0e-10);
    // Far field: 2.6e-5 at 160 Ry and 2e-6 here. The dielectric_of_potential
    // integral uses the FFT gradient of the corrected potential, as Environ,
    // which rings at the open boundary: 3.8e-5 here, a diagnostic only.
    EXPECT_NEAR(far_field, expected, 1.0e-5);
    EXPECT_NEAR(density_integral, expected, 1.0e-4);
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
        = ModuleSurchem::pw_grid_positions(basis, lattice, scale);
    const ModuleBase::Vector3<double> origin = cell_center(lattice, scale);
    ModulePcc::Pcc2dGeometry geometry = ModulePcc::pcc_2d_geometry(lattice, scale, 1, 1.0e-10);
    geometry.origin = origin.y;
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
    config.boundary = ModulePcc::Boundary::Pcc2d;
    config.max_iterations = 50;
    config.tolerance_rms = 1.0e-12;
    config.tolerance_max = 1.0e-10;
    const ModulePcc::PccGeometry pcc;
    const ModuleSurchem::SerialChargeReduction charge_reduction;
    ModuleSccs::SccsState state;
    std::string message;
    try
    {
        ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, no_core_density, electron_count,
                                     ionic_charge, 1.0e-10, positions, config, pcc,
                                     geometry, basis, charge_reduction,
                                     state);
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

    state = ModuleSccs::SccsState();
    const ModuleSccs::SccsResult anion
        = evaluate_uniform_charge(-1.0, state, basis, lattice, length);
    EXPECT_NEAR(anion.charge.net_charge, -1.0, 1.0e-12);
    EXPECT_NEAR(anion.solute_moments.charge, -1.0, 1.0e-12);
    EXPECT_NEAR(anion.polarization_moments.charge, 0.8, 1.0e-12);
    EXPECT_NEAR(anion.electrostatic.reaction_energy,
                cation.electrostatic.reaction_energy,
                2.0e-12);
}

// H3O+-like solute for the fixed-density derivative check: eight electrons
// in a Gaussian at the cell center and nine ionic charges displaced 0.4 bohr
// along y, so that a shift of the electrons changes the energy at first order.
struct CationSolute
{
    std::vector<double> electron_density;
    std::vector<double> ionic_density;
    // Zero-integral density directions: a rigid y shift and a breathing mode.
    std::vector<double> shift_mode;
    std::vector<double> breathing_mode;
};

CationSolute make_cation_solute(const std::vector<ModuleBase::Vector3<double>>& positions,
                                const ModuleBase::Vector3<double>& center,
                                const double volume_element)
{
    const double electron_count = 8.0;
    const double ionic_charge = 9.0;
    const double electron_width = 1.3;
    const double ion_width = 0.5;
    const double ion_offset_y = 0.4;
    const std::size_t size = positions.size();
    CationSolute solute;
    solute.electron_density.resize(size);
    solute.ionic_density.resize(size);
    solute.shift_mode.resize(size);
    solute.breathing_mode.resize(size);
    double electron_sum = 0.0;
    double ionic_sum = 0.0;
    for (std::size_t ir = 0; ir < size; ++ir)
    {
        const double dx = positions[ir].x - center.x;
        const double dy = positions[ir].y - center.y;
        const double dz = positions[ir].z - center.z;
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double scaled_r2 = r2 / (electron_width * electron_width);
        const double gaussian = std::exp(-scaled_r2);
        solute.electron_density[ir] = gaussian;
        solute.shift_mode[ir] = 2.0 * dy / (electron_width * electron_width) * gaussian;
        solute.breathing_mode[ir] = (scaled_r2 - 1.5) * gaussian;
        const double ion_dy = dy - ion_offset_y;
        const double ion_r2 = dx * dx + ion_dy * ion_dy + dz * dz;
        solute.ionic_density[ir] = std::exp(-ion_r2 / (ion_width * ion_width));
        electron_sum += solute.electron_density[ir] * volume_element;
        ionic_sum += solute.ionic_density[ir] * volume_element;
    }
    const double electron_scale = electron_count / electron_sum;
    const double ionic_scale = ionic_charge / ionic_sum;
    double shift_sum = 0.0;
    double breathing_sum = 0.0;
    for (std::size_t ir = 0; ir < size; ++ir)
    {
        solute.electron_density[ir] *= electron_scale;
        solute.shift_mode[ir] *= electron_scale;
        solute.breathing_mode[ir] *= electron_scale;
        solute.ionic_density[ir] *= ionic_scale;
        shift_sum += solute.shift_mode[ir] * volume_element;
        breathing_sum += solute.breathing_mode[ir] * volume_element;
    }
    // Remove the grid integral of each mode so the electron count stays fixed.
    for (std::size_t ir = 0; ir < size; ++ir)
    {
        solute.shift_mode[ir] -= shift_sum / electron_count * solute.electron_density[ir];
        solute.breathing_mode[ir] -= breathing_sum / electron_count * solute.electron_density[ir];
    }
    return solute;
}

ModuleSccs::SccsResult evaluate_cation(const std::vector<double>& electron_density,
                                       const std::vector<double>& ionic_density,
                                       const ModulePcc::Boundary boundary,
                                       const ModulePW::PW_Basis& basis,
                                       const std::vector<ModuleBase::Vector3<double>>& positions,
                                       const ModuleBase::Vector3<double>& center,
                                       const ModuleBase::Matrix3& lattice,
                                       const double scale,
                                       const double volume_element,
                                       const double lowpass_p1,
                                       const double lowpass_p2)
{
    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 2.0e-4;
    config.cavity.density_max = 3.5e-3;
    config.cavity.epsilon_bulk = 78.3;
    config.cavity.lowpass_p1 = lowpass_p1;
    config.cavity.lowpass_p2 = lowpass_p2;
    config.surface_regularization = 1.0e-8;
    config.boundary = boundary;
    config.max_iterations = 500;
    config.tolerance_rms = 1.0e-13;
    config.tolerance_max = 1.0e-11;
    ModulePcc::PccGeometry pcc;
    ModulePcc::Pcc2dGeometry pcc_2d;
    if (boundary == ModulePcc::Boundary::Pcc0d)
    {
        pcc = ModulePcc::pcc_geometry(lattice, scale, 1.0e-10);
        pcc.origin = center;
    }
    else
    {
        pcc_2d = ModulePcc::pcc_2d_geometry(lattice, scale, 1, 1.0e-10);
        pcc_2d.origin = center.y;
    }
    const ModuleSurchem::SerialChargeReduction charge_reduction;
    // A fresh state keeps every evaluation a cold start.
    ModuleSccs::SccsState state;
    return ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, no_core_density, 8.0, 9.0, 1.0e-10,
                                        positions, config, pcc, pcc_2d, basis,
                                        charge_reduction,
                                        state);
}

// Environ deriv_lowpass 10/5, validated at ecutrho 300-500 Ry.
const double test_lowpass_p1 = 10.0;
const double test_lowpass_p2 = 5.0;

void make_basis(const ModuleBase::Matrix3& lattice,
                const double scale,
                const double ecut,
                ModulePW::PW_Basis& basis)
{
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    basis.initgrids(scale, lattice, ecut);
    basis.initparameters(false, ecut, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
}

// With the switching lowpass the electronic potential must be the derivative
// of the discrete reaction energy at fixed ions: compare int v_el dn with
// central differences of E_R. Without the filter the same derivative is not
// grid-converged (pointwise 150-430 Ha at the cavity edge), so the
// continuum -eps'|grad v|^2/(8 pi) is printed and bounds the filtered one.
void check_cation_cavity_derivative(const ModulePcc::Boundary boundary,
                                    const ModuleBase::Matrix3& lattice,
                                    const double scale,
                                    const char* label)
{
    ModulePW::PW_Basis basis("cpu", "double");
    // At 400 Ry, as in production, the far-field Gauss check passes for PCC2D.
    make_basis(lattice, scale, 400.0, basis);
    const double volume = scale * scale * scale * lattice.Det();
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSurchem::pw_grid_positions(basis, lattice, scale);
    const ModuleBase::Vector3<double> center = cell_center(lattice, scale);
    const CationSolute solute = make_cation_solute(positions, center, volume_element);
    const ModuleSccs::SccsResult result
        = evaluate_cation(solute.electron_density, solute.ionic_density, boundary, basis,
                          positions, center, lattice, scale, volume_element, test_lowpass_p1,
                          test_lowpass_p2);
    double exact_maximum = 0.0;
    double continuum_maximum = 0.0;
    for (std::size_t ir = 0; ir < positions.size(); ++ir)
    {
        const double gradient_square = result.response.polarization.field.gradient[ir].norm2();
        const double continuum_cavity
            = -result.response.depsilon_drho[ir] * gradient_square / (8.0 * ModuleBase::PI);
        exact_maximum = std::max(exact_maximum, std::abs(result.response.cavity_potential[ir]));
        continuum_maximum = std::max(continuum_maximum, std::abs(continuum_cavity));
    }
    std::cout << "SCCS_CAVITY_POTENTIAL " << label << " max_exact " << exact_maximum
              << " max_continuum " << continuum_maximum << std::endl;
    // Measured 7.4 against 6.8 Ha.
    const double spike_limit = 2.0 * continuum_maximum;
    EXPECT_LT(exact_maximum, spike_limit);
    const std::vector<const std::vector<double>*> modes = {&solute.shift_mode,
                                                           &solute.breathing_mode};
    const char* mode_names[] = {"shift", "breathing"};
    // Truncation error O(step^2) is below 1e-8 here; the sqrt-CG residual adds less.
    const double step = 3.0e-5;
    for (std::size_t mode = 0; mode < modes.size(); ++mode)
    {
        const std::vector<double>& direction = *modes[mode];
        std::vector<double> plus(direction.size());
        std::vector<double> minus(direction.size());
        double exact = 0.0;
        double continuum = 0.0;
        for (std::size_t ir = 0; ir < direction.size(); ++ir)
        {
            plus[ir] = solute.electron_density[ir] + step * direction[ir];
            minus[ir] = solute.electron_density[ir] - step * direction[ir];
            exact += result.electrostatic.electron_potential[ir] * direction[ir] * volume_element;
            const double gradient_square = result.response.polarization.field.gradient[ir].norm2();
            const double continuum_potential
                = -result.electrostatic.reaction_potential[ir]
                  - result.response.depsilon_drho[ir] * gradient_square / (8.0 * ModuleBase::PI);
            continuum += continuum_potential * direction[ir] * volume_element;
        }
        const double energy_plus
            = evaluate_cation(plus, solute.ionic_density, boundary, basis, positions, center,
                              lattice, scale, volume_element, test_lowpass_p1, test_lowpass_p2)
                  .electrostatic.reaction_energy;
        const double energy_minus
            = evaluate_cation(minus, solute.ionic_density, boundary, basis, positions, center,
                              lattice, scale, volume_element, test_lowpass_p1, test_lowpass_p2)
                  .electrostatic.reaction_energy;
        const double finite_difference = (energy_plus - energy_minus) / (2.0 * step);
        std::cout << "SCCS_CAVITY_DERIVATIVE " << label << ' ' << mode_names[mode]
                  << " finite_difference " << finite_difference << " exact " << exact
                  << " error " << exact - finite_difference << " continuum_error "
                  << continuum - finite_difference << std::endl;
        // Measured errors are about 1e-9; the continuum potential misses by 1e-3.
        EXPECT_GT(std::abs(finite_difference), 1.0e-2);
        EXPECT_NEAR(exact, finite_difference, 1.0e-6);
    }
}

TEST(SccsDriver, Pcc0dLowpassElectronPotentialIsExactDerivativeOfDiscreteEnergy)
{
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    check_cation_cavity_derivative(ModulePcc::Boundary::Pcc0d, lattice, 12.0, "pcc0d");
}

TEST(SccsDriver, Pcc2dLowpassElectronPotentialIsExactDerivativeOfDiscreteEnergy)
{
    const double long_axis = 20.0 / 12.0;
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, long_axis, 0.0,
                                      0.0, 0.0, 1.0);
    check_cation_cavity_derivative(ModulePcc::Boundary::Pcc2d, lattice, 12.0, "pcc2d");
}

// Without the lowpass the PCC cavity potential is Environ's continuum
// -eps'|grad v|^2/(8 pi) with the FFT gradient of the solved potential.
TEST(SccsDriver, Pcc0dDefaultCavityPotentialIsEnvironContinuum)
{
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double scale = 12.0;
    ModulePW::PW_Basis basis("cpu", "double");
    make_basis(lattice, scale, 120.0, basis);
    const double volume = scale * scale * scale;
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSurchem::pw_grid_positions(basis, lattice, scale);
    const ModuleBase::Vector3<double> center = cell_center(lattice, scale);
    const CationSolute solute = make_cation_solute(positions, center, volume_element);
    const ModuleSccs::SccsResult result
        = evaluate_cation(solute.electron_density, solute.ionic_density,
                          ModulePcc::Boundary::Pcc0d, basis, positions, center, lattice, scale,
                          volume_element, -1.0, -1.0);
    const double tpiba = ModuleBase::TWO_PI / scale;
    const std::vector<ModuleBase::Vector3<double>> gradient
        = ModuleSccs::periodic_gradient(result.response.polarization.field.potential, basis, tpiba);
    double maximum_cavity = 0.0;
    for (std::size_t ir = 0; ir < positions.size(); ++ir)
    {
        const double gradient_square = gradient[ir].norm2();
        const double expected
            = -(result.response.depsilon_drho[ir] * gradient_square / (8.0 * ModuleBase::PI));
        const double expected_electron = -result.electrostatic.reaction_potential[ir] + expected;
        EXPECT_DOUBLE_EQ(result.response.cavity_potential[ir], expected);
        EXPECT_DOUBLE_EQ(result.electrostatic.electron_potential[ir], expected_electron);
        maximum_cavity = std::max(maximum_cavity, std::abs(expected));
    }
    EXPECT_GT(maximum_cavity, 1.0e-3);
}

TEST(SccsDriver, LowpassRequiresPccBoundary)
{
    ModuleSccs::SccsConfig config;
    config.cavity.density_min = 2.0e-4;
    config.cavity.density_max = 3.5e-3;
    config.cavity.epsilon_bulk = 78.3;
    config.cavity.lowpass_p1 = test_lowpass_p1;
    config.cavity.lowpass_p2 = test_lowpass_p2;
    config.surface_regularization = 1.0e-8;
    config.max_iterations = 10;
    config.tolerance_rms = 1.0e-10;
    config.tolerance_max = 1.0e-8;
    config.boundary = ModulePcc::Boundary::Periodic;
    EXPECT_THROW(ModuleSccs::validate_config(config), std::invalid_argument);
    config.boundary = ModulePcc::Boundary::Pcc0d;
    EXPECT_NO_THROW(ModuleSccs::validate_config(config));
    config.cavity.lowpass_p2 = -1.0;
    EXPECT_THROW(ModuleSccs::validate_config(config), std::invalid_argument);
}

// A pseudo-valence density that vanishes at the nucleus puts dielectric inside
// the atom in electronic mode. ENVIRON 'full' mode adds a core Gaussian there,
// which restores epsilon = 1 at the nucleus without changing the solute charge.
TEST(SccsDriver, FullSolventModeFillsTheNuclearCavityHole)
{
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double scale = 10.0;
    ModulePW::PW_Basis basis("cpu", "double");
    make_basis(lattice, scale, 80.0, basis);
    const double volume = scale * scale * scale;
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSurchem::pw_grid_positions(basis, lattice, scale);
    const ModuleBase::Vector3<double> center = cell_center(lattice, scale);
    const double shell_width = 1.2;
    const double core_spread = 0.5;
    std::vector<double> electron_density(basis.nrxx);
    std::vector<double> core_density(basis.nrxx);
    double electron_sum = 0.0;
    double core_sum = 0.0;
    int nucleus_index = 0;
    double nucleus_distance = 1.0e10;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const ModuleBase::Vector3<double> offset = positions[ir] - center;
        const double r2 = offset.norm2();
        electron_density[ir] = r2 * std::exp(-r2 / (shell_width * shell_width));
        core_density[ir] = std::exp(-r2 / (core_spread * core_spread));
        electron_sum += electron_density[ir] * volume_element;
        core_sum += core_density[ir] * volume_element;
        if (r2 < nucleus_distance)
        {
            nucleus_distance = r2;
            nucleus_index = ir;
        }
    }
    const double valence = 6.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        electron_density[ir] *= valence / electron_sum;
        core_density[ir] *= valence / core_sum;
    }
    // Neutral solute; the ions coincide with the core Gaussian.
    const std::vector<double>& ionic_density = core_density;
    ModuleSccs::SccsConfig config = ModuleSccs::water_preset(ModuleSccs::Preset::WaterAnion);
    config.tolerance_rms = 1.0e-12;
    config.tolerance_max = 1.0e-10;
    const ModulePcc::PccGeometry pcc;
    const ModulePcc::Pcc2dGeometry pcc_2d;
    const ModuleSurchem::SerialChargeReduction charge_reduction;
    ModuleSccs::SccsState electronic_state;
    const ModuleSccs::SccsResult electronic
        = ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, no_core_density, valence,
                                       valence, 1.0e-10, positions, config, pcc, pcc_2d,
                                       basis, charge_reduction,
                                       electronic_state);
    config.core_electrons = true;
    config.core_spread = core_spread;
    ModuleSccs::SccsState full_state;
    const ModuleSccs::SccsResult full
        = ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, core_density, valence,
                                       valence, 1.0e-10, positions, config, pcc, pcc_2d,
                                       basis, charge_reduction,
                                       full_state);
    EXPECT_GT(electronic.response.epsilon[nucleus_index], 10.0);
    EXPECT_DOUBLE_EQ(full.response.epsilon[nucleus_index], 1.0);
    // The core Gaussians shape the cavity only.
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_DOUBLE_EQ(full.charge.solute[ir], electronic.charge.solute[ir]);
    }
    config.core_electrons = false;
    ModuleSccs::SccsState mismatched_state;
    EXPECT_THROW(ModuleSccs::evaluate_pw_sccs(electron_density, ionic_density, core_density,
                                              valence, valence, 1.0e-10, positions,
                                              config, pcc, pcc_2d, basis,
                                              charge_reduction,
                                              mismatched_state),
                 std::invalid_argument);
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
    // Error-path tests throw inside timed functions and leave their
    // ModuleBase::timer entries running; production turns these exceptions
    // into WARNING_QUIT, so the timers are not under test here.
    ModuleBase::timer::disable();
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
