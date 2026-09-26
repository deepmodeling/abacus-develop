#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../surchem.h"
#include "../pcc/sccs_pcc_2d_coulomb.h"
#include "../sccs/sccs_pw_charge.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace
{

TEST(HCorrSccs, DispatchesDeferredSummaryOnlyWhenRequested)
{
    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.start_drho = 0.01;
    surchem solvent;
    solvent.set_parameters(parameters);
    std::ostringstream silent;
    solvent.write_iteration(silent, 0.1);
    EXPECT_TRUE(silent.str().empty());
    parameters.debug = 1;
    solvent.set_parameters(parameters);
    std::ostringstream summary;
    solvent.write_iteration(summary, 0.1);
    EXPECT_NE(summary.str().find("SCCS_DEFERRED DRHO 0.1"), std::string::npos);
}

TEST(HCorrSccs, DelaysActivationUntilDensityOrIterationThreshold)
{
    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.start_drho = 1.0e-2;
    parameters.start_nmax = 5;

    surchem solvent;
    solvent.set_parameters(parameters);
    EXPECT_TRUE(solvent.uses_sccs());
    EXPECT_FALSE(solvent.sccs_is_active());
    EXPECT_FALSE(solvent.try_activate_sccs(1, 2.0e-2));
    EXPECT_FALSE(solvent.sccs_is_active());
    EXPECT_TRUE(solvent.try_activate_sccs(2, 9.0e-3));
    EXPECT_TRUE(solvent.sccs_is_active());
    EXPECT_FALSE(solvent.try_activate_sccs(3, 2.0e-2));

    solvent.set_parameters(parameters);
    EXPECT_FALSE(solvent.try_activate_sccs(4, 2.0e-2));
    EXPECT_TRUE(solvent.try_activate_sccs(5, 2.0e-2));
    EXPECT_TRUE(solvent.sccs_is_active());
}

TEST(HCorrSccs, KeepsDelayedActivationLatchedWhenPotentialStorageIsCleared)
{
    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.start_drho = 1.0e-2;
    parameters.start_nmax = 5;

    surchem solvent;
    solvent.set_parameters(parameters);
    EXPECT_TRUE(solvent.try_activate_sccs(2, 9.0e-3));

    // LCAO rebuilds PotSurChem between ionic steps, whose destructor clears
    // the allocated potential storage. The one-way activation must survive it.
    solvent.clear();
    EXPECT_TRUE(solvent.sccs_is_active());
    EXPECT_FALSE(solvent.try_activate_sccs(1, 2.0e-2));
}

TEST(HCorrSccs, ConvertsHartreeResultToRydbergPotentialAndEnergy)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double length = 10.0;
    const double volume = length * length * length;
    basis.initgrids(length, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    UnitCell cell;
    cell.lat0 = length;
    cell.latvec = lattice;
    cell.omega = volume;
    cell.tpiba = ModuleBase::TWO_PI / length;
    cell.tpiba2 = cell.tpiba * cell.tpiba;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = new Atom[1];
    cell.atoms[0].na = 1;
    cell.atoms[0].mass = 1.0;
    cell.atoms[0].ncpp.zv = 1.0;
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.5, 0.5, 0.5));

    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.expected_electron_count = 0.0;
    parameters.expected_ionic_charge = 1.0;
    parameters.sccs_config.cavity.density_min = 1.0e-2;
    parameters.sccs_config.cavity.density_max = 2.0e-2;
    parameters.sccs_config.cavity.epsilon_bulk = 5.0;
    parameters.sccs_config.surface_regularization = 1.0e-6;
    parameters.sccs_config.boundary = ModuleSccs::Boundary::Pcc0d;
    parameters.pcc_boundary = ModuleSccs::Boundary::Pcc0d;
    parameters.debug = 2;
    parameters.sccs_config.max_iterations = 100;
    parameters.sccs_config.mixing = 0.7;
    parameters.sccs_config.tolerance_rms = 1.0e-14;
    parameters.sccs_config.tolerance_max = 1.0e-14;
    surchem solvent;
    solvent.set_parameters(parameters);
    std::vector<double> electron_density(basis.nrxx, 0.0);
    const double* density_channels[1] = {electron_density.data()};
    std::vector<double> local_potential(basis.nrxx, 0.0);
    ModuleBase::matrix potential;
    solvent.v_correction_sccs(cell,
                              basis,
                              1,
                              density_channels,
                              local_potential.data(),
                              potential);

    const ModuleSccs::SccsResult& result = solvent.sccs_result();
    EXPECT_NEAR(result.charge.net_charge, 1.0, 1.0e-12);
    EXPECT_NEAR(result.point_solute_moments.charge, 1.0, 1.0e-12);
    EXPECT_NEAR(result.vacuum_pcc_energy,
                0.5 * 2.837297479480619 / length,
                1.0e-12);
    EXPECT_NEAR(surchem::Ael,
                2.0 * (result.electrostatic.reaction_energy
                       + result.vacuum_pcc_energy
                       + result.ionic_shape_pcc_energy),
                1.0e-14);
    EXPECT_NEAR(surchem::Acav,
                2.0 * (result.non_electrostatic.surface_energy
                       + result.non_electrostatic.volume_energy),
                1.0e-14);
    std::ostringstream debug_output;
    solvent.write_sccs_iteration(debug_output);
    const std::string debug_text = debug_output.str();
    EXPECT_NE(debug_text.find("PCC0D_ORIGIN X/Bohr 5 Y/Bohr 5 Z/Bohr 5"),
              std::string::npos);
    EXPECT_NE(debug_text.find("PCC0D_MOMENTS SMOOTH Q/e "), std::string::npos);
    EXPECT_NE(debug_text.find("PCC0D_MOMENTS POINT Q/e "), std::string::npos);
    EXPECT_NE(debug_text.find("PCC0D_MOMENTS POLARIZATION Q/e "), std::string::npos);
    EXPECT_NE(debug_text.find("PCC0D_MOMENTS SCREENED Q/e "), std::string::npos);
    EXPECT_NE(debug_text.find("PX/eBohr "), std::string::npos);
    EXPECT_NE(debug_text.find("PY/eBohr "), std::string::npos);
    EXPECT_NE(debug_text.find("PZ/eBohr "), std::string::npos);
    EXPECT_NE(debug_text.find("QRR/eBohr2 "), std::string::npos);
    const std::size_t energy_begin = debug_text.find("PCC0D_ENERGY ");
    ASSERT_NE(energy_begin, std::string::npos);
    std::istringstream energy_stream(debug_text.substr(energy_begin));
    std::string label;
    double reaction_energy = 0.0;
    double smooth_energy = 0.0;
    double point_energy = 0.0;
    double shape_energy = 0.0;
    double used_energy = 0.0;
    energy_stream >> label >> label >> reaction_energy >> label >> smooth_energy
                  >> label >> point_energy >> label >> shape_energy
                  >> label >> used_energy;
    EXPECT_FALSE(energy_stream.fail());
    EXPECT_NEAR(reaction_energy, result.electrostatic.reaction_energy, 1.0e-11);
    EXPECT_NEAR(smooth_energy, result.smooth_vacuum_pcc_energy, 1.0e-11);
    EXPECT_NEAR(point_energy, result.vacuum_pcc_energy, 1.0e-11);
    EXPECT_NEAR(shape_energy, result.ionic_shape_pcc_energy, 1.0e-11);
    EXPECT_NEAR(used_energy, 2.0 * result.vacuum_pcc_energy, 1.0e-11);

    parameters.debug = 0;
    surchem silent_solvent;
    silent_solvent.set_parameters(parameters);
    std::ostringstream silent_output;
    silent_solvent.write_sccs_iteration(silent_output);
    silent_solvent.write_sccs_diagnostics(silent_output);
    EXPECT_TRUE(silent_output.str().empty());
    parameters.debug = 1;
    surchem quiet_solvent;
    quiet_solvent.set_parameters(parameters);
    ModuleBase::matrix quiet_potential;
    quiet_solvent.v_correction_sccs(cell,
                                    basis,
                                    1,
                                    density_channels,
                                    local_potential.data(),
                                    quiet_potential);
    std::ostringstream quiet_output;
    quiet_solvent.write_sccs_iteration(quiet_output);
    EXPECT_EQ(quiet_output.str().find("PCC0D_"), std::string::npos);

    ASSERT_EQ(potential.nr, 1);
    ASSERT_EQ(potential.nc, basis.nrxx);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_TRUE(std::isfinite(potential(0, ir)));
        EXPECT_NEAR(potential(0, ir), 2.0 * result.electron_potential_hartree[ir], 1.0e-14);
    }
}

TEST(HCorrSccs, AppliesNeutralPcc2dPointIonEnergyAndPotential)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(0.8, 0.0, 0.0,
                                      0.0, 1.2, 0.0,
                                      0.0, 0.0, 1.0);
    const double scale = 10.0;
    const double volume = 0.8 * 1.2 * scale * scale * scale;
    basis.initgrids(scale, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    UnitCell cell;
    cell.lat0 = scale;
    cell.latvec = lattice;
    cell.omega = volume;
    cell.tpiba = ModuleBase::TWO_PI / scale;
    cell.tpiba2 = cell.tpiba * cell.tpiba;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = new Atom[1];
    cell.atoms[0].na = 1;
    cell.atoms[0].mass = 1.0;
    cell.atoms[0].ncpp.zv = 1.0;
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.4, 0.78, 0.5));

    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.expected_electron_count = 1.0;
    parameters.expected_ionic_charge = 1.0;
    parameters.sccs_config.cavity.density_min = 1.0e-2;
    parameters.sccs_config.cavity.density_max = 2.0e-2;
    parameters.sccs_config.cavity.epsilon_bulk = 1.0;
    parameters.sccs_config.surface_regularization = 1.0e-6;
    parameters.sccs_config.boundary = ModuleSccs::Boundary::Pcc2d;
    parameters.pcc_boundary = ModuleSccs::Boundary::Pcc2d;
    parameters.sccs_config.max_iterations = 100;
    parameters.sccs_config.mixing = 0.7;
    parameters.sccs_config.tolerance_rms = 1.0e-14;
    parameters.sccs_config.tolerance_max = 1.0e-14;

    surchem solvent;
    solvent.set_parameters(parameters);
    const double volume_element = volume / static_cast<double>(basis.nxyz);
    const int electron_plane_y = basis.ny / 2;
    std::vector<double> electron_density(basis.nrxx, 0.0);
    for (int ix = 0; ix < basis.nx; ++ix)
    {
        for (int iz_local = 0; iz_local < basis.nplane; ++iz_local)
        {
            const int index
                = (ix * basis.ny + electron_plane_y) * basis.nplane + iz_local;
            electron_density[index]
                = 1.0
                  / (static_cast<double>(basis.nx * basis.nz) * volume_element);
        }
    }
    const double* density_channels[1] = {electron_density.data()};
    std::vector<double> local_potential(basis.nrxx, 0.0);
    ModuleBase::matrix potential;
    solvent.v_correction_sccs(cell,
                              basis,
                              1,
                              density_channels,
                              local_potential.data(),
                              potential);

    const ModuleSccs::SccsResult& result = solvent.sccs_result();
    ModuleSccs::Pcc2dGeometry geometry
        = ModuleSccs::pcc_2d_geometry(cell.latvec, cell.lat0, 1.0e-10);
    geometry.origin_y = cell.atoms[0].tau[0].y * cell.lat0;
    const double electron_y
        = geometry.parameters.cell_length_y
          * (static_cast<double>(electron_plane_y) + 0.5)
          / static_cast<double>(basis.ny);
    const double dipole_y
        = -ModuleSccs::pcc_2d_relative_y(electron_y, geometry);
    const double expected_energy = 2.0 * ModuleBase::PI * dipole_y * dipole_y / volume;
    EXPECT_NEAR(result.charge.net_charge, 0.0, 1.0e-12);
    EXPECT_NEAR(result.point_solute_moments_2d.charge, 0.0, 1.0e-12);
    EXPECT_NEAR(result.point_solute_moments_2d.dipole_y, dipole_y, 1.0e-12);
    EXPECT_NEAR(result.vacuum_pcc_energy, expected_energy, 1.0e-12);
    EXPECT_NEAR(result.electrostatic.reaction_energy, 0.0, 1.0e-14);
    EXPECT_NEAR(surchem::Ael, 2.0 * expected_energy, 1.0e-12);
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_TRUE(std::isfinite(potential(0, ir)));
        EXPECT_NEAR(potential(0, ir),
                    2.0 * result.electron_potential_hartree[ir],
                    1.0e-14);
    }
}

TEST(HCorrSccs, AppliesChargedPcc2dEnergyAndPotential)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const ModuleBase::Matrix3 lattice(0.8, 0.0, 0.0,
                                      0.0, 1.2, 0.0,
                                      0.0, 0.0, 1.0);
    const double scale = 10.0;
    const double volume = 0.8 * 1.2 * scale * scale * scale;
    basis.initgrids(scale, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();

    UnitCell cell;
    cell.lat0 = scale;
    cell.latvec = lattice;
    cell.omega = volume;
    cell.tpiba = ModuleBase::TWO_PI / scale;
    cell.tpiba2 = cell.tpiba * cell.tpiba;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = new Atom[1];
    cell.atoms[0].na = 1;
    cell.atoms[0].mass = 1.0;
    cell.atoms[0].ncpp.zv = 1.0;
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.4, 0.78, 0.5));

    SurchemParameters parameters;
    parameters.use_sccs = true;
    parameters.expected_electron_count = 0.8;
    parameters.expected_ionic_charge = 1.0;
    parameters.debug = 2;
    parameters.normalization_tolerance = 1.0e-10;
    parameters.sccs_config.cavity.density_min = 1.0e-2;
    parameters.sccs_config.cavity.density_max = 2.0e-2;
    parameters.sccs_config.cavity.epsilon_bulk = 5.0;
    parameters.sccs_config.surface_regularization = 1.0e-6;
    parameters.sccs_config.boundary = ModuleSccs::Boundary::Pcc2d;
    parameters.pcc_boundary = ModuleSccs::Boundary::Pcc2d;
    parameters.sccs_config.max_iterations = 100;
    parameters.sccs_config.mixing = 0.7;
    parameters.sccs_config.tolerance_rms = 1.0e-14;
    parameters.sccs_config.tolerance_max = 1.0e-14;

    surchem solvent;
    solvent.set_parameters(parameters);
    std::vector<double> electron_density(basis.nrxx, 0.8 / volume);
    const double* density_channels[1] = {electron_density.data()};
    std::vector<double> local_potential(basis.nrxx, 0.0);
    ModuleBase::matrix potential;
    solvent.v_correction_sccs(cell,
                              basis,
                              1,
                              density_channels,
                              local_potential.data(),
                              potential);

    const ModuleSccs::SccsResult& result = solvent.sccs_result();
    EXPECT_NEAR(result.charge.net_charge, 0.2, 1.0e-12);
    EXPECT_NEAR(result.point_solute_moments_2d.charge, 0.2, 1.0e-12);
    EXPECT_NEAR(result.solute_moments_2d.charge, 0.2, 1.0e-12);
    EXPECT_NEAR(result.polarization_moments_2d.charge, -0.16, 1.0e-12);
    EXPECT_NEAR(result.screened_moments_2d.charge, 0.04, 1.0e-12);
    EXPECT_TRUE(std::isfinite(result.vacuum_pcc_energy));
    EXPECT_TRUE(std::isfinite(result.electrostatic.reaction_energy));
    ModuleSccs::Pcc2dGeometry pcc_geometry
        = ModuleSccs::pcc_2d_geometry(cell.latvec, cell.lat0, 1.0e-10);
    pcc_geometry.origin_y = cell.atoms[0].tau[0].y * cell.lat0;
    const std::vector<ModuleBase::Vector3<double>> grid_positions
        = ModuleSccs::pw_grid_positions(basis, cell.latvec, cell.lat0);
    const ModuleSccs::SerialChargeReduction charge_reduction;
    const ModuleSccs::Pcc2dMoments smooth_ionic_moments
        = ModuleSccs::reduced_pcc_2d_density_moments(result.charge.ionic,
                                                     grid_positions,
                                                     cell.omega / basis.nxyz,
                                                     pcc_geometry,
                                                     charge_reduction);
    std::vector<ModuleSccs::PointCharge> ionic_points(1);
    ionic_points[0].charge = cell.atoms[0].ncpp.zv;
    ionic_points[0].position = cell.atoms[0].tau[0] * cell.lat0;
    const ModuleSccs::Pcc2dMoments point_ionic_moments
        = ModuleSccs::pcc_2d_point_charge_moments(ionic_points,
                                                  pcc_geometry);
    const double expected_ionic_shape_energy
        = ModuleSccs::pcc_2d_ionic_shape_energy(
            result.polarization_moments_2d.charge,
            smooth_ionic_moments,
            point_ionic_moments,
            pcc_geometry.parameters);
    EXPECT_NEAR(result.ionic_shape_pcc_energy,
                expected_ionic_shape_energy,
                1.0e-14);
    EXPECT_NEAR(surchem::Ael,
                2.0 * (result.electrostatic.reaction_energy
                       + result.vacuum_pcc_energy
                       + result.ionic_shape_pcc_energy),
                1.0e-14);
    std::ostringstream diagnostics;
    solvent.write_sccs_diagnostics(diagnostics);
    const std::string diagnostic_text = diagnostics.str();
    EXPECT_NE(diagnostic_text.find("SCCS_DIAGNOSTIC reaction_energy_hartree"),
              std::string::npos);
    EXPECT_NE(diagnostic_text.find("SCCS_DIAGNOSTIC smooth_vacuum_pcc_energy_hartree"),
              std::string::npos);
    EXPECT_NE(diagnostic_text.find("SCCS_DIAGNOSTIC point_vacuum_pcc_energy_hartree"),
              std::string::npos);
    EXPECT_NE(diagnostic_text.find("SCCS_DIAGNOSTIC ionic_shape_pcc_energy_hartree"),
              std::string::npos);
    EXPECT_NE(diagnostic_text.find("SCCS_DIAGNOSTIC electrostatic_energy_rydberg"),
              std::string::npos);
    EXPECT_NE(diagnostic_text.find("SCCS_DIAGNOSTIC screened_charge"),
              std::string::npos);
    std::ostringstream iteration_output;
    solvent.write_sccs_iteration(iteration_output);
    const std::string iteration_text = iteration_output.str();
    EXPECT_NE(iteration_text.find("SCCS_ITER "), std::string::npos);
    EXPECT_NE(iteration_text.find("SCCS_TIME/s "), std::string::npos);
    EXPECT_NE(iteration_text.find("E_SOL/Ry "), std::string::npos);
    EXPECT_NE(iteration_text.find("SCCS_MIXING VALUE "), std::string::npos);
    EXPECT_NE(iteration_text.find("RESTARTS "), std::string::npos);
    EXPECT_NE(iteration_text.find("PCC2D_MOMENTS "), std::string::npos);
    EXPECT_NE(iteration_text.find("Q_SMOOTH/e "), std::string::npos);
    EXPECT_NE(iteration_text.find("PY_SMOOTH/eBohr "), std::string::npos);
    EXPECT_NE(iteration_text.find("QYY_SMOOTH/eBohr2 "), std::string::npos);
    EXPECT_NE(iteration_text.find("Q_POINT/e "), std::string::npos);
    EXPECT_NE(iteration_text.find("PY_POINT/eBohr "), std::string::npos);
    EXPECT_NE(iteration_text.find("QYY_POINT/eBohr2 "), std::string::npos);
    EXPECT_NE(iteration_text.find("PCC2D_ENERGY "), std::string::npos);
    EXPECT_NE(iteration_text.find("REACTION/Ha "), std::string::npos);
    EXPECT_NE(iteration_text.find("PCC_SMOOTH/Ha "), std::string::npos);
    EXPECT_NE(iteration_text.find("PCC_POINT/Ha "), std::string::npos);
    EXPECT_NE(iteration_text.find("PCC_ION_SHAPE/Ha "), std::string::npos);
    EXPECT_NE(iteration_text.find("PCC_USED/Ry "), std::string::npos);
    EXPECT_NE(iteration_text.find(
                  "SCCS_ITER "
                  + std::to_string(result.response.polarization.iterations)),
              std::string::npos);
    std::istringstream iteration_stream(iteration_text);
    std::string iteration_label;
    int iteration_count = 0;
    std::string time_label;
    double elapsed_seconds = -1.0;
    std::string energy_label;
    double solvation_energy_rydberg = 0.0;
    iteration_stream >> iteration_label >> iteration_count
                     >> time_label >> elapsed_seconds
                     >> energy_label >> solvation_energy_rydberg;
    EXPECT_EQ(iteration_label, "SCCS_ITER");
    EXPECT_EQ(iteration_count, result.response.polarization.iterations);
    EXPECT_EQ(time_label, "SCCS_TIME/s");
    EXPECT_GE(elapsed_seconds, 0.0);
    const std::size_t time_begin = iteration_text.find("SCCS_TIME/s ") + 12;
    const std::size_t time_end = iteration_text.find(' ', time_begin);
    const std::string time_value = iteration_text.substr(time_begin, time_end - time_begin);
    const std::size_t decimal_point = time_value.find('.');
    ASSERT_NE(decimal_point, std::string::npos);
    EXPECT_EQ(time_value.size() - decimal_point - 1, 2);
    EXPECT_EQ(energy_label, "E_SOL/Ry");
    EXPECT_NEAR(solvation_energy_rydberg, surchem::Ael + surchem::Acav, 1.0e-7);

    parameters.debug = 0;
    surchem silent_solvent;
    silent_solvent.set_parameters(parameters);
    std::ostringstream silent_output;
    silent_solvent.write_sccs_iteration(silent_output);
    silent_solvent.write_sccs_diagnostics(silent_output);
    EXPECT_TRUE(silent_output.str().empty());
    parameters.debug = 1;
    surchem quiet_solvent;
    quiet_solvent.set_parameters(parameters);
    ModuleBase::matrix quiet_potential;
    quiet_solvent.v_correction_sccs(cell,
                                    basis,
                                    1,
                                    density_channels,
                                    local_potential.data(),
                                    quiet_potential);
    std::ostringstream quiet_iteration_output;
    quiet_solvent.write_sccs_iteration(quiet_iteration_output);
    const std::string quiet_iteration_text = quiet_iteration_output.str();
    EXPECT_NE(quiet_iteration_text.find("SCCS_ITER "), std::string::npos);
    EXPECT_EQ(quiet_iteration_text.find("SCCS_MIXING "), std::string::npos);
    EXPECT_EQ(quiet_iteration_text.find("PCC2D_MOMENTS "), std::string::npos);
    EXPECT_EQ(quiet_iteration_text.find("PCC2D_ENERGY "), std::string::npos);
    std::ostringstream quiet_diagnostics;
    quiet_solvent.write_sccs_diagnostics(quiet_diagnostics);
    EXPECT_TRUE(quiet_diagnostics.str().empty());

    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_TRUE(std::isfinite(potential(0, ir)));
        EXPECT_NEAR(potential(0, ir),
                    2.0 * result.electron_potential_hartree[ir],
                    1.0e-14);
    }
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
