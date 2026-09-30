#ifdef __MPI
#include "source_base/parallel_global.h"
#include <mpi.h>
#endif

#include "../surchem.h"
#include "../common/pw_grid.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"

#include <gtest/gtest.h>

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct StandalonePcc2dResult
{
    double energy_rydberg = 0.0;
    ModuleBase::matrix force;
};

// Standalone pcc_2d for two unit ions and a normalized Gaussian cloud of two
// electrons, all given in Cartesian bohr.
StandalonePcc2dResult run_standalone_pcc_2d(const ModuleBase::Matrix3& lattice,
                                            const std::vector<ModuleBase::Vector3<double>>& ions,
                                            const ModuleBase::Vector3<double>& cloud_center,
                                            const int axis)
{
    ModulePW::PW_Basis basis("cpu", "double");
#ifdef __MPI
    basis.initmpi(1, 0, POOL_WORLD);
#endif
    const double scale = 1.0;
    basis.initgrids(scale, lattice, 20.0);
    basis.initparameters(false, 20.0, 1, false);
    basis.setuptransform();
    basis.collect_local_pw();
    UnitCell cell;
    cell.lat0 = scale;
    cell.latvec = lattice;
    cell.omega = std::abs(lattice.Det());
    cell.ntype = 1;
    cell.nat = 2;
    cell.atoms = new Atom[1];
    cell.atoms[0].na = 2;
    cell.atoms[0].mass = 1.0;
    cell.atoms[0].ncpp.zv = 1.0;
    for (std::size_t ion = 0; ion < ions.size(); ++ion)
    {
        cell.atoms[0].tau.push_back(ions[ion]);
    }

    const std::vector<ModuleBase::Vector3<double>> positions
        = ModuleSurchem::pw_grid_positions(basis, cell.latvec, cell.lat0);
    const double volume_element = cell.omega / static_cast<double>(basis.nxyz);
    std::vector<double> electron_density(basis.nrxx, 0.0);
    double electron_sum = 0.0;
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        const ModuleBase::Vector3<double> offset = positions[ir] - cloud_center;
        electron_density[ir] = std::exp(-offset.norm2());
        electron_sum += electron_density[ir] * volume_element;
    }
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        electron_density[ir] *= 2.0 / electron_sum;
    }

    SurchemParameters parameters;
    parameters.pcc_boundary = ModulePcc::Boundary::Pcc2d;
    parameters.pcc_2d_axis = axis;
    parameters.expected_electron_count = 2.0;
    parameters.expected_ionic_charge = 2.0;
    surchem correction;
    correction.set_parameters(parameters);
    const double* density_channels[1] = {electron_density.data()};
    ModuleBase::matrix potential;
    correction.v_correction_pcc(cell, basis, 1, density_channels, potential);
    StandalonePcc2dResult result;
    result.energy_rydberg = correction.pcc_energy_rydberg();
    result.force.create(2, 3);
    correction.cal_force_pcc(cell, result.force);
    return result;
}

// Cyclic relabeling (x, y, z) -> (z, x, y) that moves the y axis onto z.
ModuleBase::Vector3<double> y_to_z(const ModuleBase::Vector3<double>& position)
{
    return ModuleBase::Vector3<double>(position.z, position.x, position.y);
}

TEST(HCorrPcc, StandalonePcc2dMatchesPointIonVacuumCorrection)
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
    UnitCell cell;
    cell.lat0 = scale;
    cell.latvec = lattice;
    cell.omega = 0.8 * 1.2 * scale * scale * scale;
    cell.ntype = 1;
    cell.nat = 1;
    cell.atoms = new Atom[1];
    cell.atoms[0].na = 1;
    cell.atoms[0].mass = 1.0;
    cell.atoms[0].ncpp.zv = 1.0;
    cell.atoms[0].tau.push_back(ModuleBase::Vector3<double>(0.4, 0.78, 0.5));

    SurchemParameters parameters;
    parameters.pcc_boundary = ModulePcc::Boundary::Pcc2d;
    parameters.pcc_2d_axis = 1;
    parameters.expected_electron_count = 1.0;
    parameters.expected_ionic_charge = 1.0;
    surchem correction;
    correction.set_parameters(parameters);
    EXPECT_TRUE(correction.uses_pcc());
    EXPECT_FALSE(correction.uses_sccs());
    const double volume_element = cell.omega / static_cast<double>(basis.nxyz);
    const int electron_plane_y = basis.ny / 2;
    std::vector<double> electron_density(basis.nrxx, 0.0);
    for (int ix = 0; ix < basis.nx; ++ix)
    {
        for (int iz_local = 0; iz_local < basis.nplane; ++iz_local)
        {
            const int index = (ix * basis.ny + electron_plane_y) * basis.nplane + iz_local;
            electron_density[index]
                = 1.0 / (static_cast<double>(basis.nx * basis.nz) * volume_element);
        }
    }
    const double* density_channels[1] = {electron_density.data()};
    ModuleBase::matrix potential;
    correction.v_correction_pcc(cell, basis, 1, density_channels, potential);
    ModulePcc::Pcc2dGeometry geometry
        = ModulePcc::pcc_2d_geometry(cell.latvec, cell.lat0, 1, 1.0e-10);
    geometry.origin = cell.atoms[0].tau[0].y * cell.lat0;
    const double electron_y
        = geometry.parameters.cell_length
          * static_cast<double>(electron_plane_y) / basis.ny;
    const double dipole_y = -ModulePcc::pcc_2d_relative_coordinate(ModuleBase::Vector3<double>(0.0, electron_y, 0.0), geometry);
    const double expected_energy = 2.0 * ModuleBase::PI * dipole_y * dipole_y / cell.omega;
    EXPECT_NEAR(correction.pcc_energy_rydberg(), 2.0 * expected_energy, 1.0e-12);
    EXPECT_TRUE(correction.validate_iteration_result());
    EXPECT_DOUBLE_EQ(surchem::Ael, 0.0);
    EXPECT_DOUBLE_EQ(surchem::Acav, 0.0);
    EXPECT_TRUE(std::isfinite(potential(0, 0)));
    // The out_pot 2 electrostatic potential takes the whole PCC potential.
    const std::vector<double>& electrostatic = correction.electrostatic_correction();
    ASSERT_EQ(electrostatic.size(), static_cast<std::size_t>(basis.nrxx));
    for (int ir = 0; ir < basis.nrxx; ++ir)
    {
        EXPECT_DOUBLE_EQ(electrostatic[ir], potential(0, ir));
    }
    ModuleBase::matrix force(1, 3);
    correction.cal_force_pcc(cell, force);
    EXPECT_TRUE(std::isfinite(force(0, 1)));
    for (int level = 0; level <= 2; ++level)
    {
        parameters.debug = level;
        correction.set_parameters(parameters);
        correction.v_correction_pcc(cell, basis, 1, density_channels, potential);
        std::ostringstream output;
        correction.write_sccs_iteration(output);
        const std::string text = output.str();
        EXPECT_EQ(text.empty(), level == 0);
        EXPECT_EQ(text.find("E_PCC/Ry") != std::string::npos, level > 0);
        EXPECT_EQ(text.find("PCC2D_MOMENTS") != std::string::npos, level == 2);
    }
    parameters.use_sccs = true;
    parameters.sccs_config.boundary = parameters.pcc_boundary;
    parameters.start_drho = 0.01;
    correction.set_parameters(parameters);
    EXPECT_FALSE(correction.sccs_is_active());
    EXPECT_TRUE(correction.uses_pcc());
    correction.v_correction_pcc(cell, basis, 1, density_channels, potential);
    EXPECT_NEAR(correction.pcc_energy_rydberg(), 2.0 * expected_energy, 1.0e-12);

    // An electron-count mismatch is reported; PCC keeps the grid charge.
    parameters.use_sccs = false;
    parameters.start_drho = 0.0;
    parameters.expected_electron_count = 1.0 + 1.0e-3;
    correction.set_parameters(parameters);
    EXPECT_NO_THROW(correction.v_correction_pcc(cell, basis, 1, density_channels, potential));
    EXPECT_NEAR(correction.pcc_energy_rydberg(), 2.0 * expected_energy, 1.0e-12);

    // Interleaved solvents must retain their own energy, including when another
    // instance is reconfigured or cleared before the electronic-state evaluation.
    const double first_energy = correction.pcc_energy_rydberg();
    surchem second;
    parameters.expected_electron_count = 0.5;
    second.set_parameters(parameters);
    EXPECT_THROW(second.pcc_energy_rydberg(), std::logic_error);
    std::vector<double> second_density = electron_density;
    for (double& density : second_density)
    {
        density *= 0.5;
    }
    const double* second_channels[1] = {second_density.data()};
    second.v_correction_pcc(cell, basis, 1, second_channels, potential);
    const double second_energy = second.pcc_energy_rydberg();
    EXPECT_GT(std::abs(second_energy - first_energy), 1.0e-6);
    EXPECT_DOUBLE_EQ(correction.pcc_energy_rydberg(), first_energy);
    correction.v_correction_pcc(cell, basis, 1, density_channels, potential);
    EXPECT_DOUBLE_EQ(second.pcc_energy_rydberg(), second_energy);
    second.clear();
    EXPECT_THROW(second.pcc_energy_rydberg(), std::logic_error);
    EXPECT_FALSE(second.validate_iteration_result());
    std::ostringstream cleared_summary;
    EXPECT_NO_THROW(second.write_iteration(cleared_summary, 0.1));
    EXPECT_TRUE(cleared_summary.str().empty());
    EXPECT_DOUBLE_EQ(correction.pcc_energy_rydberg(), first_energy);
    const SurchemParameters disabled_parameters;
    second.set_parameters(disabled_parameters);
    EXPECT_DOUBLE_EQ(second.pcc_energy_rydberg(), 0.0);
}

// The same slab open along y (pcc_2d_axis 1) and, after the cyclic
// relabeling (x, y, z) -> (z, x, y), along z (pcc_2d_axis 2) must give the
// same PCC energy and the same forces with permuted components.
TEST(HCorrPcc, Pcc2dOpenAxisIsEquivalentUnderCyclicRelabeling)
{
    const ModuleBase::Matrix3 open_y(8.0, 0.0, 0.0,
                                     0.0, 12.0, 0.0,
                                     0.0, 0.0, 10.0);
    const ModuleBase::Matrix3 open_z(10.0, 0.0, 0.0,
                                     0.0, 8.0, 0.0,
                                     0.0, 0.0, 12.0);
    const std::vector<ModuleBase::Vector3<double>> ions_y
        = {ModuleBase::Vector3<double>(4.0, 5.1, 5.0), ModuleBase::Vector3<double>(4.6, 7.3, 5.4)};
    const ModuleBase::Vector3<double> cloud_y(4.2, 6.0, 5.1);
    const std::vector<ModuleBase::Vector3<double>> ions_z = {y_to_z(ions_y[0]), y_to_z(ions_y[1])};
    const ModuleBase::Vector3<double> cloud_z = y_to_z(cloud_y);

    const StandalonePcc2dResult along_y = run_standalone_pcc_2d(open_y, ions_y, cloud_y, 1);
    const StandalonePcc2dResult along_z = run_standalone_pcc_2d(open_z, ions_z, cloud_z, 2);

    EXPECT_GT(std::abs(along_y.energy_rydberg), 1.0e-6);
    EXPECT_NEAR(along_z.energy_rydberg, along_y.energy_rydberg,
                1.0e-11 * std::abs(along_y.energy_rydberg));
    for (int ion = 0; ion < 2; ++ion)
    {
        EXPECT_GT(std::abs(along_y.force(ion, 1)), 1.0e-6);
        EXPECT_NEAR(along_z.force(ion, 2), along_y.force(ion, 1), 1.0e-11);
        EXPECT_NEAR(along_z.force(ion, 0), along_y.force(ion, 2), 1.0e-14);
        EXPECT_NEAR(along_z.force(ion, 1), along_y.force(ion, 0), 1.0e-14);
    }
}

} // namespace
