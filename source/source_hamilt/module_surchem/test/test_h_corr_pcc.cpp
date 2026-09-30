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
#include <string>
#include <vector>

namespace
{

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
    EXPECT_NEAR(surchem::Epcc, 2.0 * expected_energy, 1.0e-12);
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
    EXPECT_NEAR(surchem::Epcc, 2.0 * expected_energy, 1.0e-12);

    // An electron-count mismatch is reported; PCC keeps the grid charge.
    parameters.use_sccs = false;
    parameters.start_drho = 0.0;
    parameters.expected_electron_count = 1.0 + 1.0e-3;
    correction.set_parameters(parameters);
    EXPECT_NO_THROW(correction.v_correction_pcc(cell, basis, 1, density_channels, potential));
    EXPECT_NEAR(surchem::Epcc, 2.0 * expected_energy, 1.0e-12);
}

} // namespace
