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
    parameters.pcc_boundary = ModuleSccs::Boundary::Pcc2d;
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
    ModuleSccs::Pcc2dGeometry geometry
        = ModuleSccs::pcc_2d_geometry(cell.latvec, cell.lat0, 1.0e-10);
    geometry.origin_y = cell.atoms[0].tau[0].y * cell.lat0;
    const double electron_y
        = geometry.parameters.cell_length_y
          * (static_cast<double>(electron_plane_y) + 0.5) / basis.ny;
    const double dipole_y = -ModuleSccs::pcc_2d_relative_y(electron_y, geometry);
    const double expected_energy = 2.0 * ModuleBase::PI * dipole_y * dipole_y / cell.omega;
    EXPECT_NEAR(surchem::Ael, 2.0 * expected_energy, 1.0e-12);
    EXPECT_DOUBLE_EQ(surchem::Acav, 0.0);
    EXPECT_TRUE(std::isfinite(potential(0, 0)));
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
        EXPECT_EQ(text.find("PCC_TIME/s") != std::string::npos, level > 0);
        EXPECT_EQ(text.find("PCC2D_MOMENTS") != std::string::npos, level == 2);
    }
    parameters.use_sccs = true;
    parameters.sccs_config.boundary = parameters.pcc_boundary;
    parameters.start_drho = 0.01;
    correction.set_parameters(parameters);
    EXPECT_FALSE(correction.sccs_is_active());
    EXPECT_TRUE(correction.uses_pcc());
    correction.v_correction_pcc(cell, basis, 1, density_channels, potential);
    EXPECT_NEAR(surchem::Ael, 2.0 * expected_energy, 1.0e-12);

}

} // namespace
