#include "source_base/constants.h"
#include "source_estate/module_pot/pot_nonlocal_xc.h"
#include "source_hamilt/module_vdw/nonlocal_xc_libvdwxc.h"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

TEST(PotNonlocalXCTest, AddsCorrectionFromValenceAndCoreWithoutReplacingExistingXc)
{
    ModulePW::PW_Basis pw("cpu", "double");
    const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0, -0.5, std::sqrt(0.75), 0.0, 0.0, 0.0, 2.0);
    pw.initgrids(10.0, lattice, 6, 8, 12);
    pw.initparameters(false, 4.0, 1, true);
    pw.setuptransform();
    pw.collect_local_pw();

    UnitCell cell;
    cell.lat0 = pw.lat0;
    cell.latvec = lattice;
    cell.omega = pw.omega;
    Charge charge;
    charge.set_rhopw(&pw);
    charge.allocate(1, false, false, 0);

    // Host and adapter add the same numbers, so they agree to rounding.
    const double rounding_tolerance = 2.0e-12;
    std::vector<double> valence(pw.nrxx);
    std::vector<double> core(pw.nrxx);
    std::vector<double> total(pw.nrxx);
    ModuleBase::matrix potential(1, pw.nrxx);
    std::vector<double> initial_potential(pw.nrxx);
    for (int x = 0; x < pw.nx; ++x)
    {
        for (int y = 0; y < pw.ny; ++y)
        {
            for (int z = 0; z < pw.nz; ++z)
            {
                const int r = (x * pw.ny + y) * pw.nz + z;
                const double ax = ModuleBase::TWO_PI * x / pw.nx;
                const double ay = ModuleBase::TWO_PI * y / pw.ny;
                const double az = ModuleBase::TWO_PI * z / pw.nz;
                valence[r] = 0.007 + 0.001 * std::cos(ax) + 0.0005 * std::sin(ay + az);
                core[r] = 0.003 + 0.0008 * std::cos(ay);
                total[r] = valence[r] + core[r];
                charge.rho[0][r] = valence[r];
                charge.rho_core[r] = core[r];
                initial_potential[r] = 0.25 + 0.0001 * r;
                potential(0, r) = initial_potential[r];
            }
        }
    }

    // The separately tested adapter supplies the correction; this test checks
    // the production host's density construction and additive XC bookkeeping.
    vdw::LibvdwxcNonlocalXC reference(pw, vdw::RVV10_B_RPW86_PBE, vdw::RVV10_C_RPW86_PBE);
    const vdw::NonlocalXCResult correction = reference.evaluate(pw, total, valence);
    const double initial_energy = -1.25;
    const double initial_vtxc = 0.375;
    double energy = initial_energy;
    double vtxc = initial_vtxc;
    elecstate::PotNonlocalXC host(&pw, &energy, &vtxc, "rvv10");
    host.cal_v_eff(&charge, &cell, potential);

    const double expected_energy = initial_energy + correction.energy;
    const double expected_vtxc = initial_vtxc + correction.vtxc;
    EXPECT_NEAR(energy, expected_energy, rounding_tolerance);
    EXPECT_NEAR(vtxc, expected_vtxc, rounding_tolerance);
    for (int r = 0; r < pw.nrxx; ++r)
    {
        const double expected_potential = initial_potential[r] + correction.potential[r];
        EXPECT_NEAR(potential(0, r), expected_potential, rounding_tolerance);
        EXPECT_DOUBLE_EQ(charge.rho[0][r], valence[r]);
        EXPECT_DOUBLE_EQ(charge.rho_core[r], core[r]);
    }
}
