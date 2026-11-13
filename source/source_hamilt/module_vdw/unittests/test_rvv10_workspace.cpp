// The libvdwxc workspace (rvv10_workspace.cpp): redistribution between the ABACUS planes and the libvdwxc slabs,
// the cell orientation, and the rebuild after a change of grid, cell, processes or a release.
#include "nonlocal_xc_test_helpers.h"

#include "source_basis/module_pw/pw_basis_big.h"

#include <cmath>

namespace rvv10_test
{
TEST_F(NonlocalXCTest, DistributedGridMatchesCompleteGrid)
{
    expect_distributed_matches_complete(pw, serial);
}

TEST_F(NonlocalXCTest, BigGridDistributionMatchesCompleteGrid)
{
    // The density grid of a run is a PW_Basis_Big, which hands out whole blocks of bz planes: on two
    // processes 8 + 4 against libvdwxc's 6 + 6, on four 4 + 4 + 4 + 0 against 3 + 3 + 3 + 3, so that a
    // process may own no plane at all and one ABACUS block may span two libvdwxc slabs.
    ModulePW::PW_Basis_Big big{"cpu", "double"};
    big.setbxyz(1, 1, 4);
    big.initmpi(test_size, test_rank, MPI_COMM_WORLD);
    ModulePW::PW_Basis complete{"cpu", "double"};
    complete.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
    setup_grid(big, lattice, basis_cutoff, 5, 7, 12);
    setup_grid(complete, lattice, basis_cutoff, 5, 7, 12);
    ASSERT_EQ(big.nz, 12);
    ASSERT_EQ(big.nplane % 4, 0);
    if (test_size == 4)
    {
        EXPECT_EQ(big.numz[3], 0);
    }
    expect_distributed_matches_complete(big, complete);
    // the XC evaluator calls add_rvv10 on every process of the pool, also on one without points; v_xc_libxc
    // hands its gradient over, so its result is identical to the last bit
    UnitCell ucell;
    ucell.omega = big.omega;
    ucell.tpiba = big.tpiba;
    ucell.magnet.lsign_ = false;
    const std::vector<std::vector<double>> unpolarized(1, local_grid(big, valence_density));
    const ChargeOnGrid single(big, unpolarized, local_grid(big, core_density));
#ifdef __LIBXC
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", single.chr, ucell, false, 0, 0.0);
#else
    expect_vxc_adds_rvv10("PBE", single.chr, ucell, false, 0, rounding_tolerance);
#endif
}

TEST_F(NonlocalXCTest, LeftHandedCellGivesTheMirrorImageResult)
{
    // A mirror image has the same result; libvdwxc's signed cell volume would flip the sign of the energy.
    ModulePW::PW_Basis mirror{"cpu", "double"};
    mirror.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 mirrored = oblique_lattice(-1.0);
    setup(mirror, mirrored, basis_cutoff);
    ASSERT_LT(mirror.latvec.Det(), 0.0);
    const std::vector<double> valence = local_grid(serial, valence_density);
    const std::vector<double> core = local_grid(serial, core_density);
    const VxcResult right_handed = rvv10_result(serial, {valence}, core);
    const VxcResult left_handed = rvv10_result(mirror, {valence}, core);
    expect_same_result(left_handed, right_handed, rounding_tolerance, rounding_tolerance);
}

TEST_F(NonlocalXCTest, WorkspaceFollowsGridCellAndProcessChanges)
{
    // The libvdwxc object is reused between calls; switching the process layout, the cell, the grid and
    // back must reproduce the earlier results.
    ModulePW::PW_Basis orthogonal{"cpu", "double"};
    orthogonal.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 rectangle(std::sqrt(0.75), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0);
    setup(orthogonal, rectangle, basis_cutoff);
    ModulePW::PW_Basis finer{"cpu", "double"};
    finer.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
    setup_grid(finer, lattice, basis_cutoff, 5, 7, 11);
    const std::vector<double> valence_finer = local_grid(finer, valence_density);
    const std::vector<double> core_finer = local_grid(finer, core_density);
    const std::vector<double> valence = local_grid(serial, valence_density);
    const std::vector<double> core = local_grid(serial, core_density);
    const std::vector<double> valence_local = local_grid(pw, valence_density);
    const std::vector<double> core_local = local_grid(pw, core_density);
    const VxcResult first = rvv10_result(serial, {valence}, core);
    const VxcResult repeated = rvv10_result(serial, {valence}, core);
    const VxcResult distributed = rvv10_result(pw, {valence_local}, core_local);
    const VxcResult other_cell = rvv10_result(orthogonal, {valence}, core);
    const VxcResult other_grid = rvv10_result(finer, {valence_finer}, core_finer);
    const VxcResult again = rvv10_result(serial, {valence}, core);
    const VxcResult other_grid_again = rvv10_result(finer, {valence_finer}, core_finer);
    expect_same_result(repeated, first, 0.0, 0.0);
    expect_same_result(again, first, rounding_tolerance, rounding_tolerance);
    expect_same_result(other_grid_again, other_grid, rounding_tolerance, rounding_tolerance);
    EXPECT_NEAR(std::get<0>(distributed), std::get<0>(first), rounding_tolerance);
    EXPECT_GT(std::abs(std::get<0>(other_cell) - std::get<0>(first)), 1.0e-8);
    // the end of a run releases the libvdwxc object; a later evaluation builds it again
    vdw::finalize_rvv10();
    const VxcResult after_release = rvv10_result(serial, {valence}, core);
    expect_same_result(after_release, first, 0.0, 0.0);
    vdw::finalize_rvv10();
}

} // namespace rvv10_test
