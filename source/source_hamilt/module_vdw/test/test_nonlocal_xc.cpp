#include "test_nonlocal_xc.h"

#include "source_base/parallel_comm.h"
#include "source_base/parallel_global.h"
#include "source_basis/module_pw/pw_basis_big.h"

// The reference calls libvdwxc directly. Include its C header after the C++ headers.
#include <complex.h>
#include <vdwxc.h>
#undef I

// The target links PW/XC objects only; the tests own the charge storage.
Charge::Charge()
{
}
Charge::~Charge()
{
}
UnitCell::UnitCell()
{
}
UnitCell::~UnitCell()
{
}
Magnetism::Magnetism()
{
}
Magnetism::~Magnetism()
{
}
SepPot::SepPot()
{
}
SepPot::~SepPot()
{
}
Sep_Cell::Sep_Cell() noexcept
{
}
Sep_Cell::~Sep_Cell() noexcept
{
}

namespace rvv10_test
{
int test_rank = 0;
int test_size = 1;

TEST_F(NonlocalXCTest, MatchesDirectLibvdwxc)
{
    // unpolarized and spin-polarized layouts against libvdwxc called directly in the ABACUS axis order
    const std::vector<double> valence = local_grid(serial, valence_density);
    const std::vector<double> core = local_grid(serial, core_density);
    const VxcResult result = rvv10_result(serial, {valence}, core);

    // reference: libvdwxc on the total density in the ABACUS axis order, Hartree units
    std::vector<double> total(serial.nrxx);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        total[ir] = valence[ir] + core[ir];
    }
    std::vector<std::complex<double>> rhog(serial.npw);
    serial.real2recip(total.data(), rhog.data());
    std::vector<ModuleBase::Vector3<double>> grad(serial.nrxx);
    XC_Functional::grad_rho(rhog.data(), grad.data(), &serial, serial.tpiba);
    std::vector<double> sigma(serial.nrxx);
    std::vector<double> dedn(serial.nrxx);
    std::vector<double> dedsigma(serial.nrxx);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        sigma[ir] = grad[ir].norm2();
    }
    vdwxc_data direct = vdwxc_new(VDWXC_RVV10);
    ASSERT_NE(direct, nullptr);
    ASSERT_EQ(vdwxc_set_rvv10_parameters(direct, rvv10_b_rpw86, rvv10_c_rpw86), 0);
    const ModuleBase::Matrix3 a = serial.latvec * serial.lat0;
    const int nx = serial.nx;
    const int ny = serial.ny;
    const int nz = serial.nz;
    vdwxc_set_unit_cell(direct, nx, ny, nz, a.e11, a.e12, a.e13, a.e21, a.e22, a.e23, a.e31, a.e32, a.e33);
    vdwxc_init_serial(direct);
    const double energy_ha = vdwxc_calculate(direct, total.data(), sigma.data(), dedn.data(), dedsigma.data());
    vdwxc_finalize(&direct);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        grad[ir] *= 2.0 * dedsigma[ir];
    }
    std::vector<double> divergence(serial.nrxx);
    XC_Functional::grad_dot(grad.data(), divergence.data(), &serial, serial.tpiba);

    ModuleBase::matrix v(1, serial.nrxx);
    double vtxc = 0.0;
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        const double potential = ModuleBase::e2 * (dedn[ir] - divergence[ir]);
        v(0, ir) = potential;
        vtxc += potential * valence[ir];
    }
    const double etxc = ModuleBase::e2 * energy_ha;
    vtxc *= serial.omega / serial.nxyz;
    const VxcResult expected(etxc, vtxc, v);
    expect_same_result(result, expected, rounding_tolerance, rounding_tolerance);

    // spin channels share the total density: same energy and vtxc, and each channel gets the potential
    std::vector<double> up(serial.nrxx);
    std::vector<double> down(serial.nrxx);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        up[ir] = 0.7 * valence[ir];
        down[ir] = valence[ir] - up[ir];
    }
    const VxcResult polarized = rvv10_result(serial, {up, down}, core);
    EXPECT_NEAR(std::get<0>(polarized), etxc, rounding_tolerance);
    EXPECT_NEAR(std::get<1>(polarized), vtxc, rounding_tolerance);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        EXPECT_NEAR(std::get<2>(polarized)(0, ir), v(0, ir), rounding_tolerance);
        EXPECT_NEAR(std::get<2>(polarized)(1, ir), v(0, ir), rounding_tolerance);
    }

    // four-component layout: rho[0] is the total valence density, and only v(0) gets the potential
    std::vector<double> mx(serial.nrxx);
    std::vector<double> my(serial.nrxx);
    std::vector<double> mz(serial.nrxx);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        mx[ir] = 0.1 * valence[ir];
        my[ir] = -0.05 * valence[ir];
        mz[ir] = 0.2 * valence[ir];
    }
    const VxcResult noncollinear = rvv10_result(serial, {valence, mx, my, mz}, core);
    EXPECT_NEAR(std::get<0>(noncollinear), etxc, rounding_tolerance);
    EXPECT_NEAR(std::get<1>(noncollinear), vtxc, rounding_tolerance);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        EXPECT_NEAR(std::get<2>(noncollinear)(0, ir), v(0, ir), rounding_tolerance);
        for (int is = 1; is < 4; ++is)
        {
            EXPECT_EQ(std::get<2>(noncollinear)(is, ir), 0.0);
        }
    }
}

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
    // v_xc_libxc hands its gradient to add_rvv10 on every process of the pool, also on one without points
    UnitCell ucell;
    ucell.omega = big.omega;
    ucell.tpiba = big.tpiba;
    ucell.magnet.lsign_ = false;
    const std::vector<std::vector<double>> unpolarized(1, local_grid(big, valence_density));
    const ChargeOnGrid single(big, unpolarized, local_grid(big, core_density));
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", single.chr, ucell, false, 0, 0.0);
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

int main(int argc, char** argv)
{
    int threads = 1;
    Parallel_Global::read_pal_param(argc, argv, rvv10_test::test_size, threads, rvv10_test::test_rank);
    POOL_WORLD = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    Parallel_Global::finalize_mpi();
    return result;
}
