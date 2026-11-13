// rVV10 adapter (nonlocal_xc.cpp): libvdwxc called directly as the reference, the potential and the stress as
// derivatives of the energy, and every XC evaluator that adds the term. main() of the test target is here.
#include "nonlocal_xc_test_helpers.h"

#include "source_base/parallel_comm.h"
#include "source_base/parallel_global.h"
#ifdef __LIBXC
#include "source_hamilt/module_xc/libxc_abacus.h"
#endif

#include <cmath>

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

TEST_F(NonlocalXCTest, PotentialIsTheDerivativeOfTheEnergy)
{
    const std::vector<double> valence = local_grid(pw, valence_density);
    const std::vector<double> core = local_grid(pw, core_density);
    const std::vector<double> direction = local_grid(pw, perturbation);
    std::vector<double> plus = valence;
    std::vector<double> minus = valence;
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        plus[ir] += finite_difference_step * direction[ir];
        minus[ir] -= finite_difference_step * direction[ir];
    }
    const VxcResult reference = rvv10_result(pw, {valence}, core);
    const VxcResult result_plus = rvv10_result(pw, {plus}, core);
    const VxcResult result_minus = rvv10_result(pw, {minus}, core);
    const double energy_difference = std::get<0>(result_plus) - std::get<0>(result_minus);
    const double finite_difference = energy_difference / (2.0 * finite_difference_step);

    double contraction = 0.0;
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        contraction += std::get<2>(reference)(0, ir) * direction[ir];
    }
    Parallel_Common::reduce_data(&contraction, 1, pw.pool_world);
    contraction *= pw.omega / pw.nxyz;
    EXPECT_NEAR(finite_difference, contraction, finite_difference_tolerance);
}

TEST_F(NonlocalXCTest, StressIsTheStrainDerivativeOfTheEnergy)
{
    // A strain keeps the grid and every G component of the density, so the stress is the exact derivative
    // of the discrete energy: sigma_ab = -(1/omega) dE/d(eps_ab) at a fixed electron number.
    ModulePW::PW_Basis box{"cpu", "double"};
    box.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
    setup(box, lattice, largest_sphere_cutoff);
    const std::vector<double> density = local_grid(box, valence_density);
    const std::vector<double> no_core(box.nrxx, 0.0);
    const VxcResult reference = rvv10_result(box, {density}, no_core);
    const ModuleBase::matrix extra = rvv10_stress(box, {density}, no_core);
    // the XC diagonal term that ABACUS adds from etxc and vtxc
    const double diagonal = -(std::get<0>(reference) - std::get<1>(reference)) / box.omega;

    const int components[4][2] = {{0, 0}, {2, 2}, {0, 1}, {1, 2}};
    for (const int* ab: components)
    {
        const int a = ab[0];
        const int b = ab[1];
        double energy[2] = {0.0, 0.0};
        for (int side = 0; side < 2; ++side)
        {
            const double h = (side == 0) ? strain_step : -strain_step;
            ModulePW::PW_Basis strained{"cpu", "double"};
            strained.initmpi(1, 0, MPI_COMM_SELF);
            const ModuleBase::Matrix3 strained_cell = strained_lattice(lattice, a, b, h);
            setup(strained, strained_cell, largest_sphere_cutoff);
            // keep the electron number: the density scales with the inverse volume ratio
            const double volume_ratio = strained.omega / box.omega;
            std::vector<double> scaled = density;
            for (double& value: scaled)
            {
                value /= volume_ratio;
            }
            const VxcResult strained_result = rvv10_result(strained, {scaled}, no_core);
            energy[side] = std::get<0>(strained_result);
        }
        const double finite_difference = -(energy[0] - energy[1]) / (2.0 * strain_step * box.omega);
        const double analytic = extra(a, b) + ((a == b) ? diagonal : 0.0);
        const double tolerance = strain_tolerance * std::abs(analytic) + 1.0e-12;
        SCOPED_TRACE(testing::Message() << "component " << a << b);
        EXPECT_NEAR(analytic, finite_difference, tolerance);
    }
}

TEST_F(NonlocalXCTest, VxcAddsRvv10ToEveryFunctionalPath)
{
    // v_xc and v_xc_meta with rVV10 equal the semilocal result plus add_rvv10, for native PBE (final
    // return), noncollinear PBE with gga_grad 0 (collinear-mapped branch) and 2 (NCGGA return), and with
    // Libxc a libxc GGA (v_xc_libxc) in the same layouts and SCAN (v_xc_meta).
    const std::vector<double> valence = local_grid(pw, valence_density);
    const std::vector<double> core = local_grid(pw, core_density);
    const std::vector<double> zeros(pw.nrxx);
    std::vector<std::vector<double>> collinear(2, zeros);
    std::vector<std::vector<double>> noncollinear(4, zeros);
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        collinear[0][ir] = 0.7 * valence[ir];
        collinear[1][ir] = valence[ir] - collinear[0][ir];
        noncollinear[0][ir] = valence[ir];
        noncollinear[1][ir] = 0.1 * valence[ir];
        noncollinear[2][ir] = -0.05 * valence[ir];
        noncollinear[3][ir] = 0.2 * valence[ir];
    }
    UnitCell ucell;
    ucell.omega = pw.omega;
    ucell.tpiba = pw.tpiba;
    ucell.magnet.lsign_ = false;
    const ChargeOnGrid polarized(pw, collinear, core);
    const ChargeOnGrid magnetized(pw, noncollinear, core);
    expect_vxc_adds_rvv10("PBE", polarized.chr, ucell, false, 0, rounding_tolerance);
    expect_vxc_adds_rvv10("PBE", magnetized.chr, ucell, true, 2, rounding_tolerance);
    expect_vxc_adds_rvv10("PBE", magnetized.chr, ucell, true, 0, rounding_tolerance);
#ifdef __LIBXC
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", polarized.chr, ucell, false, 0, rounding_tolerance);
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", magnetized.chr, ucell, true, 0, rounding_tolerance);
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", magnetized.chr, ucell, true, 2, rounding_tolerance);
    // with one spin channel, v_xc_libxc and v_xc_meta hand their density gradient to add_rvv10, so the
    // results are identical to the last bit
    const std::vector<std::vector<double>> unpolarized(1, valence);
    const ChargeOnGrid single(pw, unpolarized, core);
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", single.chr, ucell, false, 0, 0.0);

    XC_Functional::ScopedType scan("SCAN");
    std::vector<double> kinetic = local_grid(pw, kinetic_density);
    std::vector<double*> kinetic_rows(1, kinetic.data());
    ChargeOnGrid charge(pw, unpolarized, core);
    charge.chr.kin_r = kinetic_rows.data();
    const std::vector<int> func_id = XC_Functional::get_func_id();
    const int nspin = 1;
    const double hybrid_alpha = 0.0;
    const double hse_omega = 0.0;
    VxcResult results[2];
    for (int with_rvv10 = 0; with_rvv10 < 2; ++with_rvv10)
    {
        const bool rvv10 = with_rvv10 == 1;
        RuntimeParameters parameters(rvv10);
        const MetaResult meta = XC_Functional_Libxc::v_xc_meta(func_id,
                                                               pw.nrxx,
                                                               pw.omega,
                                                               pw.tpiba,
                                                               &charge.chr,
                                                               nspin,
                                                               hybrid_alpha,
                                                               hse_omega);
        results[with_rvv10] = VxcResult(std::get<0>(meta), std::get<1>(meta), std::get<2>(meta));
    }
    expect_rvv10_added(results[1], results[0], charge.chr, 0.0);
#endif
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
