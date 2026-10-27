#include "../nonlocal_xc.h"
#include "source_base/constants.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_base/parallel_global.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <gtest/gtest.h>
#include <tuple>
#include <vector>

// The reference calls libvdwxc directly. Include its C header after the C++ headers.
#include <complex.h>
#include <complex>
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

namespace
{
int test_rank = 0;
int test_size = 1;
typedef std::tuple<double, double, ModuleBase::matrix> VxcResult;
// rVV10 parameters that set_xc_type assigns to RVV10
const double RVV10_B = 6.3;
const double RVV10_C = 0.0093;
// The adapter and the reference share the plane-wave gradient, so they agree to rounding.
const double ROUNDING_TOLERANCE = 2.0e-12;
// Central difference on a smooth density: the truncation error is of order step^2.
const double FINITE_DIFFERENCE_STEP = 1.0e-5;
const double FINITE_DIFFERENCE_TOLERANCE = 2.0e-6;

typedef double (*GridFunction)(const ModulePW::PW_Basis&, const int, const int, const int);

double valence_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ax = ModuleBase::TWO_PI * ix / pw.nx;
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    const double az = ModuleBase::TWO_PI * iz / pw.nz;
    return 0.004 + 0.006 * (1.0 + 0.25 * std::sin(ax) + 0.15 * std::cos(ay + az));
}

double core_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double az = ModuleBase::TWO_PI * iz / pw.nz;
    return 0.002 * (1.0 + std::cos(az));
}

double perturbation(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ax = ModuleBase::TWO_PI * ix / pw.nx;
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    return std::sin(ax + ay) + 0.3 * std::cos(2.0 * ay);
}

double x_only_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    return 0.01 + 0.003 * std::cos(ModuleBase::TWO_PI * ix / pw.nx);
}

double zero_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    return 0.0;
}

/// Values of f on the z planes of this process, in the ABACUS layout [ixy * nplane + iz].
std::vector<double> local_grid(const ModulePW::PW_Basis& pw, GridFunction f)
{
    std::vector<double> values(pw.nrxx);
    for (int ix = 0; ix < pw.nx; ++ix)
    {
        for (int iy = 0; iy < pw.ny; ++iy)
        {
            for (int iz = 0; iz < pw.nplane; ++iz)
            {
                values[(ix * pw.ny + iy) * pw.nplane + iz] = f(pw, ix, iy, pw.startz_current + iz);
            }
        }
    }
    return values;
}

/// A Charge on pw whose spin channels and core density point into the given vectors.
void attach(Charge& chr,
            ModulePW::PW_Basis& pw,
            std::vector<std::vector<double>>& channels,
            std::vector<double*>& rows,
            std::vector<double>& core)
{
    rows.clear();
    for (std::vector<double>& channel: channels)
    {
        rows.push_back(channel.data());
    }
    chr.rhopw = &pw;
    chr.nrxx = pw.nrxx;
    chr.nspin = static_cast<int>(channels.size());
    chr.rho = rows.data();
    chr.rho_core = core.data();
}

/// rVV10 contribution for the given valence channels and core density on pw.
VxcResult rvv10_result(ModulePW::PW_Basis& pw, std::vector<std::vector<double>> channels, std::vector<double> core)
{
    Charge chr;
    std::vector<double*> rows;
    attach(chr, pw, channels, rows, core);
    VxcResult result(0.0, 0.0, ModuleBase::matrix(chr.nspin, pw.nrxx));
    vdw::add_rvv10(RVV10_B, RVV10_C, chr, result);
    chr.rho = nullptr;
    chr.rho_core = nullptr;
    return result;
}

/// Oblique lattice; z_sign = -1 gives its left-handed mirror image.
ModuleBase::Matrix3 oblique_lattice(const double z_sign)
{
    const double c3 = 2.0 * z_sign;
    return ModuleBase::Matrix3(1.0, 0.0, 0.0, -0.5, std::sqrt(0.75), 0.0, 0.0, 0.0, c3);
}

class NonlocalXCTest : public testing::Test
{
  protected:
    /// 6x8x9 grid on an oblique cell, distributed over the processes of the test (5 + 4 planes on two)
    ModulePW::PW_Basis pw{"cpu", "double"};
    /// the same grid, complete on every process
    ModulePW::PW_Basis serial{"cpu", "double"};

    void SetUp() override
    {
        // The oblique cell catches row/column mistakes in the cell handed to libvdwxc.
        const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
        pw.initmpi(test_size, test_rank, MPI_COMM_WORLD);
        serial.initmpi(1, 0, MPI_COMM_SELF);
        setup(pw, lattice);
        setup(serial, lattice);
    }

    static void setup(ModulePW::PW_Basis& basis, const ModuleBase::Matrix3& lattice)
    {
        basis.initgrids(10.0, lattice, 6, 8, 9);
        basis.initparameters(false, 4.0, 1, true);
        basis.setuptransform();
        basis.collect_local_pw();
    }
};

TEST(NonlocalXCSelectionTest, ShortNameSelectsRvv10)
{
    XC_Functional::set_xc_type("rVV10");
    EXPECT_EQ(XC_Functional::get_func_type(), 2);
    const std::vector<int> expected_id = {XC_GGA_X_RPW86, XC_GGA_C_PBE};
    EXPECT_EQ(XC_Functional::get_func_id(), expected_id);
    EXPECT_DOUBLE_EQ(XC_Functional::get_nlc_b(), RVV10_B);
    {
        XC_Functional::ScopedType semilocal("PBE");
        EXPECT_EQ(XC_Functional::get_nlc_b(), 0.0);
    }
    EXPECT_DOUBLE_EQ(XC_Functional::get_nlc_b(), RVV10_B);
    XC_Functional::set_xc_type("PBE");
    EXPECT_EQ(XC_Functional::get_nlc_b(), 0.0);
}

TEST_F(NonlocalXCTest, MatchesDirectLibvdwxc)
{
    const std::vector<double> valence = local_grid(serial, valence_density);
    const std::vector<double> core = local_grid(serial, core_density);
    const VxcResult result = rvv10_result(serial, {valence}, core);

    // reference: libvdwxc on the total density with the plane-wave gradient, in Hartree units
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
    ASSERT_EQ(vdwxc_set_rvv10_parameters(direct, RVV10_B, RVV10_C), 0);
    const ModuleBase::Matrix3 a = serial.latvec * serial.lat0;
    vdwxc_set_unit_cell(direct,
                        serial.nx,
                        serial.ny,
                        serial.nz,
                        a.e11,
                        a.e12,
                        a.e13,
                        a.e21,
                        a.e22,
                        a.e23,
                        a.e31,
                        a.e32,
                        a.e33);
    vdwxc_init_serial(direct);
    const double energy_ha = vdwxc_calculate(direct, total.data(), sigma.data(), dedn.data(), dedsigma.data());
    vdwxc_finalize(&direct);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        grad[ir] *= 2.0 * dedsigma[ir];
    }
    std::vector<double> divergence(serial.nrxx);
    XC_Functional::grad_dot(grad.data(), divergence.data(), &serial, serial.tpiba);

    const double expected_energy = ModuleBase::e2 * energy_ha;
    EXPECT_NEAR(std::get<0>(result), expected_energy, ROUNDING_TOLERANCE);
    double vtxc = 0.0;
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        const double potential = ModuleBase::e2 * (dedn[ir] - divergence[ir]);
        EXPECT_NEAR(std::get<2>(result)(0, ir), potential, ROUNDING_TOLERANCE);
        vtxc += potential * valence[ir];
    }
    const double expected_vtxc = vtxc * serial.omega / serial.nxyz;
    EXPECT_NEAR(std::get<1>(result), expected_vtxc, ROUNDING_TOLERANCE);
}

TEST_F(NonlocalXCTest, SpinChannelsShareTheTotalDensity)
{
    const std::vector<double> valence = local_grid(serial, valence_density);
    const std::vector<double> core = local_grid(serial, core_density);
    std::vector<double> up(serial.nrxx);
    std::vector<double> down(serial.nrxx);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        up[ir] = 0.7 * valence[ir];
        down[ir] = valence[ir] - up[ir];
    }
    const VxcResult unpolarized = rvv10_result(serial, {valence}, core);
    const VxcResult polarized = rvv10_result(serial, {up, down}, core);
    EXPECT_NEAR(std::get<0>(polarized), std::get<0>(unpolarized), ROUNDING_TOLERANCE);
    EXPECT_NEAR(std::get<1>(polarized), std::get<1>(unpolarized), ROUNDING_TOLERANCE);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        EXPECT_NEAR(std::get<2>(polarized)(0, ir), std::get<2>(unpolarized)(0, ir), ROUNDING_TOLERANCE);
        EXPECT_NEAR(std::get<2>(polarized)(1, ir), std::get<2>(unpolarized)(0, ir), ROUNDING_TOLERANCE);
    }
}

TEST_F(NonlocalXCTest, DistributedGridMatchesCompleteGrid)
{
    const std::vector<double> valence = local_grid(pw, valence_density);
    const std::vector<double> core = local_grid(pw, core_density);
    const std::vector<double> valence_complete = local_grid(serial, valence_density);
    const std::vector<double> core_complete = local_grid(serial, core_density);
    const VxcResult distributed = rvv10_result(pw, {valence}, core);
    const VxcResult complete = rvv10_result(serial, {valence_complete}, core_complete);
    EXPECT_NEAR(std::get<0>(distributed), std::get<0>(complete), ROUNDING_TOLERANCE);
    EXPECT_NEAR(std::get<1>(distributed), std::get<1>(complete), ROUNDING_TOLERANCE);
    for (int ixy = 0; ixy < pw.nxy; ++ixy)
    {
        for (int iz = 0; iz < pw.nplane; ++iz)
        {
            const int local = ixy * pw.nplane + iz;
            const int global = ixy * pw.nz + pw.startz_current + iz;
            EXPECT_NEAR(std::get<2>(distributed)(0, local), std::get<2>(complete)(0, global), ROUNDING_TOLERANCE);
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
        plus[ir] += FINITE_DIFFERENCE_STEP * direction[ir];
        minus[ir] -= FINITE_DIFFERENCE_STEP * direction[ir];
    }
    const VxcResult reference = rvv10_result(pw, {valence}, core);
    const VxcResult result_plus = rvv10_result(pw, {plus}, core);
    const VxcResult result_minus = rvv10_result(pw, {minus}, core);
    const double energy_difference = std::get<0>(result_plus) - std::get<0>(result_minus);
    const double finite_difference = energy_difference / (2.0 * FINITE_DIFFERENCE_STEP);

    double contraction = 0.0;
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        contraction += std::get<2>(reference)(0, ir) * direction[ir];
    }
    Parallel_Common::reduce_data(&contraction, 1, pw.pool_world);
    contraction *= pw.omega / pw.nxyz;
    EXPECT_NEAR(finite_difference, contraction, FINITE_DIFFERENCE_TOLERANCE);
}

TEST_F(NonlocalXCTest, ObliqueCellMatchesEquivalentOrthogonalCell)
{
    // Both cells have the volume sqrt(3) lat0^3 and the plane spacing sqrt(3)/2 lat0 along the
    // first lattice direction, so a density that depends only on the first fractional coordinate
    // has the same nonlocal energy and potential in both. A transposed cell matrix breaks this.
    ModulePW::PW_Basis orthogonal{"cpu", "double"};
    orthogonal.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 rectangle(std::sqrt(0.75), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0);
    setup(orthogonal, rectangle);
    ASSERT_NEAR(serial.omega, orthogonal.omega, ROUNDING_TOLERANCE);

    const std::vector<double> density = local_grid(serial, x_only_density);
    const std::vector<double> no_core = local_grid(serial, zero_density);
    const VxcResult oblique = rvv10_result(serial, {density}, no_core);
    const VxcResult rectangular = rvv10_result(orthogonal, {density}, no_core);
    EXPECT_NEAR(std::get<0>(oblique), std::get<0>(rectangular), ROUNDING_TOLERANCE);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        EXPECT_NEAR(std::get<2>(oblique)(0, ir), std::get<2>(rectangular)(0, ir), ROUNDING_TOLERANCE);
    }
}

TEST_F(NonlocalXCTest, LeftHandedCellGivesTheMirrorImageResult)
{
    // A mirror image has the same nonlocal energy and potential; libvdwxc's signed cell volume would
    // otherwise flip the sign of the energy.
    ModulePW::PW_Basis mirror{"cpu", "double"};
    mirror.initmpi(1, 0, MPI_COMM_SELF);
    setup(mirror, oblique_lattice(-1.0));
    ASSERT_LT(mirror.latvec.Det(), 0.0);

    const std::vector<double> valence = local_grid(serial, valence_density);
    const std::vector<double> core = local_grid(serial, core_density);
    const VxcResult right_handed = rvv10_result(serial, {valence}, core);
    const VxcResult left_handed = rvv10_result(mirror, {valence}, core);
    EXPECT_NEAR(std::get<0>(left_handed), std::get<0>(right_handed), ROUNDING_TOLERANCE);
    EXPECT_NEAR(std::get<1>(left_handed), std::get<1>(right_handed), ROUNDING_TOLERANCE);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        EXPECT_NEAR(std::get<2>(left_handed)(0, ir), std::get<2>(right_handed)(0, ir), ROUNDING_TOLERANCE);
    }
}

TEST_F(NonlocalXCTest, VxcOfRvv10AddsTheNonlocalCorrelation)
{
    // XC_Functional::v_xc for dft_functional RVV10 equals the semilocal rPW86+PBE result plus add_rvv10.
    const std::vector<double> valence = local_grid(pw, valence_density);
    std::vector<double> core = local_grid(pw, core_density);
    std::vector<std::vector<double>> channels(2, std::vector<double>(pw.nrxx));
    for (int ir = 0; ir < pw.nrxx; ++ir)
    {
        channels[0][ir] = 0.7 * valence[ir];
        channels[1][ir] = valence[ir] - channels[0][ir];
    }
    Charge chr;
    std::vector<double*> rows;
    attach(chr, pw, channels, rows, core);
    UnitCell ucell;
    ucell.omega = pw.omega;
    ucell.tpiba = pw.tpiba;
    const int nspin = 2;
    const bool domag = false;
    const int gga_grad = 0;
    const double hybrid_alpha = 0.0;
    const double hse_omega = 0.0;

    XC_Functional::set_xc_type("RVV10");
    const VxcResult total
        = XC_Functional::v_xc(pw.nrxx, &chr, &ucell, nspin, domag, domag, gga_grad, hybrid_alpha, hse_omega);
    VxcResult expected;
    {
        XC_Functional::ScopedType semilocal("GGA_X_RPW86+GGA_C_PBE");
        expected = XC_Functional::v_xc(pw.nrxx, &chr, &ucell, nspin, domag, domag, gga_grad, hybrid_alpha, hse_omega);
    }
    vdw::add_rvv10(RVV10_B, RVV10_C, chr, expected);
    XC_Functional::set_xc_type("PBE");
    chr.rho = nullptr;
    chr.rho_core = nullptr;

    EXPECT_NEAR(std::get<0>(total), std::get<0>(expected), ROUNDING_TOLERANCE);
    EXPECT_NEAR(std::get<1>(total), std::get<1>(expected), ROUNDING_TOLERANCE);
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < pw.nrxx; ++ir)
        {
            EXPECT_NEAR(std::get<2>(total)(is, ir), std::get<2>(expected)(is, ir), ROUNDING_TOLERANCE);
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    int threads = 1;
    Parallel_Global::read_pal_param(argc, argv, test_size, threads, test_rank);
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
