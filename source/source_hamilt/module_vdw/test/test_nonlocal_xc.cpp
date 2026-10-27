#include "../nonlocal_xc.h"
#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_base/parallel_global.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_xc/libxc_abacus.h"
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
// rVV10 parameters of the rPW86-PBE pairing
const double RVV10_B = 6.3;
const double RVV10_C = 0.0093;
// Results that share the plane-wave gradient agree to rounding.
const double ROUNDING_TOLERANCE = 2.0e-12;
// Central differences on smooth densities: the truncation error is of order step^2.
const double FINITE_DIFFERENCE_STEP = 1.0e-5;
const double FINITE_DIFFERENCE_TOLERANCE = 2.0e-6;
const double STRAIN_STEP = 1.0e-4;
const double STRAIN_TOLERANCE = 1.0e-5;
const double BASIS_CUTOFF = 4.0;
// PW_Basis lowers a larger cutoff to the largest G sphere inside the FFT box. The test densities have
// only first-shell G components, well inside that sphere, so a small strain changes none of them.
const double LARGEST_SPHERE_CUTOFF = 1.0e6;

typedef double (*GridFunction)(const ModulePW::PW_Basis&, const int, const int, const int);

double valence_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ax = ModuleBase::TWO_PI * ix / pw.nx;
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    const double az = ModuleBase::TWO_PI * iz / pw.nz;
    const double ayz = ay + az;
    return 0.004 + 0.006 * (1.0 + 0.25 * std::sin(ax) + 0.15 * std::cos(ayz));
}

double core_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double az = ModuleBase::TWO_PI * iz / pw.nz;
    return 0.002 * (1.0 + std::cos(az));
}

double kinetic_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    return 0.003 * (1.0 + 0.2 * std::cos(ay));
}

double perturbation(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ax = ModuleBase::TWO_PI * ix / pw.nx;
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    const double axy = ax + ay;
    const double ay2 = 2.0 * ay;
    return std::sin(axy) + 0.3 * std::cos(ay2);
}

double x_only_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ax = ModuleBase::TWO_PI * ix / pw.nx;
    return 0.01 + 0.003 * std::cos(ax);
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

/// Oblique lattice; z_sign = -1 gives its left-handed mirror image.
ModuleBase::Matrix3 oblique_lattice(const double z_sign)
{
    const double c3 = 2.0 * z_sign;
    return ModuleBase::Matrix3(1.0, 0.0, 0.0, -0.5, std::sqrt(0.75), 0.0, 0.0, 0.0, c3);
}

/// The lattice after the strain eps_ab = h: every lattice vector gets a_i[a] += h * a_i[b].
ModuleBase::Matrix3 strained_lattice(const ModuleBase::Matrix3& lattice, const int a, const int b, const double h)
{
    double m[3][3] = {{lattice.e11, lattice.e12, lattice.e13},
                      {lattice.e21, lattice.e22, lattice.e23},
                      {lattice.e31, lattice.e32, lattice.e33}};
    for (int i = 0; i < 3; ++i)
    {
        m[i][a] += h * m[i][b];
    }
    return ModuleBase::Matrix3(m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]);
}

/// Basis on the nx x ny x nz FFT grid of the given cell.
void setup_grid(ModulePW::PW_Basis& basis,
                const ModuleBase::Matrix3& lattice,
                const double cutoff,
                const int nx,
                const int ny,
                const int nz)
{
    basis.initgrids(10.0, lattice, nx, ny, nz);
    basis.initparameters(false, cutoff, 1, true);
    basis.setuptransform();
    basis.collect_local_pw();
}

/// 5x7x9 grid: on two processes the z planes split 5 + 4, on four processes 3 + 2 + 2 + 2 against
/// libvdwxc's 3 + 3 + 3 + 0, which exercises the redistribution and an empty slab. Odd nx and ny leave
/// no Nyquist plane, so libvdwxc's result does not depend on the order in which it gets the axes.
void setup(ModulePW::PW_Basis& basis, const ModuleBase::Matrix3& lattice, const double cutoff)
{
    setup_grid(basis, lattice, cutoff, 5, 7, 9);
}

/// Sets the rVV10 parameters of the XC runtime parameters for the lifetime of the object.
class Rvv10Parameters
{
  public:
    Rvv10Parameters(const double b, const double c)
    {
        saved_ = XC_Functional::get_runtime_parameters();
        XCFunctionalParameters parameters = saved_;
        parameters.rvv10_b = b;
        parameters.rvv10_c = c;
        XC_Functional::set_runtime_parameters(parameters);
    }
    ~Rvv10Parameters()
    {
        XC_Functional::set_runtime_parameters(saved_);
    }

  private:
    XCFunctionalParameters saved_;
};

/// A Charge on pw that owns its density components, core density and the G-space buffers of gradcorr.
class ChargeOnGrid
{
  public:
    ChargeOnGrid(ModulePW::PW_Basis& pw,
                 const std::vector<std::vector<double>>& channels,
                 const std::vector<double>& core)
    {
        channels_ = channels;
        core_ = core;
        const std::size_t ncomponents = channels.size();
        const std::vector<std::complex<double>> reciprocal(pw.npw);
        rhog_.assign(ncomponents, reciprocal);
        core_g_ = reciprocal;
        for (std::size_t is = 0; is < ncomponents; ++is)
        {
            rows_.push_back(channels_[is].data());
            rhog_rows_.push_back(rhog_[is].data());
        }
        chr.rhopw = &pw;
        chr.nrxx = pw.nrxx;
        chr.nspin = static_cast<int>(ncomponents);
        chr.rho = rows_.data();
        chr.rhog = rhog_rows_.data();
        chr.rho_core = core_.data();
        chr.rhog_core = core_g_.data();
    }
    ~ChargeOnGrid()
    {
        chr.rho = nullptr;
        chr.rhog = nullptr;
        chr.rho_core = nullptr;
        chr.rhog_core = nullptr;
        chr.kin_r = nullptr;
    }
    ChargeOnGrid(const ChargeOnGrid&) = delete;
    ChargeOnGrid& operator=(const ChargeOnGrid&) = delete;

    Charge chr;

  private:
    std::vector<std::vector<double>> channels_;
    std::vector<double*> rows_;
    std::vector<double> core_;
    std::vector<std::vector<std::complex<double>>> rhog_;
    std::vector<std::complex<double>*> rhog_rows_;
    std::vector<std::complex<double>> core_g_;
};

/// rVV10 contribution for the given valence channels and core density on pw.
VxcResult rvv10_result(ModulePW::PW_Basis& pw,
                       const std::vector<std::vector<double>>& channels,
                       const std::vector<double>& core)
{
    const ChargeOnGrid charge(pw, channels, core);
    double etxc = 0.0;
    double vtxc = 0.0;
    ModuleBase::matrix v(charge.chr.nspin, pw.nrxx);
    vdw::add_rvv10(charge.chr, etxc, vtxc, v);
    return VxcResult(etxc, vtxc, v);
}

/// rVV10 stress beyond the XC diagonal term for the given density on pw.
ModuleBase::matrix rvv10_stress(ModulePW::PW_Basis& pw,
                                const std::vector<std::vector<double>>& channels,
                                const std::vector<double>& core)
{
    const ChargeOnGrid charge(pw, channels, core);
    ModuleBase::matrix sigma(3, 3);
    vdw::add_rvv10_stress(charge.chr, sigma);
    return sigma;
}

void expect_same_result(const VxcResult& actual, const VxcResult& expected, const double tolerance)
{
    EXPECT_NEAR(std::get<0>(actual), std::get<0>(expected), tolerance);
    EXPECT_NEAR(std::get<1>(actual), std::get<1>(expected), tolerance);
    const ModuleBase::matrix& v = std::get<2>(actual);
    const ModuleBase::matrix& w = std::get<2>(expected);
    ASSERT_EQ(v.nr, w.nr);
    ASSERT_EQ(v.nc, w.nc);
    for (int is = 0; is < v.nr; ++is)
    {
        for (int ir = 0; ir < v.nc; ++ir)
        {
            EXPECT_NEAR(v(is, ir), w(is, ir), tolerance);
        }
    }
}

/// A functional evaluated with rVV10 equals the evaluation without it plus add_rvv10.
void expect_rvv10_added(const VxcResult& total, VxcResult semilocal, const Charge& chr)
{
    double& etxc = std::get<0>(semilocal);
    double& vtxc = std::get<1>(semilocal);
    ModuleBase::matrix& v = std::get<2>(semilocal);
    vdw::add_rvv10(chr, etxc, vtxc, v);
    expect_same_result(total, semilocal, ROUNDING_TOLERANCE);
}

class NonlocalXCTest : public testing::Test
{
  protected:
    Rvv10Parameters rvv10{RVV10_B, RVV10_C};
    /// grid on an oblique cell, distributed over the processes of the test
    ModulePW::PW_Basis pw{"cpu", "double"};
    /// the same grid, complete on every process
    ModulePW::PW_Basis serial{"cpu", "double"};

    void SetUp() override
    {
        // The oblique cell catches row/column mistakes in the cell handed to libvdwxc.
        const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
        pw.initmpi(test_size, test_rank, MPI_COMM_WORLD);
        serial.initmpi(1, 0, MPI_COMM_SELF);
        setup(pw, lattice, BASIS_CUTOFF);
        setup(serial, lattice, BASIS_CUTOFF);
    }
};

TEST_F(NonlocalXCTest, MatchesDirectLibvdwxc)
{
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
    expect_same_result(result, expected, ROUNDING_TOLERANCE);
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
    setup(orthogonal, rectangle, BASIS_CUTOFF);
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
    const ModuleBase::Matrix3 mirrored = oblique_lattice(-1.0);
    setup(mirror, mirrored, BASIS_CUTOFF);
    ASSERT_LT(mirror.latvec.Det(), 0.0);

    const std::vector<double> valence = local_grid(serial, valence_density);
    const std::vector<double> core = local_grid(serial, core_density);
    const VxcResult right_handed = rvv10_result(serial, {valence}, core);
    const VxcResult left_handed = rvv10_result(mirror, {valence}, core);
    expect_same_result(left_handed, right_handed, ROUNDING_TOLERANCE);
}

TEST_F(NonlocalXCTest, WorkspaceFollowsGridCellAndProcessChanges)
{
    // The libvdwxc object is reused between calls; switching the process layout, the cell, the grid and
    // back must reproduce the earlier results.
    ModulePW::PW_Basis orthogonal{"cpu", "double"};
    orthogonal.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 rectangle(std::sqrt(0.75), 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0);
    setup(orthogonal, rectangle, BASIS_CUTOFF);
    ModulePW::PW_Basis finer{"cpu", "double"};
    finer.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
    setup_grid(finer, lattice, BASIS_CUTOFF, 5, 7, 11);
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
    expect_same_result(repeated, first, 0.0);
    expect_same_result(again, first, ROUNDING_TOLERANCE);
    expect_same_result(other_grid_again, other_grid, ROUNDING_TOLERANCE);
    EXPECT_NEAR(std::get<0>(distributed), std::get<0>(first), ROUNDING_TOLERANCE);
    EXPECT_GT(std::abs(std::get<0>(other_cell) - std::get<0>(first)), 1.0e-8);
}

TEST_F(NonlocalXCTest, DistributedStressMatchesCompleteGrid)
{
    const std::vector<double> valence = local_grid(pw, valence_density);
    const std::vector<double> core = local_grid(pw, core_density);
    const std::vector<double> valence_complete = local_grid(serial, valence_density);
    const std::vector<double> core_complete = local_grid(serial, core_density);
    const ModuleBase::matrix distributed = rvv10_stress(pw, {valence}, core);
    const ModuleBase::matrix complete = rvv10_stress(serial, {valence_complete}, core_complete);
    EXPECT_GT(std::abs(complete(0, 0)), 1.0e-10);
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            EXPECT_NEAR(distributed(a, b), complete(a, b), ROUNDING_TOLERANCE);
        }
    }
}

TEST_F(NonlocalXCTest, StressIsTheStrainDerivativeOfTheEnergy)
{
    // A strain keeps the grid and every G component of the density, so the stress is the exact derivative
    // of the discrete energy: sigma_ab = -(1/omega) dE/d(eps_ab) at a fixed electron number.
    ModulePW::PW_Basis box{"cpu", "double"};
    box.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
    setup(box, lattice, LARGEST_SPHERE_CUTOFF);
    const std::vector<double> density = local_grid(box, valence_density);
    const std::vector<double> no_core = local_grid(box, zero_density);
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
            const double h = (side == 0) ? STRAIN_STEP : -STRAIN_STEP;
            ModulePW::PW_Basis strained{"cpu", "double"};
            strained.initmpi(1, 0, MPI_COMM_SELF);
            const ModuleBase::Matrix3 strained_cell = strained_lattice(lattice, a, b, h);
            setup(strained, strained_cell, LARGEST_SPHERE_CUTOFF);
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
        const double finite_difference = -(energy[0] - energy[1]) / (2.0 * STRAIN_STEP * box.omega);
        const double analytic = extra(a, b) + ((a == b) ? diagonal : 0.0);
        const double tolerance = STRAIN_TOLERANCE * std::abs(analytic) + 1.0e-12;
        SCOPED_TRACE(testing::Message() << "component " << a << b);
        EXPECT_NEAR(analytic, finite_difference, tolerance);
    }
}

TEST_F(NonlocalXCTest, VxcAddsRvv10ToEveryFunctionalPath)
{
    // v_xc and v_xc_meta with vdw_method rvv10 equal the semilocal result plus add_rvv10, for native PBE
    // (final return), a libxc GGA (libxc branch), noncollinear PBE (NCGGA exit) and SCAN (v_xc_meta).
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
    const bool domag_z = false;
    const double hybrid_alpha = 0.0;
    const double hse_omega = 0.0;

    const int ncases = 3;
    const char* functionals[ncases] = {"PBE", "GGA_X_RPW86+GGA_C_PBE", "PBE"};
    const std::vector<std::vector<double>>* densities[ncases] = {&collinear, &collinear, &noncollinear};
    const bool domags[ncases] = {false, false, true};
    const int gga_grads[ncases] = {0, 0, 2};
    for (int icase = 0; icase < ncases; ++icase)
    {
        SCOPED_TRACE(testing::Message() << functionals[icase] << " with nspin " << densities[icase]->size());
        XC_Functional::ScopedType scope(functionals[icase]);
        const ChargeOnGrid charge(pw, *densities[icase], core);
        const int nspin = charge.chr.nspin;
        const bool domag = domags[icase];
        const int gga_grad = gga_grads[icase];
        const VxcResult total = XC_Functional::v_xc(pw.nrxx,
                                                    &charge.chr,
                                                    &ucell,
                                                    nspin,
                                                    domag,
                                                    domag_z,
                                                    gga_grad,
                                                    hybrid_alpha,
                                                    hse_omega);
        VxcResult semilocal;
        {
            Rvv10Parameters off(0.0, 0.0);
            semilocal = XC_Functional::v_xc(pw.nrxx,
                                            &charge.chr,
                                            &ucell,
                                            nspin,
                                            domag,
                                            domag_z,
                                            gga_grad,
                                            hybrid_alpha,
                                            hse_omega);
        }
        expect_rvv10_added(total, semilocal, charge.chr);
    }

    XC_Functional::ScopedType scan("SCAN");
    const std::vector<std::vector<double>> unpolarized(1, valence);
    std::vector<double> kinetic = local_grid(pw, kinetic_density);
    std::vector<double*> kinetic_rows(1, kinetic.data());
    ChargeOnGrid charge(pw, unpolarized, core);
    charge.chr.kin_r = kinetic_rows.data();
    const std::vector<int> func_id = XC_Functional::get_func_id();
    const int nspin = 1;
    typedef std::tuple<double, double, ModuleBase::matrix, ModuleBase::matrix> MetaResult;
    const MetaResult total = XC_Functional_Libxc::v_xc_meta(func_id,
                                                            pw.nrxx,
                                                            pw.omega,
                                                            pw.tpiba,
                                                            &charge.chr,
                                                            nspin,
                                                            hybrid_alpha,
                                                            hse_omega);
    MetaResult semilocal;
    {
        Rvv10Parameters off(0.0, 0.0);
        semilocal = XC_Functional_Libxc::v_xc_meta(func_id,
                                                   pw.nrxx,
                                                   pw.omega,
                                                   pw.tpiba,
                                                   &charge.chr,
                                                   nspin,
                                                   hybrid_alpha,
                                                   hse_omega);
    }
    const VxcResult total_xc(std::get<0>(total), std::get<1>(total), std::get<2>(total));
    const VxcResult semilocal_xc(std::get<0>(semilocal), std::get<1>(semilocal), std::get<2>(semilocal));
    expect_rvv10_added(total_xc, semilocal_xc, charge.chr);
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
