// nonlocal_xc.cpp, the rVV10 adapter to libvdwxc: libvdwxc called directly as the reference, the potential and the
// stress as derivatives of the energy, every XC evaluator that adds the term, and the exchange between the ABACUS
// planes and the libvdwxc slabs on the processes of the test.
#include "source_hamilt/module_vdw/nonlocal_xc.h"

#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_device.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_basis/module_pw/pw_basis_big.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_xc/xc_functional.h"
#ifdef __LIBXC
#include "source_hamilt/module_xc/libxc_abacus.h"
#endif

#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <tuple>
#include <vector>

// The reference calls libvdwxc directly. Include its C header after the C++ headers.
#include <complex.h>
#include <vdwxc.h>
#undef I

// The target links PW/XC objects only; the tests own the charge storage.
Charge::Charge() {}
Charge::~Charge() {}
UnitCell::UnitCell() {}
UnitCell::~UnitCell() {}
Magnetism::Magnetism() {}
Magnetism::~Magnetism() {}
SepPot::SepPot() {}
SepPot::~SepPot() {}
Sep_Cell::Sep_Cell() noexcept {}
Sep_Cell::~Sep_Cell() noexcept {}

namespace
{
typedef std::tuple<double, double, ModuleBase::matrix> VxcResult;
// rVV10 parameters of the rPW86-PBE pairing
const double rvv10_b_rpw86 = 6.3;
const double rvv10_c_rpw86 = 0.0093;
// Independent evaluations of one density differ by rounding only (FFT layout, order of the sums).
const double rounding_tolerance = 2.0e-12;
const double basis_cutoff = 4.0;

/// A smooth test field; the arguments are the grid coordinates as phases 2 pi i / n.
typedef double (*GridFunction)(const double ax, const double ay, const double az);

double valence_density(const double ax, const double ay, const double az)
{
    const double ayz = ay + az;
    return 0.004 + 0.006 * (1.0 + 0.25 * std::sin(ax) + 0.15 * std::cos(ayz));
}

double core_density(const double, const double, const double az) { return 0.002 * (1.0 + std::cos(az)); }

#ifdef __LIBXC
double kinetic_density(const double, const double ay, const double) { return 0.003 * (1.0 + 0.2 * std::cos(ay)); }
#endif

double perturbation(const double ax, const double ay, const double)
{
    const double axy = ax + ay;
    const double ay2 = 2.0 * ay;
    return std::sin(axy) + 0.3 * std::cos(ay2);
}

/// Values of f on the z planes of this process, in the ABACUS layout [ixy * nplane + iz].
std::vector<double> local_grid(const ModulePW::PW_Basis &pw, GridFunction f)
{
    std::vector<double> values(pw.nrxx);
    for (int ix = 0; ix < pw.nx; ++ix)
    {
        const double ax = ModuleBase::TWO_PI * ix / pw.nx;
        for (int iy = 0; iy < pw.ny; ++iy)
        {
            const double ay = ModuleBase::TWO_PI * iy / pw.ny;
            for (int iz = 0; iz < pw.nplane; ++iz)
            {
                const int global_z = pw.startz_current + iz;
                const double az = ModuleBase::TWO_PI * global_z / pw.nz;
                values[(ix * pw.ny + iy) * pw.nplane + iz] = f(ax, ay, az);
            }
        }
    }
    return values;
}

/// The valence density as two spin channels (nspin = 2) or as the density and a magnetization (nspin = 4).
std::vector<std::vector<double>> spin_components(const std::vector<double> &valence, const int nspin)
{
    std::vector<std::vector<double>> components(nspin, valence);
    for (std::size_t ir = 0; ir < valence.size(); ++ir)
    {
        if (nspin == 2)
        {
            components[0][ir] = 0.7 * valence[ir];
            components[1][ir] = valence[ir] - components[0][ir];
        }
        else
        {
            components[1][ir] = 0.1 * valence[ir];
            components[2][ir] = -0.05 * valence[ir];
            components[3][ir] = 0.2 * valence[ir];
        }
    }
    return components;
}

/// Oblique lattice; z_sign = -1 gives its left-handed mirror image.
ModuleBase::Matrix3 oblique_lattice(const double z_sign)
{
    const double c3 = 2.0 * z_sign;
    const double b2 = std::sqrt(0.75);
    return ModuleBase::Matrix3(1.0, 0.0, 0.0, -0.5, b2, 0.0, 0.0, 0.0, c3);
}

/// The lattice after the strain eps_ab = h: every lattice vector gets a_i[a] += h * a_i[b].
ModuleBase::Matrix3 strained_lattice(
    const ModuleBase::Matrix3 &lattice,
    const int a,
    const int b,
    const double h)
{
    double m[3][3] = {{lattice.e11, lattice.e12, lattice.e13}, {lattice.e21, lattice.e22, lattice.e23},
        {lattice.e31, lattice.e32, lattice.e33}};
    for (int i = 0; i < 3; ++i)
    {
        m[i][a] += h * m[i][b];
    }
    return ModuleBase::Matrix3(m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]);
}

/// Basis on the 5 x 7 x nz FFT grid of the given cell. With nz = 9 the z planes split 5 + 4 on two processes and
/// 3 + 2 + 2 + 2 on four, against libvdwxc's 3 + 3 + 3 + 0, which exercises the redistribution and an empty slab.
/// Odd ny and nz leave no Nyquist plane on libvdwxc's halved axis (ABACUS y in the adapter's (z, x, y) order,
/// ABACUS z in the direct call's (x, y, z) order), so libvdwxc's result does not depend on the axis order.
void setup_grid(
    ModulePW::PW_Basis &basis,
    const ModuleBase::Matrix3 &lattice,
    const double cutoff,
    const int nz)
{
    basis.initgrids(10.0, lattice, 5, 7, nz);
    basis.initparameters(false, cutoff, 1, true);
    basis.setuptransform();
    basis.collect_local_pw();
}

/// XC runtime parameters with the rVV10 switch set as given and the parameters of the rPW86-PBE pairing.
XCFunctionalParameters rvv10_parameters(const bool rvv10)
{
    XCFunctionalParameters parameters = XC_Functional::get_runtime_parameters();
    parameters.rvv10 = rvv10;
    parameters.rvv10_b = rvv10_b_rpw86;
    parameters.rvv10_c = rvv10_c_rpw86;
    return parameters;
}

/// Sets the XC runtime parameters for the lifetime of the object, as ESolver_FP does.
class RuntimeParameters
{
  public:
    explicit RuntimeParameters(const bool rvv10)
    {
        saved_ = XC_Functional::get_runtime_parameters();
        const XCFunctionalParameters parameters = rvv10_parameters(rvv10);
        XC_Functional::set_runtime_parameters(parameters);
    }
    ~RuntimeParameters() { XC_Functional::set_runtime_parameters(saved_); }

  private:
    XCFunctionalParameters saved_;
};

/// A Charge on pw that owns its density components, core density and the G-space buffers of gradcorr.
class ChargeOnGrid : public Charge
{
  public:
    ChargeOnGrid(
        ModulePW::PW_Basis &pw, const std::vector<std::vector<double>> &channels, const std::vector<double> &core)
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
        rhopw = &pw;
        nrxx = pw.nrxx;
        nspin = static_cast<int>(ncomponents);
        rho = rows_.data();
        rhog = rhog_rows_.data();
        rho_core = core_.data();
        rhog_core = core_g_.data();
    }

  private:
    std::vector<std::vector<double>> channels_;
    std::vector<double *> rows_;
    std::vector<double> core_;
    std::vector<std::vector<std::complex<double>>> rhog_;
    std::vector<std::complex<double> *> rhog_rows_;
    std::vector<std::complex<double>> core_g_;
};

/// rVV10 contribution for the given valence channels and core density on pw.
VxcResult rvv10_result(
    ModulePW::PW_Basis &pw, const std::vector<std::vector<double>> &channels, const std::vector<double> &core)
{
    const ChargeOnGrid charge(pw, channels, core);
    const XCFunctionalParameters parameters = rvv10_parameters(true);
    double etxc = 0.0;
    double vtxc = 0.0;
    ModuleBase::matrix v(charge.nspin, pw.nrxx);
    vdw::add_rvv10(parameters, charge, nullptr, etxc, vtxc, v);
    return VxcResult(etxc, vtxc, v);
}

/// rVV10 stress beyond the XC diagonal term for the given density on pw.
ModuleBase::matrix rvv10_stress(
    ModulePW::PW_Basis &pw, const std::vector<std::vector<double>> &channels, const std::vector<double> &core)
{
    const ChargeOnGrid charge(pw, channels, core);
    const XCFunctionalParameters parameters = rvv10_parameters(true);
    ModuleBase::matrix sigma(3, 3);
    vdw::add_rvv10_stress(parameters, charge, sigma);
    return sigma;
}

/// The energies, vtxc and the potentials on every point agree to rounding.
void expect_same_result(const VxcResult &actual, const VxcResult &expected)
{
    EXPECT_NEAR(std::get<0>(actual), std::get<0>(expected), rounding_tolerance);
    EXPECT_NEAR(std::get<1>(actual), std::get<1>(expected), rounding_tolerance);
    const ModuleBase::matrix &v = std::get<2>(actual);
    const ModuleBase::matrix &w = std::get<2>(expected);
    ASSERT_EQ(v.nr, w.nr);
    ASSERT_EQ(v.nc, w.nc);
    for (int is = 0; is < v.nr; ++is)
    {
        for (int ir = 0; ir < v.nc; ++ir)
        {
            EXPECT_NEAR(v(is, ir), w(is, ir), rounding_tolerance);
        }
    }
}

/// A functional evaluated with rVV10 equals the evaluation without it plus add_rvv10.
void expect_rvv10_added(const VxcResult &total, VxcResult semilocal, const Charge &chr)
{
    const XCFunctionalParameters parameters = rvv10_parameters(true);
    double &etxc = std::get<0>(semilocal);
    double &vtxc = std::get<1>(semilocal);
    ModuleBase::matrix &v = std::get<2>(semilocal);
    vdw::add_rvv10(parameters, chr, nullptr, etxc, vtxc, v);
    expect_same_result(total, semilocal);
}

/// XC_Functional::v_xc for the given functional and density with and without rVV10; a four-component density is
/// evaluated as a magnetized one.
void expect_vxc_adds_rvv10(const char *functional, const Charge &chr, const int gga_grad)
{
    XC_Functional::ScopedType scope(functional);
    UnitCell ucell;
    ucell.omega = chr.rhopw->omega;
    ucell.tpiba = chr.rhopw->tpiba;
    ucell.magnet.lsign_ = false;
    const int nrxx = chr.nrxx;
    const int nspin = chr.nspin;
    const bool domag = nspin == 4;
    const bool domag_z = false;
    const double hybrid_alpha = 0.0;
    const double hse_omega = 0.0;
    VxcResult total;
    VxcResult semilocal;
    for (int with_rvv10 = 0; with_rvv10 < 2; ++with_rvv10)
    {
        const bool rvv10 = with_rvv10 == 1;
        RuntimeParameters parameters(rvv10);
        VxcResult &result = rvv10 ? total : semilocal;
        result = XC_Functional::v_xc(nrxx, &chr, &ucell, nspin, domag, domag_z, gga_grad, hybrid_alpha, hse_omega);
    }
    expect_rvv10_added(total, semilocal, chr);
}

/// The result on the distributed grid equals the result on the complete grid, on the planes of this process, for
/// the energy, vtxc, the potential and the stress.
void expect_distributed_matches_complete(ModulePW::PW_Basis &distributed, ModulePW::PW_Basis &complete)
{
    const std::vector<double> valence = local_grid(distributed, valence_density);
    const std::vector<double> core = local_grid(distributed, core_density);
    const std::vector<double> valence_complete = local_grid(complete, valence_density);
    const std::vector<double> core_complete = local_grid(complete, core_density);
    const VxcResult result = rvv10_result(distributed, {valence}, core);
    const VxcResult reference = rvv10_result(complete, {valence_complete}, core_complete);
    EXPECT_NEAR(std::get<0>(result), std::get<0>(reference), rounding_tolerance);
    EXPECT_NEAR(std::get<1>(result), std::get<1>(reference), rounding_tolerance);
    for (int ixy = 0; ixy < distributed.nxy; ++ixy)
    {
        for (int iz = 0; iz < distributed.nplane; ++iz)
        {
            const int local = ixy * distributed.nplane + iz;
            const int global = ixy * distributed.nz + distributed.startz_current + iz;
            EXPECT_NEAR(std::get<2>(result)(0, local), std::get<2>(reference)(0, global), rounding_tolerance);
        }
    }
    const ModuleBase::matrix stress = rvv10_stress(distributed, {valence}, core);
    const ModuleBase::matrix stress_complete = rvv10_stress(complete, {valence_complete}, core_complete);
    // a lower bound that shows the compared stress is not trivially zero
    const double nonzero_stress = 1.0e-10;
    const double stress_magnitude = std::abs(stress_complete(0, 0));
    EXPECT_GT(stress_magnitude, nonzero_stress);
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            EXPECT_NEAR(stress(a, b), stress_complete(a, b), rounding_tolerance);
        }
    }
}

class NonlocalXCTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        // the XC evaluators reduce over POOL_WORLD, which the shared main of the MPI tests leaves unset
        const ModuleBase::CommunicationDomain world = ModuleBase::world_comm_domain();
        const int nproc = world.size();
        const int rank = world.rank();
        POOL_WORLD = world.communicator();
        // The oblique cell catches row/column mistakes in the cell handed to libvdwxc.
        const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
        pw.initmpi(nproc, rank, POOL_WORLD);
        serial.initmpi(1, 0, MPI_COMM_SELF);
        setup_grid(pw, lattice, basis_cutoff, 9);
        setup_grid(serial, lattice, basis_cutoff, 9);
    }

    /// grid on an oblique cell, distributed over the processes of the test
    ModulePW::PW_Basis pw{"cpu", "double"};
    /// the same grid, complete on every process
    ModulePW::PW_Basis serial{"cpu", "double"};
};

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
    expect_same_result(result, expected);

    // b and C reach libvdwxc: other values give the energy that libvdwxc gives for them
    const double other_b = 15.7;
    const double other_c = 0.02;
    // a lower bound that shows the energy depends on the parameters
    const double distinct_energy = 1.0e-8;
    XCFunctionalParameters other = rvv10_parameters(true);
    other.rvv10_b = other_b;
    other.rvv10_c = other_c;
    const ChargeOnGrid charge(serial, {valence}, core);
    double other_etxc = 0.0;
    double other_vtxc = 0.0;
    ModuleBase::matrix other_v(1, serial.nrxx);
    vdw::add_rvv10(other, charge, nullptr, other_etxc, other_vtxc, other_v);
    direct = vdwxc_new(VDWXC_RVV10);
    ASSERT_EQ(vdwxc_set_rvv10_parameters(direct, other_b, other_c), 0);
    vdwxc_set_unit_cell(direct, nx, ny, nz, a.e11, a.e12, a.e13, a.e21, a.e22, a.e23, a.e31, a.e32, a.e33);
    vdwxc_init_serial(direct);
    const double other_energy_ha = vdwxc_calculate(direct, total.data(), sigma.data(), dedn.data(), dedsigma.data());
    vdwxc_finalize(&direct);
    const double other_expected = ModuleBase::e2 * other_energy_ha;
    EXPECT_NEAR(other_etxc, other_expected, rounding_tolerance);
    const double energy_change = std::abs(other_etxc - etxc);
    EXPECT_GT(energy_change, distinct_energy);

    // spin channels share the total density: same energy and vtxc, and each channel gets the potential
    const std::vector<std::vector<double>> collinear = spin_components(valence, 2);
    const VxcResult polarized = rvv10_result(serial, collinear, core);
    EXPECT_NEAR(std::get<0>(polarized), etxc, rounding_tolerance);
    EXPECT_NEAR(std::get<1>(polarized), vtxc, rounding_tolerance);
    for (int ir = 0; ir < serial.nrxx; ++ir)
    {
        EXPECT_NEAR(std::get<2>(polarized)(0, ir), v(0, ir), rounding_tolerance);
        EXPECT_NEAR(std::get<2>(polarized)(1, ir), v(0, ir), rounding_tolerance);
    }

    // four-component layout: rho[0] is the total valence density, and only v(0) gets the potential
    const std::vector<std::vector<double>> magnetized = spin_components(valence, 4);
    const VxcResult noncollinear = rvv10_result(serial, magnetized, core);
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
    // central difference of the energy along a smooth perturbation: the truncation error is of order step^2
    const double finite_difference_step = 1.0e-5;
    const double finite_difference_tolerance = 2.0e-6;
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
    const double strain_step = 1.0e-4;
    const double strain_tolerance = 1.0e-5;
    // PW_Basis lowers a larger cutoff to the largest G sphere inside the FFT box. The test densities have
    // only first-shell G components, well inside that sphere, so a small strain changes none of them.
    const double largest_sphere_cutoff = 1.0e6;
    ModulePW::PW_Basis box{"cpu", "double"};
    box.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
    setup_grid(box, lattice, largest_sphere_cutoff, 9);
    const std::vector<double> density = local_grid(box, valence_density);
    const std::vector<double> no_core(box.nrxx, 0.0);
    const VxcResult reference = rvv10_result(box, {density}, no_core);
    const ModuleBase::matrix extra = rvv10_stress(box, {density}, no_core);
    // the XC diagonal term that ABACUS adds from etxc and vtxc
    const double diagonal = -(std::get<0>(reference) - std::get<1>(reference)) / box.omega;

    const int components[4][2] = {{0, 0}, {2, 2}, {0, 1}, {1, 2}};
    for (const int *ab : components)
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
            setup_grid(strained, strained_cell, largest_sphere_cutoff, 9);
            // keep the electron number: the density scales with the inverse volume ratio
            const double volume_ratio = strained.omega / box.omega;
            std::vector<double> scaled = density;
            for (double &value : scaled)
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
    // Libxc a libxc GGA (v_xc_libxc) in the same layouts and unpolarized, and SCAN (v_xc_meta).
    const std::vector<double> valence = local_grid(pw, valence_density);
    const std::vector<double> core = local_grid(pw, core_density);
    const std::vector<std::vector<double>> collinear = spin_components(valence, 2);
    const std::vector<std::vector<double>> noncollinear = spin_components(valence, 4);
    const ChargeOnGrid polarized(pw, collinear, core);
    const ChargeOnGrid magnetized(pw, noncollinear, core);
    expect_vxc_adds_rvv10("PBE", polarized, 0);
    expect_vxc_adds_rvv10("PBE", magnetized, 2);
    expect_vxc_adds_rvv10("PBE", magnetized, 0);
#ifdef __LIBXC
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", polarized, 0);
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", magnetized, 0);
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", magnetized, 2);
    const std::vector<std::vector<double>> unpolarized(1, valence);
    const ChargeOnGrid single(pw, unpolarized, core);
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", single, 0);

    XC_Functional::ScopedType scan("SCAN");
    std::vector<double> kinetic = local_grid(pw, kinetic_density);
    std::vector<double *> kinetic_rows(1, kinetic.data());
    ChargeOnGrid charge(pw, unpolarized, core);
    charge.kin_r = kinetic_rows.data();
    const std::vector<int> func_id = XC_Functional::get_func_id();
    const int nspin = 1;
    const double hybrid_alpha = 0.0;
    const double hse_omega = 0.0;
    typedef std::tuple<double, double, ModuleBase::matrix, ModuleBase::matrix> MetaResult;
    VxcResult results[2];
    for (int with_rvv10 = 0; with_rvv10 < 2; ++with_rvv10)
    {
        const bool rvv10 = with_rvv10 == 1;
        RuntimeParameters parameters(rvv10);
        const MetaResult meta = XC_Functional_Libxc::v_xc_meta(
            func_id, pw.nrxx, pw.omega, pw.tpiba, &charge, nspin, hybrid_alpha, hse_omega);
        results[with_rvv10] = VxcResult(std::get<0>(meta), std::get<1>(meta), std::get<2>(meta));
    }
    expect_rvv10_added(results[1], results[0], charge);
#endif
}

TEST_F(NonlocalXCTest, DistributedGridMatchesCompleteGrid)
{
    // the exchange between the ABACUS planes and the libvdwxc slabs, uneven and with an empty slab on four processes
    expect_distributed_matches_complete(pw, serial);
}

TEST_F(NonlocalXCTest, BigGridDistributionMatchesCompleteGrid)
{
    // The density grid of a run is a PW_Basis_Big, which hands out whole blocks of bz planes: on two
    // processes 8 + 4 against libvdwxc's 6 + 6, on four 4 + 4 + 4 + 0 against 3 + 3 + 3 + 3, so that a
    // process may own no plane at all and one ABACUS block may span two libvdwxc slabs.
    ModulePW::PW_Basis_Big big{"cpu", "double"};
    big.setbxyz(1, 1, 4);
    big.initmpi(pw.poolnproc, pw.poolrank, pw.pool_world);
    ModulePW::PW_Basis complete{"cpu", "double"};
    complete.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
    setup_grid(big, lattice, basis_cutoff, 12);
    setup_grid(complete, lattice, basis_cutoff, 12);
    ASSERT_EQ(big.nz, 12);
    ASSERT_EQ(big.nplane % 4, 0);
    if (pw.poolnproc == 4)
    {
        EXPECT_EQ(big.numz[3], 0);
    }
    expect_distributed_matches_complete(big, complete);
    // the XC evaluator calls add_rvv10 on every process of the pool, also on one without points
    const std::vector<double> valence = local_grid(big, valence_density);
    const std::vector<double> core = local_grid(big, core_density);
    const std::vector<std::vector<double>> unpolarized(1, valence);
    const ChargeOnGrid single(big, unpolarized, core);
#ifdef __LIBXC
    expect_vxc_adds_rvv10("GGA_X_RPW86+GGA_C_PBE", single, 0);
#else
    expect_vxc_adds_rvv10("PBE", single, 0);
#endif
}

TEST_F(NonlocalXCTest, LeftHandedCellGivesTheMirrorImageResult)
{
    // A mirror image has the same result; libvdwxc's signed cell volume would flip the sign of the energy.
    ModulePW::PW_Basis mirror{"cpu", "double"};
    mirror.initmpi(1, 0, MPI_COMM_SELF);
    const ModuleBase::Matrix3 mirrored = oblique_lattice(-1.0);
    setup_grid(mirror, mirrored, basis_cutoff, 9);
    ASSERT_LT(mirror.latvec.Det(), 0.0);
    const std::vector<double> valence = local_grid(serial, valence_density);
    const std::vector<double> core = local_grid(serial, core_density);
    const VxcResult right_handed = rvv10_result(serial, {valence}, core);
    const VxcResult left_handed = rvv10_result(mirror, {valence}, core);
    expect_same_result(left_handed, right_handed);
}
} // namespace
