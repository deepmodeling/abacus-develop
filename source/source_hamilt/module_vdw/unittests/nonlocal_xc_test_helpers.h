#ifndef NONLOCAL_XC_TEST_HELPERS_H
#define NONLOCAL_XC_TEST_HELPERS_H

// Fixture and helpers shared by test_nonlocal_xc.cpp and test_rvv10_workspace.cpp.

#include "../nonlocal_xc.h"
#include "source_base/constants.h"
#include "source_base/matrix.h"
#include "source_base/matrix3.h"
#include "source_base/parallel_device.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_xc/xc_functional.h"

#include <cmath>
#include <complex>
#include <gtest/gtest.h>
#include <tuple>
#include <vector>

namespace rvv10_test
{
/// number of processes and rank of the test, set in main
extern int test_rank;
extern int test_size;
typedef std::tuple<double, double, ModuleBase::matrix> VxcResult;
typedef std::tuple<double, double, ModuleBase::matrix, ModuleBase::matrix> MetaResult;
// rVV10 parameters of the rPW86-PBE pairing
const double rvv10_b_rpw86 = 6.3;
const double rvv10_c_rpw86 = 0.0093;
// Results that share the plane-wave gradient agree to rounding.
const double rounding_tolerance = 2.0e-12;
// Central differences on smooth densities: the truncation error is of order step^2.
const double finite_difference_step = 1.0e-5;
const double finite_difference_tolerance = 2.0e-6;
const double strain_step = 1.0e-4;
const double strain_tolerance = 1.0e-5;
const double basis_cutoff = 4.0;
// PW_Basis lowers a larger cutoff to the largest G sphere inside the FFT box. The test densities have
// only first-shell G components, well inside that sphere, so a small strain changes none of them.
const double largest_sphere_cutoff = 1.0e6;

typedef double (*GridFunction)(const ModulePW::PW_Basis&, const int, const int, const int);

inline double valence_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz);

inline double core_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz);

inline double kinetic_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz);

inline double perturbation(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz);

/// Values of f on the z planes of this process, in the ABACUS layout [ixy * nplane + iz].
inline std::vector<double> local_grid(const ModulePW::PW_Basis& pw, GridFunction f);

/// Oblique lattice; z_sign = -1 gives its left-handed mirror image.
inline ModuleBase::Matrix3 oblique_lattice(const double z_sign);

/// The lattice after the strain eps_ab = h: every lattice vector gets a_i[a] += h * a_i[b].
inline ModuleBase::Matrix3 strained_lattice(const ModuleBase::Matrix3& lattice,
                                     const int a,
                                     const int b,
                                     const double h);

/// Basis on the nx x ny x nz FFT grid of the given cell.
inline void setup_grid(ModulePW::PW_Basis& basis,
                const ModuleBase::Matrix3& lattice,
                const double cutoff,
                const int nx,
                const int ny,
                const int nz);

/// 5x7x9 grid: on two processes the z planes split 5 + 4, on four processes 3 + 2 + 2 + 2 against
/// libvdwxc's 3 + 3 + 3 + 0, which exercises the redistribution and an empty slab. Odd ny and nz leave no
/// Nyquist plane on libvdwxc's halved axis (ABACUS y in the adapter's (z, x, y) order, ABACUS z in the
/// direct call's (x, y, z) order), so libvdwxc's result does not depend on the axis order.
inline void setup(ModulePW::PW_Basis& basis, const ModuleBase::Matrix3& lattice, const double cutoff);

/// XC runtime parameters with the rVV10 switch set as given and the parameters of the rPW86-PBE pairing.
inline XCFunctionalParameters rvv10_parameters(const bool rvv10);

/// Sets the XC runtime parameters for the lifetime of the object, as ESolver_FP does.
class RuntimeParameters
{
  public:
    explicit RuntimeParameters(const bool rvv10);
    ~RuntimeParameters();

  private:
    XCFunctionalParameters saved_;
};

/// A Charge on pw that owns its density components, core density and the G-space buffers of gradcorr.
class ChargeOnGrid
{
  public:
    ChargeOnGrid(ModulePW::PW_Basis& pw,
                 const std::vector<std::vector<double>>& channels,
                 const std::vector<double>& core);
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
inline VxcResult rvv10_result(ModulePW::PW_Basis& pw,
                       const std::vector<std::vector<double>>& channels,
                       const std::vector<double>& core);

/// rVV10 stress beyond the XC diagonal term for the given density on pw.
inline ModuleBase::matrix rvv10_stress(ModulePW::PW_Basis& pw,
                                const std::vector<std::vector<double>>& channels,
                                const std::vector<double>& core);

/// etxc and vtxc get scalar_tolerance, the potential point_tolerance: the scalars are sums whose order
/// may depend on the thread count, the potential is point-wise.
inline void expect_same_result(const VxcResult& actual,
                        const VxcResult& expected,
                        const double scalar_tolerance,
                        const double point_tolerance);

/// A functional evaluated with rVV10 equals the evaluation without it plus add_rvv10; point_tolerance 0.0
/// where the functional hands its gradient to add_rvv10.
inline void expect_rvv10_added(const VxcResult& total,
                        VxcResult semilocal,
                        const Charge& chr,
                        const double point_tolerance);

/// XC_Functional::v_xc for the given functional and density with and without rVV10.
inline void expect_vxc_adds_rvv10(const char* functional,
                           const Charge& chr,
                           const UnitCell& ucell,
                           const bool domag,
                           const int gga_grad,
                           const double point_tolerance);

/// The result on the distributed grid equals the result on the complete grid, on the planes of this
/// process, for the energy, vtxc, the potential and the stress.
inline void expect_distributed_matches_complete(ModulePW::PW_Basis& distributed, ModulePW::PW_Basis& complete);

class NonlocalXCTest : public testing::Test
{
  protected:
    /// grid on an oblique cell, distributed over the processes of the test
    ModulePW::PW_Basis pw{"cpu", "double"};
    /// the same grid, complete on every process
    ModulePW::PW_Basis serial{"cpu", "double"};

    void SetUp() override;
};

// ---- definitions (header-only: shared test helpers live in headers, AGENTS.md rule 11) ----

inline double valence_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ax = ModuleBase::TWO_PI * ix / pw.nx;
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    const double az = ModuleBase::TWO_PI * iz / pw.nz;
    const double ayz = ay + az;
    return 0.004 + 0.006 * (1.0 + 0.25 * std::sin(ax) + 0.15 * std::cos(ayz));
}

inline double core_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double az = ModuleBase::TWO_PI * iz / pw.nz;
    return 0.002 * (1.0 + std::cos(az));
}

inline double kinetic_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    return 0.003 * (1.0 + 0.2 * std::cos(ay));
}

inline double perturbation(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz)
{
    const double ax = ModuleBase::TWO_PI * ix / pw.nx;
    const double ay = ModuleBase::TWO_PI * iy / pw.ny;
    const double axy = ax + ay;
    const double ay2 = 2.0 * ay;
    return std::sin(axy) + 0.3 * std::cos(ay2);
}

inline std::vector<double> local_grid(const ModulePW::PW_Basis& pw, GridFunction f)
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

inline ModuleBase::Matrix3 oblique_lattice(const double z_sign)
{
    const double c3 = 2.0 * z_sign;
    return ModuleBase::Matrix3(1.0, 0.0, 0.0, -0.5, std::sqrt(0.75), 0.0, 0.0, 0.0, c3);
}

inline ModuleBase::Matrix3 strained_lattice(const ModuleBase::Matrix3& lattice,
                                     const int a,
                                     const int b,
                                     const double h)
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

inline void setup_grid(ModulePW::PW_Basis& basis,
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

inline void setup(ModulePW::PW_Basis& basis, const ModuleBase::Matrix3& lattice, const double cutoff)
{
    setup_grid(basis, lattice, cutoff, 5, 7, 9);
}

inline XCFunctionalParameters rvv10_parameters(const bool rvv10)
{
    XCFunctionalParameters parameters = XC_Functional::get_runtime_parameters();
    parameters.rvv10 = rvv10;
    parameters.rvv10_b = rvv10_b_rpw86;
    parameters.rvv10_c = rvv10_c_rpw86;
    return parameters;
}

inline VxcResult rvv10_result(ModulePW::PW_Basis& pw,
                       const std::vector<std::vector<double>>& channels,
                       const std::vector<double>& core)
{
    const ChargeOnGrid charge(pw, channels, core);
    const XCFunctionalParameters parameters = rvv10_parameters(true);
    double etxc = 0.0;
    double vtxc = 0.0;
    ModuleBase::matrix v(charge.chr.nspin, pw.nrxx);
    vdw::add_rvv10(parameters, charge.chr, etxc, vtxc, v);
    return VxcResult(etxc, vtxc, v);
}

inline ModuleBase::matrix rvv10_stress(ModulePW::PW_Basis& pw,
                                const std::vector<std::vector<double>>& channels,
                                const std::vector<double>& core)
{
    const ChargeOnGrid charge(pw, channels, core);
    const XCFunctionalParameters parameters = rvv10_parameters(true);
    ModuleBase::matrix sigma(3, 3);
    vdw::add_rvv10_stress(parameters, charge.chr, sigma);
    return sigma;
}

inline void expect_same_result(const VxcResult& actual,
                        const VxcResult& expected,
                        const double scalar_tolerance,
                        const double point_tolerance)
{
    EXPECT_NEAR(std::get<0>(actual), std::get<0>(expected), scalar_tolerance);
    EXPECT_NEAR(std::get<1>(actual), std::get<1>(expected), scalar_tolerance);
    const ModuleBase::matrix& v = std::get<2>(actual);
    const ModuleBase::matrix& w = std::get<2>(expected);
    ASSERT_EQ(v.nr, w.nr);
    ASSERT_EQ(v.nc, w.nc);
    for (int is = 0; is < v.nr; ++is)
    {
        for (int ir = 0; ir < v.nc; ++ir)
        {
            EXPECT_NEAR(v(is, ir), w(is, ir), point_tolerance);
        }
    }
}

inline void expect_rvv10_added(const VxcResult& total,
                        VxcResult semilocal,
                        const Charge& chr,
                        const double point_tolerance)
{
    const XCFunctionalParameters parameters = rvv10_parameters(true);
    double& etxc = std::get<0>(semilocal);
    double& vtxc = std::get<1>(semilocal);
    ModuleBase::matrix& v = std::get<2>(semilocal);
    vdw::add_rvv10(parameters, chr, etxc, vtxc, v);
    expect_same_result(total, semilocal, rounding_tolerance, point_tolerance);
}

inline void expect_vxc_adds_rvv10(const char* functional,
                           const Charge& chr,
                           const UnitCell& ucell,
                           const bool domag,
                           const int gga_grad,
                           const double point_tolerance)
{
    XC_Functional::ScopedType scope(functional);
    const int nrxx = chr.nrxx;
    const int nspin = chr.nspin;
    const bool domag_z = false;
    const double hybrid_alpha = 0.0;
    const double hse_omega = 0.0;
    VxcResult total;
    VxcResult semilocal;
    for (int with_rvv10 = 0; with_rvv10 < 2; ++with_rvv10)
    {
        const bool rvv10 = with_rvv10 == 1;
        RuntimeParameters parameters(rvv10);
        VxcResult& result = rvv10 ? total : semilocal;
        result = XC_Functional::v_xc(nrxx, &chr, &ucell, nspin, domag, domag_z, gga_grad, hybrid_alpha, hse_omega);
    }
    expect_rvv10_added(total, semilocal, chr, point_tolerance);
}

inline void expect_distributed_matches_complete(ModulePW::PW_Basis& distributed, ModulePW::PW_Basis& complete)
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
    EXPECT_GT(std::abs(stress_complete(0, 0)), 1.0e-10);
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < 3; ++b)
        {
            EXPECT_NEAR(stress(a, b), stress_complete(a, b), rounding_tolerance);
        }
    }
}

inline ChargeOnGrid::ChargeOnGrid(ModulePW::PW_Basis& pw,
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

inline RuntimeParameters::RuntimeParameters(const bool rvv10)
{
    saved_ = XC_Functional::get_runtime_parameters();
    const XCFunctionalParameters parameters = rvv10_parameters(rvv10);
    XC_Functional::set_runtime_parameters(parameters);
}

inline RuntimeParameters::~RuntimeParameters()
{
    XC_Functional::set_runtime_parameters(saved_);
}

inline void NonlocalXCTest::SetUp()
{
    // The oblique cell catches row/column mistakes in the cell handed to libvdwxc.
    const ModuleBase::Matrix3 lattice = oblique_lattice(1.0);
    pw.initmpi(test_size, test_rank, MPI_COMM_WORLD);
    serial.initmpi(1, 0, MPI_COMM_SELF);
    setup(pw, lattice, basis_cutoff);
    setup(serial, lattice, basis_cutoff);
}
} // namespace rvv10_test

#endif // NONLOCAL_XC_TEST_HELPERS_H
