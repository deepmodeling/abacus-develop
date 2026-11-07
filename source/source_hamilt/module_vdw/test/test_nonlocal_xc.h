#ifndef TEST_NONLOCAL_XC_H
#define TEST_NONLOCAL_XC_H

// Fixture and helpers shared by test_nonlocal_xc.cpp and test_nonlocal_xc_paths.cpp; bodies in
// test_nonlocal_xc_helpers.cpp.

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

double valence_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz);

double core_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz);

double kinetic_density(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz);

double perturbation(const ModulePW::PW_Basis& pw, const int ix, const int iy, const int iz);

/// Values of f on the z planes of this process, in the ABACUS layout [ixy * nplane + iz].
std::vector<double> local_grid(const ModulePW::PW_Basis& pw, GridFunction f);

/// Oblique lattice; z_sign = -1 gives its left-handed mirror image.
ModuleBase::Matrix3 oblique_lattice(const double z_sign);

/// The lattice after the strain eps_ab = h: every lattice vector gets a_i[a] += h * a_i[b].
ModuleBase::Matrix3 strained_lattice(const ModuleBase::Matrix3& lattice,
                                     const int a,
                                     const int b,
                                     const double h);

/// Basis on the nx x ny x nz FFT grid of the given cell.
void setup_grid(ModulePW::PW_Basis& basis,
                const ModuleBase::Matrix3& lattice,
                const double cutoff,
                const int nx,
                const int ny,
                const int nz);

/// 5x7x9 grid: on two processes the z planes split 5 + 4, on four processes 3 + 2 + 2 + 2 against
/// libvdwxc's 3 + 3 + 3 + 0, which exercises the redistribution and an empty slab. Odd ny and nz leave no
/// Nyquist plane on libvdwxc's halved axis (ABACUS y in the adapter's (z, x, y) order, ABACUS z in the
/// direct call's (x, y, z) order), so libvdwxc's result does not depend on the axis order.
void setup(ModulePW::PW_Basis& basis, const ModuleBase::Matrix3& lattice, const double cutoff);

/// XC runtime parameters with the rVV10 switch set as given and the parameters of the rPW86-PBE pairing.
XCFunctionalParameters rvv10_parameters(const bool rvv10);

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
VxcResult rvv10_result(ModulePW::PW_Basis& pw,
                       const std::vector<std::vector<double>>& channels,
                       const std::vector<double>& core);

/// rVV10 stress beyond the XC diagonal term for the given density on pw.
ModuleBase::matrix rvv10_stress(ModulePW::PW_Basis& pw,
                                const std::vector<std::vector<double>>& channels,
                                const std::vector<double>& core);

/// etxc and vtxc get scalar_tolerance, the potential point_tolerance: the scalars are sums whose order
/// may depend on the thread count, the potential is point-wise.
void expect_same_result(const VxcResult& actual,
                        const VxcResult& expected,
                        const double scalar_tolerance,
                        const double point_tolerance);

/// A functional evaluated with rVV10 equals the evaluation without it plus add_rvv10; point_tolerance 0.0
/// where the functional hands its gradient to add_rvv10.
void expect_rvv10_added(const VxcResult& total,
                        VxcResult semilocal,
                        const Charge& chr,
                        const double point_tolerance);

/// XC_Functional::v_xc for the given functional and density with and without rVV10.
void expect_vxc_adds_rvv10(const char* functional,
                           const Charge& chr,
                           const UnitCell& ucell,
                           const bool domag,
                           const int gga_grad,
                           const double point_tolerance);

/// The result on the distributed grid equals the result on the complete grid, on the planes of this
/// process, for the energy, vtxc, the potential and the stress.
void expect_distributed_matches_complete(ModulePW::PW_Basis& distributed, ModulePW::PW_Basis& complete);

class NonlocalXCTest : public testing::Test
{
  protected:
    /// grid on an oblique cell, distributed over the processes of the test
    ModulePW::PW_Basis pw{"cpu", "double"};
    /// the same grid, complete on every process
    ModulePW::PW_Basis serial{"cpu", "double"};

    void SetUp() override;
};
} // namespace rvv10_test

#endif
