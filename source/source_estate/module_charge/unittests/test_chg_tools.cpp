#include "gtest/gtest.h"

#include "source_cell/unitcell.h"
#include "source_estate/module_charge/chg_tools.h"

#include <algorithm>
#include <cmath>
#include <vector>

// chg_tools.cpp references UnitCell (set_rho_core), so the test binary links
// the cell_info objects whose unitcell.cpp needs Magnetism symbols. Provide
// the same lightweight mocks as test_charge.cpp.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of module_charge/chg_tools.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - cal_rho2ne: integrate a single spin channel over the grid and scale
 *     by omega / nxyz to obtain the electron number
 *   - check_rho: nspin == 1/4 total-density check, nspin == 2 spin-up/down
 *     checks, mismatch warning path and negative-channel abort path
 *   - core_charge_form_factor: the shared per-point radial kernel against
 *     the closed-form Gaussian transform -- this pins the ABSOLUTE radial
 *     normalization (4pi/Omega and the Simpson mesh convention), which the
 *     DFPT localization/phase regressions do not constrain independently
 */

class ChgToolsTest : public ::testing::Test
{
  protected:
    const int nrxx = 8;       ///< local real-space grid points
    const int nxyz = 8;       ///< global real-space grid points
    const double omega = 2.0; ///< cell volume
    const double nelec = 4.0; ///< target electron number

    std::vector<double> rho_up;
    std::vector<double> rho_dn;
    std::vector<double*> rho;

    void SetUp() override
    {
        // uniform value 2.0 integrates to 8 * 2.0 * 2.0 / 8 = 4.0 electrons
        rho_up.assign(nrxx, 2.0);
        rho_dn.assign(nrxx, 2.0);
        rho.resize(2);
        rho[0] = rho_up.data();
        rho[1] = rho_dn.data();
    }
};

TEST_F(ChgToolsTest, CalRho2ne)
{
    EXPECT_NEAR(module_charge::cal_rho2ne(rho[0], nrxx, omega, nxyz), nelec, 1e-12);
}

TEST_F(ChgToolsTest, CheckRhoNonSpinMatched)
{
    module_charge::check_rho(rho.data(), 1, nrxx, omega, nxyz, nelec);
}

TEST_F(ChgToolsTest, CheckRhoSocTreatedAsTotal)
{
    module_charge::check_rho(rho.data(), 4, nrxx, omega, nxyz, nelec);
}

TEST_F(ChgToolsTest, CheckRhoNonSpinMismatchWarns)
{
    // total 4.0 differs from the target 4.5: a warning is emitted but the
    // call returns normally
    module_charge::check_rho(rho.data(), 1, nrxx, omega, nxyz, 4.5);
}

TEST_F(ChgToolsTest, CheckRhoSpin2Matched)
{
    // 2.0 spin-up + 2.0 spin-down electrons
    module_charge::check_rho(rho.data(), 2, nrxx, omega, nxyz, nelec);
}

TEST_F(ChgToolsTest, CheckRhoSpin2MismatchWarns)
{
    // spin-down integrates to 1.0 electron, total 3.0 vs target 4.0
    std::fill(rho_dn.begin(), rho_dn.end(), 0.5);
    module_charge::check_rho(rho.data(), 2, nrxx, omega, nxyz, nelec);
}

TEST_F(ChgToolsTest, CheckRhoNegativeSpinUpAborts)
{
    std::fill(rho_up.begin(), rho_up.end(), -0.5);
    EXPECT_DEATH(module_charge::check_rho(rho.data(), 2, nrxx, omega, nxyz, nelec), "");
}

TEST_F(ChgToolsTest, CheckRhoNegativeSpinDownAborts)
{
    std::fill(rho_dn.begin(), rho_dn.end(), -0.5);
    EXPECT_DEATH(module_charge::check_rho(rho.data(), 2, nrxx, omega, nxyz, nelec), "");
}

TEST_F(ChgToolsTest, CheckRhoUnsupportedNspinNoop)
{
    // nspin values other than 1/2/4 are silently skipped, as in the
    // original Charge::check_rho
    module_charge::check_rho(rho.data(), 3, nrxx, omega, nxyz, nelec);
}

TEST_F(ChgToolsTest, CoreChargeFormFactorMatchesAnalyticGaussian)
{
    // rho_atc(r) = exp(-alpha r^2) on a uniform radial grid; the closed
    // form of its spherical transform is
    //   Int r^2 exp(-alpha r^2) j0(k r) dr = sqrt(pi)/(4 alpha^{3/2}) * exp(-k^2/(4 alpha))
    // so the form factor is 4pi/omega times that. Checks K = 0 and several
    // finite K values with ONE scratch buffer reused across the calls.
    const double pi = 3.14159265358979323846;
    const double alpha = 0.5;
    const double dr = 0.005;
    const int mesh = 2001; // r_max = 10: exp(-50) ~ 2e-22, truncated safely
    std::vector<double> r(mesh);
    std::vector<double> rab(mesh, dr);
    std::vector<double> rho_atc(mesh);
    for (int i = 0; i < mesh; ++i)
    {
        r[i] = i * dr;
        rho_atc[i] = std::exp(-alpha * r[i] * r[i]);
    }
    const double omega = 137.0; // arbitrary cell volume for the 4pi/omega scale
    const double gaussian_norm = std::sqrt(pi) / (4.0 * std::pow(alpha, 1.5));

    std::vector<double> scratch(mesh); // one buffer reused across the K values
    for (const double k : {0.0, 0.5, 1.0, 2.0, 4.0})
    {
        const double radial = gaussian_norm * std::exp(-k * k / (4.0 * alpha));
        const double expect = 4.0 * pi / omega * radial;
        const double got = module_charge::core_charge_form_factor(mesh, r.data(), rab.data(), rho_atc.data(), omega, k, scratch);
        // the Simpson sum on this grid converges to ~1e-11 relative; the
        // tolerance leaves architectural margin
        EXPECT_NEAR(got, expect, 1e-9 * std::abs(expect)) << "k = " << k;
    }
}
