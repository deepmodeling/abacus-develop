#include "gtest/gtest.h"

#include "source_cell/unitcell.h"
#include "source_estate/module_charge/chg_tools.h"

#include <algorithm>
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
