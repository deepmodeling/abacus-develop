#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"
#include "prepare_unitcell.h"
// mock functions for UnitCell

Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

// mock functions for Charge
int XC_Functional::func_type = 1;
bool XC_Functional::ked_flag = false;
namespace elecstate
{
double tmp_ucell_omega = 500.0;
double tmp_gridecut = 80.0;
} // namespace elecstate

/// The Charge methods under test take nspin, test_charge and nelec as
/// arguments, so the fixture owns them instead of writing the global singleton.
namespace
{
constexpr int test_nspin = 1;
constexpr int test_charge_verbosity = 0;
constexpr double test_nelec = 8.0;
} // namespace

/************************************************
 *  unit test of module_charge/charge.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Constructor: Charge::Charge() and Charge::~Charge()
 *     - this is a trivial test
 *   - Allocate: Charge::set_rhopw(), Charge::allocate(), Charge::destroy()
 *     - allocate rho, rhog, rho_save, rhog_save, kin_r, kin_r_save
 *     - using rhopw and PARAM.input.nspin
 *   - SumRho: Charge::sum_rho()
 *     - calculate \sum_{is}^nspin \sum_{ir}^nrxx rho[is][ir]
 *   - RenormalizeRho: Charge::renormalize_rho()
 *     - renormalize rho so as to ensure the sum of rho equals to total number of electrons
 *   - CheckNe: Charge::cal_rho2ne()
 *     - check the total number of electrons summed from rho[is]
 *   - SaveRhoBeforeSumBand: Charge::save_rho_before_sum_band()
 *     - meaning as the function name
 *   - InitFinalScf:: Charge::init_final_scf()
 *     - similar to Charge::allocate(), but for final scf
 */

class ChargeTest : public ::testing::Test
{
  protected:
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::unique_ptr<UnitCell> ucell;
    Charge* charge;
    ModulePW::PW_Basis* rhopw;
    std::string output;
    void SetUp() override
    {
        ucell = utp.SetUcellInfo();
        charge = new Charge;
        rhopw = new ModulePW::PW_Basis;
        rhopw->initgrids(ucell->lat0, ucell->latvec, elecstate::tmp_gridecut);
        rhopw->initparameters(false, elecstate::tmp_gridecut);
        // the documented public sequence (see pw_basis.h): setuptransform()
        // performs distribute_r() and distribute_g() itself
        rhopw->setuptransform();
    }
    void TearDown() override
    {
        delete charge;
        delete rhopw;
    }
};

TEST_F(ChargeTest, Constructor)
{
    EXPECT_FALSE(charge->get_allocate_rho());
    EXPECT_FALSE(charge->get_allocate_rho_final_scf());
}

TEST_F(ChargeTest, Allocate)
{
    // ucell info
    EXPECT_DOUBLE_EQ(ucell->omega, 265.302);
    // rhopw info
    EXPECT_DOUBLE_EQ(rhopw->lat0, 10.2);
    EXPECT_EQ(rhopw->nx, 24);
    EXPECT_EQ(rhopw->ny, 24);
    EXPECT_EQ(rhopw->nz, 24);
    EXPECT_EQ(rhopw->nxyz, 13824);
    EXPECT_EQ(rhopw->nrxx, 13824);
    EXPECT_EQ(rhopw->npw, 3143);
    EXPECT_EQ(rhopw->npwtot, 3143);
    // call Charge::allocate() with the verbose diagnostics branch enabled
    XC_Functional::set_func_type(3);
    XC_Functional::set_ked_flag(true);
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->get_allocate_rho());
    const bool kin_den = charge->kin_density();
    charge->allocate(test_nspin, kin_den, /*test_charge=*/2, test_nelec);
    EXPECT_TRUE(charge->get_allocate_rho());
    // test if Charge::allocate() be called twice
    EXPECT_NO_THROW(charge->allocate(test_nspin, kin_den, /*test_charge=*/2, test_nelec));
    EXPECT_TRUE(charge->get_allocate_rho());
}

TEST_F(ChargeTest, SumRho)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->get_allocate_rho());
    const bool kin_den = charge->kin_density();
    charge->allocate(test_nspin, kin_den, test_charge_verbosity, test_nelec);
    EXPECT_TRUE(charge->get_allocate_rho());
    int nspin = (test_nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    charge->set_omega(&ucell->omega);;
    EXPECT_NEAR(charge->sum_rho(), 0.1 * nspin * rhopw->nrxx * ucell->omega / rhopw->nxyz, 1E-10);
}

TEST_F(ChargeTest, RenormalizeRho)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->get_allocate_rho());
    const bool kin_den = charge->kin_density();
    charge->allocate(test_nspin, kin_den, test_charge_verbosity, test_nelec);
    EXPECT_TRUE(charge->get_allocate_rho());
    int nspin = (test_nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    // guards the 8.0 the assertions below expect from renormalize_rho()
    EXPECT_DOUBLE_EQ(test_nelec, 8.0);
    charge->set_omega(&ucell->omega);;
    charge->renormalize_rho();
    EXPECT_NEAR(charge->sum_rho(), 8.0, 1e-10);
}

TEST_F(ChargeTest, CheckNe)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->get_allocate_rho());
    const bool kin_den = charge->kin_density();
    charge->allocate(test_nspin, kin_den, test_charge_verbosity, test_nelec);
    EXPECT_TRUE(charge->get_allocate_rho());
    int nspin = (test_nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    // guards the 8.0 the assertions below expect from renormalize_rho()
    EXPECT_DOUBLE_EQ(test_nelec, 8.0);
    charge->set_omega(&ucell->omega);;
    charge->renormalize_rho();
    EXPECT_NEAR(charge->sum_rho(), 8.0, 1e-10);
    EXPECT_NEAR(charge->cal_rho2ne(charge->rho[0]), 8.0, 1e-10);
}

TEST_F(ChargeTest, SaveRhoBeforeSumBand)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->get_allocate_rho());
    const bool kin_den = charge->kin_density();
    charge->allocate(test_nspin, kin_den, test_charge_verbosity, test_nelec);
    EXPECT_TRUE(charge->get_allocate_rho());
    int nspin = (test_nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    // guards the 8.0 the assertions below expect from renormalize_rho()
    EXPECT_DOUBLE_EQ(test_nelec, 8.0);
    XC_Functional::set_func_type(3);
    XC_Functional::set_ked_flag(true);
    charge->set_omega(&ucell->omega);;
    charge->renormalize_rho();
    charge->save_rho_before_sum_band(test_nspin);
    EXPECT_NEAR(charge->cal_rho2ne(charge->rho_save[0]), 8.0, 1e-10);
}

TEST_F(ChargeTest, InitFinalScf)
{
    charge->set_rhopw(rhopw);
    XC_Functional::set_func_type(1);
    XC_Functional::set_ked_flag(false);
    charge->init_final_scf(test_nspin, /*test_charge=*/2);
    EXPECT_TRUE(charge->get_allocate_rho_final_scf());
}

