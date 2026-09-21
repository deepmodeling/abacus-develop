#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "source_cell/unitcell.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_tools.h"
#include "source_hamilt/module_xc/xc_functional.h"
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
// xc_functional.cpp is not linked into this target, so the private statics
// need a definition here. Defining them out of line does not require access
// to the class, only changing them does - that goes through the setters.
int XC_Functional::func_type = 1;
bool XC_Functional::ked_flag = false;
namespace elecstate
{
double tmp_ucell_omega = 500.0;
double tmp_gridecut = 80.0;
} // namespace elecstate

/************************************************
 *  unit test of module_charge/charge.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Constructor: Charge::Charge() and Charge::~Charge()
 *     - this is a trivial test
 *   - Allocate: Charge::set_rhopw(), Charge::allocate(), Charge::destroy()
 *     - allocate rho, rhog, rho_save, rhog_save, kin_r, kin_r_save
 *     - using rhopw and the nspin passed in
 *   - SumRho: Charge::sum_rho()
 *     - calculate \sum_{is}^nspin \sum_{ir}^nrxx rho[is][ir]
 *   - RenormalizeRho: Charge::renormalize_rho()
 *     - renormalize rho so as to ensure the sum of rho equals to total number of electrons
 *   - CheckNe: module_charge::cal_rho2ne()
 *     - check the total number of electrons summed from rho[is]
 *   - SaveRhoBeforeSumBand: Charge::save_rho_before_sum_band()
 *     - meaning as the function name
 */

class ChargeTest : public ::testing::Test
{
  protected:
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::unique_ptr<UnitCell> ucell;
    Charge* charge;
    ModulePW::PW_Basis* rhopw;
    std::string output;
    /// Charge::allocate() and Charge::renormalize_rho() take these explicitly,
    /// so the fixture owns them instead of writing the global parameter
    /// singleton. The values mirror the Input_para defaults the test relied on.
    int nspin = 1;
    int test_charge = 0;
    double nelec = 8;
    bool out_elf_on = false;
    void SetUp() override
    {
        ucell = utp.SetUcellInfo();
        charge = new Charge;
        rhopw = new ModulePW::PW_Basis;
        rhopw->initgrids(ucell->lat0, ucell->latvec, elecstate::tmp_gridecut);
        rhopw->initparameters(false, elecstate::tmp_gridecut);
        // setuptransform() runs distribute_r() then distribute_g(), both of
        // which are protected; this is the public route to the same state.
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
    // call Charge::allocate()
    test_charge = 2;
    XC_Functional::set_func_type(3);
    XC_Functional::set_ked_flag(true);
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->get_allocate_rho());
    const bool kin_den = XC_Functional::get_ked_flag() || out_elf_on;
    charge->allocate(nspin, kin_den, XC_Functional::get_ked_flag(), test_charge);
    EXPECT_TRUE(charge->get_allocate_rho());
    // test if Charge::allocate() be called twice
    EXPECT_NO_THROW(charge->allocate(nspin, kin_den, XC_Functional::get_ked_flag(),
                                     test_charge));
    EXPECT_TRUE(charge->get_allocate_rho());
}

TEST_F(ChargeTest, SumRho)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->get_allocate_rho());
    const bool kin_den = XC_Functional::get_ked_flag() || out_elf_on;
    charge->allocate(nspin, kin_den, XC_Functional::get_ked_flag(), test_charge);
    EXPECT_TRUE(charge->get_allocate_rho());
    int nspin_rho = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin_rho; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    EXPECT_NEAR(charge->sum_rho(ucell->omega), 0.1 * nspin_rho * rhopw->nrxx * ucell->omega / rhopw->nxyz, 1E-10);
}

TEST_F(ChargeTest, RenormalizeRho)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->get_allocate_rho());
    const bool kin_den = XC_Functional::get_ked_flag() || out_elf_on;
    charge->allocate(nspin, kin_den, XC_Functional::get_ked_flag(), test_charge);
    EXPECT_TRUE(charge->get_allocate_rho());
    int nspin_rho = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin_rho; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    EXPECT_EQ(nelec, 8);
    charge->renormalize_rho(nelec, ucell->omega);
    EXPECT_NEAR(charge->sum_rho(ucell->omega), 8.0, 1e-10);
}

TEST_F(ChargeTest, CheckNe)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->get_allocate_rho());
    const bool kin_den = XC_Functional::get_ked_flag() || out_elf_on;
    charge->allocate(nspin, kin_den, XC_Functional::get_ked_flag(), test_charge);
    EXPECT_TRUE(charge->get_allocate_rho());
    int nspin_rho = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin_rho; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    EXPECT_EQ(nelec, 8);
    charge->renormalize_rho(nelec, ucell->omega);
    EXPECT_NEAR(charge->sum_rho(ucell->omega), 8.0, 1e-10);
    EXPECT_NEAR(module_charge::cal_rho2ne(charge->rho[0], rhopw->nrxx, ucell->omega, rhopw->nxyz),
                8.0, 1e-10);
}

TEST_F(ChargeTest, SaveRhoBeforeSumBand)
{
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->get_allocate_rho());
    const bool kin_den = XC_Functional::get_ked_flag() || out_elf_on;
    charge->allocate(nspin, kin_den, XC_Functional::get_ked_flag(), test_charge);
    EXPECT_TRUE(charge->get_allocate_rho());
    int nspin_rho = (nspin == 2) ? 2 : 1;
    for (int is = 0; is < nspin_rho; is++)
    {
        for (int ir = 0; ir < rhopw->nrxx; ir++)
        {
            charge->rho[is][ir] = 0.1;
        }
    }
    EXPECT_EQ(nelec, 8);
    XC_Functional::set_func_type(3);
    XC_Functional::set_ked_flag(true);
    charge->renormalize_rho(nelec, ucell->omega);
    charge->save_rho_before_sum_band();
    EXPECT_NEAR(module_charge::cal_rho2ne(charge->rho_save[0], rhopw->nrxx, ucell->omega, rhopw->nxyz),
                8.0, 1e-10);
}
