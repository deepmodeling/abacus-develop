#include "gtest/gtest.h"

#define private public
#define protected public
#include "module_cell/unitcell.h"
#include "module_elecstate/module_charge/charge.h"
#include "prepare_unitcell.h"

// mock functions for UnitCell
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
#endif
Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = nullptr;
}
Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}

// mock functions for Charge
namespace elecstate
{
int tmp_xc_func_type = 1;
int get_xc_func_type()
{
    return tmp_xc_func_type;
}
double tmp_ucell_omega = 500.0;
double get_ucell_omega()
{
    return tmp_ucell_omega;
}
double tmp_gridecut = 80.0;
void Set_GlobalV_Default()
{
    GlobalV::NSPIN = 1;
    GlobalV::test_charge = 0;
}
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
 *     - using rhopw and GlobalV::NSPIN
 *   - SumRho: Charge::sum_rho()
 *     - calculate \sum_{spin}\sum_{nrxx} rho[is][ir]
 */

class ChargeTest : public ::testing::Test
{
  protected:
    UcellTestPrepare utp = UcellTestLib["Si"];
    std::unique_ptr<UnitCell> ucell;
    Charge* charge;
    ModulePW::PW_Basis* rhopw;
    void SetUp() override
    {
        elecstate::Set_GlobalV_Default();
        ucell = utp.SetUcellInfo();
        charge = new Charge;
        rhopw = new ModulePW::PW_Basis;
    }
    void TearDown() override
    {
        delete charge;
        delete rhopw;
    }
};

TEST_F(ChargeTest, Constructor)
{
    EXPECT_FALSE(charge->allocate_rho);
    EXPECT_FALSE(charge->allocate_rho_final_scf);
}

TEST_F(ChargeTest, Allocate)
{
    // UcellTestPrepare utp = UcellTestLib["Si"];
    // ucell = utp.SetUcellInfo();
    //  init ucell
    EXPECT_DOUBLE_EQ(ucell->omega, 265.302);
    // init rhopw

    rhopw->initgrids(ucell->lat0, ucell->latvec, elecstate::tmp_gridecut);
    EXPECT_DOUBLE_EQ(rhopw->lat0, 10.2);
    EXPECT_EQ(rhopw->nx, 24);
    EXPECT_EQ(rhopw->ny, 24);
    EXPECT_EQ(rhopw->nz, 24);
    EXPECT_EQ(rhopw->nxyz, 13824);
    rhopw->distribute_r();
    EXPECT_EQ(rhopw->nrxx, 13824);
    rhopw->initparameters(false, elecstate::tmp_gridecut);
    rhopw->distribute_g();
    EXPECT_EQ(rhopw->npw, 3143);
    EXPECT_EQ(rhopw->npwtot, 3143);
    // call Charge::allocate()
    charge = new Charge;
    GlobalV::test_charge = 2;
    elecstate::tmp_xc_func_type = 3;
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->allocate_rho);
    charge->allocate(GlobalV::NSPIN);
    EXPECT_TRUE(charge->allocate_rho);
    // test if Charge::allocate() be called twice

    EXPECT_NO_THROW(charge->allocate(GlobalV::NSPIN));
    EXPECT_TRUE(charge->allocate_rho);
    // delete rhopw;
    //charge->destroy();
    /*
    if (charge->rhopw != nullptr)
    {
        delete charge->rhopw;
        charge->rhopw = nullptr;
    }
    delete charge;
    */
}

TEST_F(ChargeTest, SumRho)
{
    EXPECT_DOUBLE_EQ(ucell->omega, 265.302);

    rhopw->initgrids(ucell->lat0, ucell->latvec, elecstate::tmp_gridecut);
    EXPECT_DOUBLE_EQ(rhopw->lat0, 10.2);
    EXPECT_EQ(rhopw->nx, 24);
    EXPECT_EQ(rhopw->ny, 24);
    EXPECT_EQ(rhopw->nz, 24);
    EXPECT_EQ(rhopw->nxyz, 13824);
    rhopw->distribute_r();
    EXPECT_EQ(rhopw->nrxx, 13824);
    rhopw->initparameters(false, elecstate::tmp_gridecut);
    rhopw->distribute_g();
    EXPECT_EQ(rhopw->npw, 3143);
    EXPECT_EQ(rhopw->npwtot, 3143);
    // call Charge::allocate()
    charge = new Charge;
    GlobalV::test_charge = 2;
    elecstate::tmp_xc_func_type = 3;
    charge->set_rhopw(rhopw);
    EXPECT_FALSE(charge->allocate_rho);
    charge->allocate(GlobalV::NSPIN);
    EXPECT_TRUE(charge->allocate_rho);
    // test if Charge::allocate() be called twice

    EXPECT_NO_THROW(charge->allocate(GlobalV::NSPIN));
    EXPECT_TRUE(charge->allocate_rho);
}

#undef protected
#undef private