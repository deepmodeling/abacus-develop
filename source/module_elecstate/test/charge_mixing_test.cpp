#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#define private public
#include "../module_charge/charge_mixing.h"

// mock function
Magnetism::~Magnetism(){}
Magnetism::Magnetism(){}
int XC_Functional::get_func_type(){return 1;}
#ifdef __LCAO
InfoNonlocal::InfoNonlocal(){}
InfoNonlocal::~InfoNonlocal(){}
#endif
// mock class cell
namespace GlobalC{
  UnitCell ucell;
};

/************************************************
 *  unit test of charge_mixing.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - SetMixingTest: Charge_Mixing::set_mixing(mixing_mode_in,mixing_beta_in,mixing_ndim_in,mixing_gg0_in,mixing_tau_in)
 *      - set the basic parameters of class charge_mixing
 */

class ChargeMixingTest : public ::testing::Test
{
    
};

TEST_F(ChargeMixingTest,SetMixingTest)
{
    Charge_Mixing CMtest;
    std::string mode ="1";
    double beta=1.0;
    int dim=1;
	double gg0=1;
	bool tau= true;
    CMtest.set_mixing(mode,beta,dim,gg0,tau);
    EXPECT_EQ(CMtest.mixing_mode, "1");
    EXPECT_EQ(CMtest.mixing_beta, 1.0);
    EXPECT_EQ(CMtest.mixing_ndim, 1);
    EXPECT_EQ(CMtest.mixing_gg0, 1);
    EXPECT_EQ(CMtest.mixing_tau, true);
}

#undef private
