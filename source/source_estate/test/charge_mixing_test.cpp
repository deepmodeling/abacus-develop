#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "../module_charge/charge_mixing.h"
#include "source_base/module_mixing/broyden_mixing.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int XC_Functional::func_type = 1;
bool XC_Functional::ked_flag = false;

// mock function
Magnetism::~Magnetism()
{
}
Magnetism::Magnetism()
{
}
Charge::~Charge()
{
}
Charge::Charge()
{
}

void Charge::set_rhopw(ModulePW::PW_Basis* rhopw_in)
{
    this->rhopw = rhopw_in;
}
#ifdef __LCAO
InfoNonlocal::InfoNonlocal()
{
}
InfoNonlocal::~InfoNonlocal()
{
}
#endif
// mock class cell
/************************************************
 *  unit test of charge_mixing.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - SetMixingTest:
 * Charge_Mixing::set_mixing()
 *                    Charge_Mixing::init_mixing()
 *                    Charge_Mixing::set_rhopw(rhopw_in)
 *                    Charge_Mixing::get_mixing_mode()
 *                    Charge_Mixing::get_mixing_beta()
 *                    Charge_Mixing::get_mixing_ndim()
 *                    Charge_Mixing::get_mixing_gg0()
 *      - set the basic parameters of class charge_mixing
 *   - KerkerScreenTest: Charge_Mixing::Kerker_screen_recip(drhog)
 *                       Charge_Mixing::Kerker_screen_real(drhog)
 *      - screen drho with Kerker method
 *   - InnerDotTest: Charge_Mixing::inner_product_recip_hartree(rhog1, rhog2)
 *                   Charge_Mixing::inner_product_recip_rho(rhog1, rhog2)
 *                   Charge_Mixing::inner_product_recip_simple(rhog1, rhog2)
 *                   Charge_Mixing::inner_product_real(rho1, rho2)
 *      - calculate the inner product of two vectors
 *   - MixRhoTest: Charge_Mixing::mix_rho(chr)
 *                 Charge_Mixing::mix_rho_recip(chr)
 *                 Charge_Mixing::mix_rho_real(chr)
 *      - mix rho with different methods
 *   - MixDivCombTest: Charge_Mixing::divide_data
 *                     Charge_Mixing::combine_data
 *                     Charge_Mixing::clean_data
 *    - divide and combine data
 *
 */

class ChargeMixingTest : public ::testing::Test
{
  public:
    UnitCell ucell;
    ChargeMixingTest()
    {
        // Init pw_basis
        pw_basis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 20);
        pw_basis.initparameters(false, 20);
        pw_basis.setuptransform();
        pw_basis.collect_local_pw();
        pw_dbasis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 40);
        pw_dbasis.initparameters(false, 40);
        pw_dbasis.setuptransform(&pw_basis);
        pw_dbasis.collect_local_pw();
        // default mixing parameters
        PARAM.input.mixing_mode = "broyden";
        PARAM.input.mixing_beta = 0.8;
        PARAM.input.mixing_ndim = 8;
        PARAM.input.mixing_gg0  = 1.0;
        PARAM.input.mixing_tau  = false;
        PARAM.input.mixing_beta_mag = 1.6;
        PARAM.input.mixing_gg0_mag = 0.0;
        PARAM.input.mixing_gg0_min = 0.1;
        PARAM.input.mixing_angle = -10.0;
        PARAM.input.mixing_dmr = false;
        ucell.omega = 1.0;
        ucell.tpiba = 1.0;
    }
    ModulePW::PW_Basis pw_basis;
    ModulePW::PW_Basis_Sup pw_dbasis;
    Charge charge;
};

TEST_F(ChargeMixingTest, SetMixingTest)
{
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    PARAM.input.nspin = 1;
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    PARAM.input.mixing_beta = 1.0;
    PARAM.input.mixing_ndim = 1;
    PARAM.input.mixing_gg0 = 1.0;

    CMtest.set_mixing(PARAM.input.mixing_mode,
                    PARAM.input.mixing_beta,
                    PARAM.input.mixing_ndim,
                    PARAM.input.mixing_gg0,
                    PARAM.input.mixing_tau,
                    PARAM.input.mixing_beta_mag,
                    PARAM.input.mixing_gg0_mag,
                    PARAM.input.mixing_gg0_min,
                    PARAM.input.mixing_angle,
                    PARAM.input.mixing_dmr,
                    ucell.omega,
                    ucell.tpiba);
    EXPECT_EQ(CMtest.get_mixing_mode(), "broyden");
    EXPECT_EQ(CMtest.get_mixing_beta(), 1.0);
    EXPECT_EQ(CMtest.get_mixing_ndim(), 1);
    EXPECT_EQ(CMtest.get_mixing_gg0(), 1.0);
    EXPECT_EQ(CMtest.mixing_tau, false);
    EXPECT_EQ(CMtest.mixing_beta_mag, 1.6);
    EXPECT_EQ(CMtest.mixing_gg0_mag, 0.0);
    EXPECT_EQ(CMtest.mixing_gg0_min, 0.1);
    EXPECT_EQ(CMtest.mixing_angle, -10.0);
    EXPECT_EQ(CMtest.mixing_dmr, false);

    PARAM.input.mixing_tau = true;
    PARAM.input.mixing_mode = "plain";
    CMtest.set_mixing(PARAM.input.mixing_mode,
                    PARAM.input.mixing_beta,
                    PARAM.input.mixing_ndim,
                    PARAM.input.mixing_gg0,
                    PARAM.input.mixing_tau,
                    PARAM.input.mixing_beta_mag,
                    PARAM.input.mixing_gg0_mag,
                    PARAM.input.mixing_gg0_min,
                    PARAM.input.mixing_angle,
                    PARAM.input.mixing_dmr,
                    ucell.omega,
                    ucell.tpiba);
    EXPECT_EQ(CMtest.mixing_mode, "plain");
    EXPECT_EQ(CMtest.mixing_tau, true);

    PARAM.input.mixing_beta = 1.1;
    std::string output;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(PARAM.input.mixing_mode,
                                PARAM.input.mixing_beta,
                                PARAM.input.mixing_ndim,
                                PARAM.input.mixing_gg0,
                                PARAM.input.mixing_tau,
                                PARAM.input.mixing_beta_mag,
                                PARAM.input.mixing_gg0_mag,
                                PARAM.input.mixing_gg0_min,
                                PARAM.input.mixing_angle,
                                PARAM.input.mixing_dmr,
                    ucell.omega,
                    ucell.tpiba);, ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("You'd better set mixing_beta to [0.0, 1.0]!"));

    PARAM.input.mixing_beta = 0.7;
    PARAM.input.mixing_beta_mag = -0.1;
    PARAM.input.nspin = 2;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(PARAM.input.mixing_mode,
                                PARAM.input.mixing_beta,
                                PARAM.input.mixing_ndim,
                                PARAM.input.mixing_gg0,
                                PARAM.input.mixing_tau,
                                PARAM.input.mixing_beta_mag,
                                PARAM.input.mixing_gg0_mag,
                                PARAM.input.mixing_gg0_min,
                                PARAM.input.mixing_angle,
                                PARAM.input.mixing_dmr,
                    ucell.omega,
                    ucell.tpiba);, ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("You'd better set mixing_beta_mag >= 0.0!"));

    PARAM.input.nspin = 1;
    PARAM.input.mixing_beta = 0.7;
    PARAM.input.mixing_beta_mag = 1.6;
    PARAM.input.mixing_mode = "nothing";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(PARAM.input.mixing_mode,
                                PARAM.input.mixing_beta,
                                PARAM.input.mixing_ndim,
                                PARAM.input.mixing_gg0,
                                PARAM.input.mixing_tau,
                                PARAM.input.mixing_beta_mag,
                                PARAM.input.mixing_gg0_mag,
                                PARAM.input.mixing_gg0_min,
                                PARAM.input.mixing_angle,
                                PARAM.input.mixing_dmr,
                    ucell.omega,
                    ucell.tpiba);, ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("This Mixing mode is not implemended yet,coming soon."));
}

TEST_F(ChargeMixingTest, InitMixingTest)
{
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    PARAM.input.nspin = 1;
    XC_Functional::func_type = 1;
    XC_Functional::ked_flag = false;
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);

    CMtest.set_mixing(PARAM.input.mixing_mode,
                    PARAM.input.mixing_beta,
                    PARAM.input.mixing_ndim,
                    PARAM.input.mixing_gg0,
                    PARAM.input.mixing_tau,
                    PARAM.input.mixing_beta_mag,
                    PARAM.input.mixing_gg0_mag,
                    PARAM.input.mixing_gg0_min,
                    PARAM.input.mixing_angle,
                    PARAM.input.mixing_dmr,
                    ucell.omega,
                    ucell.tpiba);

    PARAM.input.scf_thr_type= 1;
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, pw_basis.npw);

    PARAM.input.scf_thr_type= 2;
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, pw_basis.nrxx);

    PARAM.input.nspin = 4;
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, 4 * pw_basis.nrxx);

    PARAM.input.nspin = 1;
    PARAM.input.mixing_tau = true;
    CMtest.set_mixing(PARAM.input.mixing_mode,
                    PARAM.input.mixing_beta,
                    PARAM.input.mixing_ndim,
                    PARAM.input.mixing_gg0,
                    PARAM.input.mixing_tau,
                    PARAM.input.mixing_beta_mag,
                    PARAM.input.mixing_gg0_mag,
                    PARAM.input.mixing_gg0_min,
                    PARAM.input.mixing_angle,
                    PARAM.input.mixing_dmr,
                    ucell.omega,
                    ucell.tpiba);
    XC_Functional::func_type = 3;
    XC_Functional::ked_flag = true;
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.tau_mdata.length, pw_basis.nrxx);

    PARAM.input.nspin = 4;
    PARAM.input.mixing_angle = 1.0;
    CMtest.set_mixing(PARAM.input.mixing_mode,
                    PARAM.input.mixing_beta,
                    PARAM.input.mixing_ndim,
                    PARAM.input.mixing_gg0,
                    PARAM.input.mixing_tau,
                    PARAM.input.mixing_beta_mag,
                    PARAM.input.mixing_gg0_mag,
                    PARAM.input.mixing_gg0_min,
                    PARAM.input.mixing_angle,
                    PARAM.input.mixing_dmr,
                    ucell.omega,
                    ucell.tpiba);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, 2 * pw_basis.nrxx);
}

