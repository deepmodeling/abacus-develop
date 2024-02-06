#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/module_mixing/broyden_mixing.h"
#define private public
#include "../module_charge/charge_mixing.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_hamilt_general/module_xc/xc_functional.h"
#include <omp.h>

int FUNC_TYPE = 1;

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
int XC_Functional::get_func_type()
{
    return FUNC_TYPE;
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
namespace GlobalC
{
UnitCell ucell;
};
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
        GlobalV::MIXING_MODE = "broyden";
        GlobalV::MIXING_BETA = 0.8;
        GlobalV::MIXING_NDIM = 8;
        GlobalV::MIXING_GG0  = 1.0;
        GlobalV::MIXING_TAU  = false;
        GlobalV::MIXING_BETA_MAG = 1.6;
        GlobalV::MIXING_GG0_MAG = 0.0;
        GlobalV::MIXING_GG0_MIN = 0.1;
        GlobalV::MIXING_ANGLE = -10.0;
        GlobalV::MIXING_DMR = false;
    }
    ModulePW::PW_Basis pw_basis;
    ModulePW::PW_Basis_Sup pw_dbasis;
    Charge charge;    
};

TEST_F(ChargeMixingTest, SetMixingTest)
{
    omp_set_num_threads(1);
    GlobalV::NSPIN = 1;
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    GlobalV::MIXING_BETA = 1.0;
    GlobalV::MIXING_NDIM = 1;
    GlobalV::MIXING_GG0 = 1.0;

    CMtest.set_mixing(GlobalV::MIXING_MODE,
                    GlobalV::MIXING_BETA,
                    GlobalV::MIXING_NDIM,
                    GlobalV::MIXING_GG0,
                    GlobalV::MIXING_TAU,
                    GlobalV::MIXING_BETA_MAG,
                    GlobalV::MIXING_GG0_MAG,
                    GlobalV::MIXING_GG0_MIN,
                    GlobalV::MIXING_ANGLE,
                    GlobalV::MIXING_DMR);
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

    GlobalV::MIXING_TAU = true;
    GlobalV::MIXING_MODE = "plain";
    CMtest.set_mixing(GlobalV::MIXING_MODE,
                    GlobalV::MIXING_BETA,
                    GlobalV::MIXING_NDIM,
                    GlobalV::MIXING_GG0,
                    GlobalV::MIXING_TAU,
                    GlobalV::MIXING_BETA_MAG,
                    GlobalV::MIXING_GG0_MAG,
                    GlobalV::MIXING_GG0_MIN,
                    GlobalV::MIXING_ANGLE,
                    GlobalV::MIXING_DMR);
    EXPECT_EQ(CMtest.mixing_mode, "plain");
    EXPECT_EQ(CMtest.mixing_tau, true);

    GlobalV::MIXING_BETA = 1.1;
    std::string output;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(GlobalV::MIXING_MODE,
                                GlobalV::MIXING_BETA,
                                GlobalV::MIXING_NDIM,
                                GlobalV::MIXING_GG0,
                                GlobalV::MIXING_TAU,
                                GlobalV::MIXING_BETA_MAG,
                                GlobalV::MIXING_GG0_MAG,
                                GlobalV::MIXING_GG0_MIN,
                                GlobalV::MIXING_ANGLE,
                                GlobalV::MIXING_DMR);, ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("You'd better set mixing_beta to [0.0, 1.0]!"));

    GlobalV::MIXING_BETA = 0.7;
    GlobalV::MIXING_BETA_MAG = -0.1;
    GlobalV::NSPIN = 2;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(GlobalV::MIXING_MODE,
                                GlobalV::MIXING_BETA,
                                GlobalV::MIXING_NDIM,
                                GlobalV::MIXING_GG0,
                                GlobalV::MIXING_TAU,
                                GlobalV::MIXING_BETA_MAG,
                                GlobalV::MIXING_GG0_MAG,
                                GlobalV::MIXING_GG0_MIN,
                                GlobalV::MIXING_ANGLE,
                                GlobalV::MIXING_DMR);, ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("You'd better set mixing_beta_mag >= 0.0!"));

    GlobalV::NSPIN = 1;
    GlobalV::MIXING_BETA = 0.7;
    GlobalV::MIXING_BETA_MAG = 1.6;
    GlobalV::MIXING_MODE = "nothing";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(GlobalV::MIXING_MODE,
                                GlobalV::MIXING_BETA,
                                GlobalV::MIXING_NDIM,
                                GlobalV::MIXING_GG0,
                                GlobalV::MIXING_TAU,
                                GlobalV::MIXING_BETA_MAG,
                                GlobalV::MIXING_GG0_MAG,
                                GlobalV::MIXING_GG0_MIN,
                                GlobalV::MIXING_ANGLE,
                                GlobalV::MIXING_DMR);, ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("This Mixing mode is not implemended yet,coming soon."));
}

TEST_F(ChargeMixingTest, InitMixingTest)
{
    omp_set_num_threads(1);
    GlobalV::NSPIN = 1;
    FUNC_TYPE = 1;
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);

    CMtest.set_mixing(GlobalV::MIXING_MODE,
                    GlobalV::MIXING_BETA,
                    GlobalV::MIXING_NDIM,
                    GlobalV::MIXING_GG0,
                    GlobalV::MIXING_TAU,
                    GlobalV::MIXING_BETA_MAG,
                    GlobalV::MIXING_GG0_MAG,
                    GlobalV::MIXING_GG0_MIN,
                    GlobalV::MIXING_ANGLE,
                    GlobalV::MIXING_DMR);
    
    GlobalV::SCF_THR_TYPE = 1;
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, pw_basis.npw);
    
    GlobalV::SCF_THR_TYPE = 2;
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, pw_basis.nrxx);

    GlobalV::NSPIN = 4;
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, 4 * pw_basis.nrxx);

    GlobalV::NSPIN = 1;
    GlobalV::MIXING_TAU = true;
    CMtest.set_mixing(GlobalV::MIXING_MODE,
                    GlobalV::MIXING_BETA,
                    GlobalV::MIXING_NDIM,
                    GlobalV::MIXING_GG0,
                    GlobalV::MIXING_TAU,
                    GlobalV::MIXING_BETA_MAG,
                    GlobalV::MIXING_GG0_MAG,
                    GlobalV::MIXING_GG0_MIN,
                    GlobalV::MIXING_ANGLE,
                    GlobalV::MIXING_DMR);
    FUNC_TYPE = 3;
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.tau_mdata.length, pw_basis.nrxx);

    GlobalV::NSPIN = 4;
    GlobalV::MIXING_ANGLE = 1.0;
    CMtest.set_mixing(GlobalV::MIXING_MODE,
                    GlobalV::MIXING_BETA,
                    GlobalV::MIXING_NDIM,
                    GlobalV::MIXING_GG0,
                    GlobalV::MIXING_TAU,
                    GlobalV::MIXING_BETA_MAG,
                    GlobalV::MIXING_GG0_MAG,
                    GlobalV::MIXING_GG0_MIN,
                    GlobalV::MIXING_ANGLE,
                    GlobalV::MIXING_DMR);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, 2 * pw_basis.nrxx);
}

TEST_F(ChargeMixingTest, InnerDotRealTest)
{
    Charge_Mixing CMtest;
    // non mixing angle case
    CMtest.set_mixing(GlobalV::MIXING_MODE,
                    GlobalV::MIXING_BETA,
                    GlobalV::MIXING_NDIM,
                    GlobalV::MIXING_GG0,
                    GlobalV::MIXING_TAU,
                    GlobalV::MIXING_BETA_MAG,
                    GlobalV::MIXING_GG0_MAG,
                    GlobalV::MIXING_GG0_MIN,
                    GlobalV::MIXING_ANGLE,
                    GlobalV::MIXING_DMR);
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    GlobalV::NSPIN = 4;

    // a simple sum for inner product
    std::vector<double> drho1(pw_basis.nrxx * GlobalV::NSPIN);
    std::vector<double> drho2(pw_basis.nrxx * GlobalV::NSPIN);
    for (int i = 0; i < pw_basis.nrxx * GlobalV::NSPIN; ++i)
    {
        drho1[i] = 1.0;
        drho2[i] = double(i);
    }
    double inner = CMtest.inner_product_real(drho1.data(), drho2.data());
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * GlobalV::NSPIN  * (pw_basis.nrxx * GlobalV::NSPIN - 1), 1e-8);

    // mixing angle case
    GlobalV::MIXING_ANGLE = 1.0;
    CMtest.set_mixing(GlobalV::MIXING_MODE,
                    GlobalV::MIXING_BETA,
                    GlobalV::MIXING_NDIM,
                    GlobalV::MIXING_GG0,
                    GlobalV::MIXING_TAU,
                    GlobalV::MIXING_BETA_MAG,
                    GlobalV::MIXING_GG0_MAG,
                    GlobalV::MIXING_GG0_MIN,
                    GlobalV::MIXING_ANGLE,
                    GlobalV::MIXING_DMR);
    GlobalV::NSPIN = 4;

    // a simple sum for inner product
    drho1.resize(pw_basis.nrxx * 2);
    drho2.resize(pw_basis.nrxx * 2);
    for (int i = 0; i < pw_basis.nrxx * 2; ++i)
    {
        drho1[i] = 1.0;
        drho2[i] = double(i);
    }
    inner = CMtest.inner_product_real(drho1.data(), drho2.data());
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * 2  * (pw_basis.nrxx * 2 - 1), 1e-8);
}

TEST_F(ChargeMixingTest, InnerDotRecipSimpleTest)
{
    Charge_Mixing CMtest;
    // non mixing angle case
    CMtest.set_mixing(GlobalV::MIXING_MODE,
                    GlobalV::MIXING_BETA,
                    GlobalV::MIXING_NDIM,
                    GlobalV::MIXING_GG0,
                    GlobalV::MIXING_TAU,
                    GlobalV::MIXING_BETA_MAG,
                    GlobalV::MIXING_GG0_MAG,
                    GlobalV::MIXING_GG0_MIN,
                    GlobalV::MIXING_ANGLE,
                    GlobalV::MIXING_DMR);
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    GlobalV::NSPIN = 2;

    // a simple sum for inner product
    std::vector<std::complex<double>> drhog1(pw_basis.npw * GlobalV::NSPIN);
    std::vector<std::complex<double>> drhog2(pw_basis.npw * GlobalV::NSPIN);
    for (int i = 0; i < pw_basis.npw * GlobalV::NSPIN; ++i)
    {
        drhog1[i] = 1.0;
        drhog2[i] = double(i);
    }
    double inner = CMtest.inner_product_recip_simple(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 0.5 * pw_basis.npw * GlobalV::NSPIN * (pw_basis.npw * GlobalV::NSPIN - 1), 1e-8);
}

TEST_F(ChargeMixingTest, InnerDotRecipHartreeTest)
{
    // REAL
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    const int npw = pw_basis.npw;
    const int nrxx = pw_basis.nrxx;
    GlobalV::NSPIN = 1;
    std::vector<double> drhor1(pw_basis.nrxx);
    std::vector<double> drhor2(pw_basis.nrxx);
    for (int i = 0; i < pw_basis.nrxx; ++i)
    {
        drhor1[i] = 1.0;
        drhor2[i] = double(i);
    }
    double inner = CMtest.inner_product_real(drhor1.data(), drhor2.data());
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * (pw_basis.nrxx - 1), 1e-8);

    // RECIPROCAL NSPIN=1
    GlobalC::ucell.tpiba2 = 1.0;
    GlobalC::ucell.omega = 2.0;

    GlobalV::NSPIN = 1;
    std::vector<std::complex<double>> drhog1(pw_basis.npw);
    std::vector<std::complex<double>> drhog2(pw_basis.npw);
    for (int i = 0; i < pw_basis.nrxx; ++i)
    {
        drhor1[i] = 0.0;
    }
    drhor1[2] = 1.0;
    pw_basis.real2recip(drhor1.data(), drhog1.data());
    pw_basis.real2recip(drhor2.data(), drhog2.data());

    inner = CMtest.inner_product_recip_hartree(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, -0.3 * ModuleBase::e2 * ModuleBase::FOUR_PI, 1e-8);

    // RECIPROCAL NSPIN=2
    GlobalV::NSPIN = 2;
    drhog1.resize(pw_basis.npw * GlobalV::NSPIN);
    drhog2.resize(pw_basis.npw * GlobalV::NSPIN);
    std::vector<std::complex<double>> drhog1_mag(pw_basis.npw * GlobalV::NSPIN);
    std::vector<std::complex<double>> drhog2_mag(pw_basis.npw * GlobalV::NSPIN);
    for (int i = 0; i < pw_basis.npw * GlobalV::NSPIN; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }
    // set mag
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        drhog1_mag[i] = drhog1[i] + drhog1[i+pw_basis.npw];
        drhog1_mag[i+pw_basis.npw] = drhog1[i] - drhog1[i+pw_basis.npw];
        drhog2_mag[i] = drhog2[i] + drhog2[i+pw_basis.npw];
        drhog2_mag[i+pw_basis.npw] = drhog2[i] - drhog2[i+pw_basis.npw];
    }
    GlobalV::GAMMA_ONLY_PW = false;
    inner = CMtest.inner_product_recip_hartree(drhog1_mag.data(), drhog2_mag.data());
    EXPECT_NEAR(inner, 236763.82650318215, 1e-8);
    GlobalV::GAMMA_ONLY_PW = true;
    inner = CMtest.inner_product_recip_hartree(drhog1_mag.data(), drhog2_mag.data());
    EXPECT_NEAR(inner, 236763.82650318215 * 2, 1e-8);

    // RECIPROCAL NSPIN=4 without mixing_angle
    GlobalV::NSPIN = 4;
    drhog1.resize(pw_basis.npw * GlobalV::NSPIN);
    drhog2.resize(pw_basis.npw * GlobalV::NSPIN);
    for (int i = 0; i < pw_basis.npw * GlobalV::NSPIN; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }

    GlobalV::DOMAG = false;
    GlobalV::DOMAG_Z = false;
    inner = CMtest.inner_product_recip_hartree(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 28260.091995611871, 1e-8);
    GlobalV::GAMMA_ONLY_PW = true;
    GlobalV::DOMAG = true;
    GlobalV::DOMAG_Z = true;
    inner = CMtest.inner_product_recip_hartree(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 110668.61166927818, 1e-8);

    // RECIPROCAL NSPIN=4 with mixing_angle
    GlobalV::NSPIN = 4;
    GlobalV::MIXING_ANGLE = 1.0;
    CMtest.set_mixing(GlobalV::MIXING_MODE,
                    GlobalV::MIXING_BETA,
                    GlobalV::MIXING_NDIM,
                    GlobalV::MIXING_GG0,
                    GlobalV::MIXING_TAU,
                    GlobalV::MIXING_BETA_MAG,
                    GlobalV::MIXING_GG0_MAG,
                    GlobalV::MIXING_GG0_MIN,
                    GlobalV::MIXING_ANGLE,
                    GlobalV::MIXING_DMR);
    drhog1.resize(pw_basis.npw * 2);
    drhog2.resize(pw_basis.npw * 2);
    for (int i = 0; i < pw_basis.npw * 2; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }
    GlobalV::GAMMA_ONLY_PW = false;
    inner = CMtest.inner_product_recip_hartree(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 36548.881431837777, 1e-8);
    GlobalV::GAMMA_ONLY_PW = true;
    inner = CMtest.inner_product_recip_hartree(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 44776.555369916401, 1e-8);
}

TEST_F(ChargeMixingTest, InnerDotRecipRhoTest)
{
    // REAL
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    GlobalV::NSPIN = 1;
    std::vector<double> drhor1(pw_basis.nrxx);
    std::vector<double> drhor2(pw_basis.nrxx);
    for (int i = 0; i < pw_basis.nrxx; ++i)
    {
        drhor1[i] = 1.0;
        drhor2[i] = double(i);
    }
    double inner = CMtest.inner_product_real(drhor1.data(), drhor2.data());
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * (pw_basis.nrxx - 1), 1e-8);

    // RECIPROCAL
    GlobalC::ucell.tpiba2 = 1.0;
    GlobalC::ucell.omega = 2.0;

    GlobalV::NSPIN = 1;
    std::vector<std::complex<double>> drhog1(pw_basis.npw);
    std::vector<std::complex<double>> drhog2(pw_basis.npw);
    for (int i = 0; i < pw_basis.nrxx; ++i)
    {
        drhor1[i] = 0.0;
    }
    drhor1[2] = 1.0;
    pw_basis.real2recip(drhor1.data(), drhog1.data());
    pw_basis.real2recip(drhor2.data(), drhog2.data());

    inner = CMtest.inner_product_recip_rho(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, -0.3 * ModuleBase::e2 * ModuleBase::FOUR_PI, 1e-8);

    GlobalV::NSPIN = 2;
    drhog1.resize(pw_basis.npw * GlobalV::NSPIN);
    drhog2.resize(pw_basis.npw * GlobalV::NSPIN);
    for (int i = 0; i < pw_basis.npw * GlobalV::NSPIN; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }
    GlobalV::GAMMA_ONLY_PW = false;
    inner = CMtest.inner_product_recip_rho(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 236763.82650318215, 1e-8);
    GlobalV::GAMMA_ONLY_PW = true;
    inner = CMtest.inner_product_recip_rho(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 236763.82650318215 * 2, 1e-8);

    GlobalV::NSPIN = 4;
    drhog1.resize(pw_basis.npw * GlobalV::NSPIN);
    drhog2.resize(pw_basis.npw * GlobalV::NSPIN);
    for (int i = 0; i < pw_basis.npw * GlobalV::NSPIN; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }

    GlobalV::DOMAG = false;
    GlobalV::DOMAG_Z = false;
    inner = CMtest.inner_product_recip_rho(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 28260.091995611871, 1e-8);
    GlobalV::GAMMA_ONLY_PW = true;
    GlobalV::DOMAG = true;
    GlobalV::DOMAG_Z = true;
    inner = CMtest.inner_product_recip_rho(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 110668.61166927818, 1e-8);
}

TEST_F(ChargeMixingTest, KerkerScreenRecipNewTest)
{
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    GlobalC::ucell.tpiba = 1.0;

    // for new kerker
    GlobalV::NSPIN = 2;
    CMtest.mixing_gg0 = 0.0;
    GlobalV::MIXING_GG0_MAG = 0.0;
    std::complex<double>* drhog = new std::complex<double>[GlobalV::NSPIN*pw_basis.npw];
    std::complex<double>* drhog_old = new std::complex<double>[GlobalV::NSPIN*pw_basis.npw];
    double* drhor = new double[GlobalV::NSPIN*pw_basis.nrxx];
    double* drhor_ref = new double[GlobalV::NSPIN*pw_basis.nrxx];
    for (int i = 0; i < GlobalV::NSPIN*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    CMtest.Kerker_screen_recip_new(drhog);
    for (int i = 0; i < GlobalV::NSPIN*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }

    // RECIPROCAL
    CMtest.mixing_gg0 = 1.0;
    GlobalV::MIXING_GG0_MAG = 0.0;
    CMtest.Kerker_screen_recip_new(drhog);
    const double gg0 = std::pow(0.529177, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        std::complex<double> ration = drhog[i] / drhog[i+pw_basis.npw];
        double gg = this->pw_basis.gg[i];
        double ration_ref = std::max(gg / (gg + gg0), 0.1 / CMtest.mixing_beta);
        EXPECT_NEAR(ration.real(), ration_ref, 1e-10);
        EXPECT_NEAR(ration.imag(), 0, 1e-10);
    }

    // REAL
    pw_basis.recip2real(drhog, drhor_ref);
    pw_basis.recip2real(drhog_old, drhor);

    CMtest.mixing_gg0 = 0.0;
    GlobalV::MIXING_GG0_MAG = 0.0;
    // nothing happens
    CMtest.Kerker_screen_real(drhor);

    CMtest.mixing_gg0 = 1.0;
    CMtest.Kerker_screen_real(drhor);
    for (int i = 0; i < pw_basis.nrxx; ++i)
    {
        EXPECT_NEAR(drhor[i], drhor_ref[i], 1e-8);
    }

    delete[] drhog;
    delete[] drhog_old;
    delete[] drhor;
    delete[] drhor_ref;
}

TEST_F(ChargeMixingTest, MixRhoTest)
{
    GlobalV::double_grid = false;
    charge.set_rhopw(&pw_basis);
    const int nspin = GlobalV::NSPIN = 1;
    GlobalV::DOMAG_Z = false;
    FUNC_TYPE = 3;
    GlobalV::MIXING_BETA = 0.7;
    GlobalV::MIXING_NDIM = 1;
    GlobalV::MIXING_GG0 = 0.0;
    GlobalV::MIXING_TAU = true;
    GlobalV::MIXING_MODE = "plain";
    const int nrxx = pw_basis.nrxx;
    const int npw = pw_basis.npw;
    charge._space_rho = new double[nspin * nrxx];
    charge._space_rho_save = new double[nspin * nrxx];
    charge._space_rhog = new std::complex<double>[nspin * npw];
    charge._space_rhog_save = new std::complex<double>[nspin * npw];
    charge._space_kin_r = new double[nspin * nrxx];
    charge._space_kin_r_save = new double[nspin * nrxx];
    charge.rho = new double*[nspin];
    charge.rhog = new std::complex<double>*[nspin];
    charge.rho_save = new double*[nspin];
    charge.rhog_save = new std::complex<double>*[nspin];
    charge.kin_r = new double*[nspin];
    charge.kin_r_save = new double*[nspin];
    for (int is = 0; is < nspin; is++)
    {
        charge.rho[is] = charge._space_rho + is * nrxx;
        charge.rhog[is] = charge._space_rhog + is * npw;
        charge.rho_save[is] = charge._space_rho_save + is * nrxx;
        charge.rhog_save[is] = charge._space_rhog_save + is * npw;
        charge.kin_r[is] = charge._space_kin_r + is * nrxx;
        charge.kin_r_save[is] = charge._space_kin_r_save + is * nrxx;
    }
    std::vector<double> real_ref(nspin * nrxx);
    std::vector<double> real_save_ref(nspin * nrxx);
    std::vector<std::complex<double>> recip_ref(nspin * npw);
    std::vector<std::complex<double>> recip_save_ref(nspin * npw);
    for(int i = 0 ; i < nspin * npw; ++i)
    {
       recip_ref[i] = std::complex<double>(double(i), 1.0);
       recip_save_ref[i] = std::complex<double>(double(i), 0.0);
    }
    for(int i = 0 ; i < nspin ; ++i)
    {
        pw_basis.recip2real(recip_ref.data() + i * npw, real_ref.data() + i * nrxx);
        pw_basis.recip2real(recip_save_ref.data() + i * npw, real_save_ref.data() + i * nrxx);
    }
    //--------------------------------MAIN BODY--------------------------------
    // RECIPROCAL
    Charge_Mixing CMtest_recip;
    CMtest_recip.set_rhopw(&pw_basis, &pw_basis);
    GlobalV::SCF_THR_TYPE = 1;
    CMtest_recip.set_mixing(GlobalV::MIXING_MODE,
                            GlobalV::MIXING_BETA,
                            GlobalV::MIXING_NDIM,
                            GlobalV::MIXING_GG0,
                            GlobalV::MIXING_TAU,
                            GlobalV::MIXING_BETA_MAG,
                            GlobalV::MIXING_GG0_MAG,
                            GlobalV::MIXING_GG0_MIN,
                            GlobalV::MIXING_ANGLE,
                            GlobalV::MIXING_DMR);
    CMtest_recip.init_mixing();
    for(int i = 0 ; i < nspin * npw; ++i)
    {
        charge._space_rhog[i] = recip_ref[i];
        charge._space_rhog_save[i] = recip_save_ref[i];
    }
    for(int i = 0 ; i < nspin * nrxx; ++i)
    {
        charge._space_rho[i] = real_ref[i];
        charge._space_rho_save[i] = real_save_ref[i];
    }
    CMtest_recip.mix_rho(&charge);
    for(int is = 0 ; is < nspin; ++is)
    {
        for(int ir = 0 ; ir < nrxx ; ++ir)
        {
            EXPECT_NEAR(charge.rho_save[is][ir], real_ref[is*nrxx + ir], 1e-8);
        }
        for(int ig = 0; ig < npw ; ++ig)
        {
            EXPECT_NEAR(charge.rhog[is][ig].real(), recip_save_ref[is*npw + ig].real(), 1e-8);
            EXPECT_NEAR(charge.rhog[is][ig].imag(), recip_save_ref[is*npw + ig].imag() + 0.7, 1e-8);
        }
    }

    // REAL
    Charge_Mixing CMtest_real;
    GlobalV::SCF_THR_TYPE = 2;
    CMtest_real.set_rhopw(&pw_basis, &pw_basis);
    CMtest_real.set_mixing(GlobalV::MIXING_MODE,
                        GlobalV::MIXING_BETA,
                        GlobalV::MIXING_NDIM,
                        GlobalV::MIXING_GG0,
                        GlobalV::MIXING_TAU,
                        GlobalV::MIXING_BETA_MAG,
                        GlobalV::MIXING_GG0_MAG,
                        GlobalV::MIXING_GG0_MIN,
                        GlobalV::MIXING_ANGLE,
                        GlobalV::MIXING_DMR);
    CMtest_real.init_mixing();
    for(int i = 0 ; i < nspin * nrxx; ++i)
    {
        charge._space_rho[i] = real_ref[i];
        charge._space_rho_save[i] = real_save_ref[i];
    }
    CMtest_recip.mix_rho(&charge);
    for(int is = 0 ; is < nspin; ++is)
    {
        for(int ir = 0 ; ir < nrxx ; ++ir)
        {
            EXPECT_NEAR(charge.rho_save[is][ir], real_ref[is*nrxx + ir], 1e-8);
            EXPECT_NEAR(charge.rho[is][ir], 0.3*real_save_ref[is*nrxx+ir] + 0.7*real_ref[is*nrxx+ir], 1e-8);
        }
    }

    //-------------------------------------------------------------------------
    delete[] charge._space_rho;
    delete[] charge._space_rho_save;
    delete[] charge._space_rhog;
    delete[] charge._space_rhog_save;
    delete[] charge._space_kin_r;
    delete[] charge._space_kin_r_save;
    delete[] charge.rho;
    delete[] charge.rhog;
    delete[] charge.rho_save;
    delete[] charge.rhog_save;
    delete[] charge.kin_r;
    delete[] charge.kin_r_save;
}

TEST_F(ChargeMixingTest, MixDoubleGridRhoTest)
{
    GlobalV::double_grid = true;
    charge.set_rhopw(&pw_dbasis);
    const int nspin = GlobalV::NSPIN = 1;
    GlobalV::DOMAG_Z = false;
    FUNC_TYPE = 3;
    GlobalV::MIXING_BETA = 0.7;
    GlobalV::MIXING_NDIM = 1;
    GlobalV::MIXING_GG0 = 0.0;
    GlobalV::MIXING_TAU = true;
    GlobalV::MIXING_MODE = "plain";
    const int nrxx = pw_dbasis.nrxx;
    const int npw = pw_dbasis.npw;
    charge._space_rho = new double[nspin * nrxx];
    charge._space_rho_save = new double[nspin * nrxx];
    charge._space_rhog = new std::complex<double>[nspin * npw];
    charge._space_rhog_save = new std::complex<double>[nspin * npw];
    charge._space_kin_r = new double[nspin * nrxx];
    charge._space_kin_r_save = new double[nspin * nrxx];
    charge.rho = new double*[nspin];
    charge.rhog = new std::complex<double>*[nspin];
    charge.rho_save = new double*[nspin];
    charge.rhog_save = new std::complex<double>*[nspin];
    charge.kin_r = new double*[nspin];
    charge.kin_r_save = new double*[nspin];
    for (int is = 0; is < nspin; is++)
    {
        charge.rho[is] = charge._space_rho + is * nrxx;
        charge.rhog[is] = charge._space_rhog + is * npw;
        charge.rho_save[is] = charge._space_rho_save + is * nrxx;
        charge.rhog_save[is] = charge._space_rhog_save + is * npw;
        charge.kin_r[is] = charge._space_kin_r + is * nrxx;
        charge.kin_r_save[is] = charge._space_kin_r_save + is * nrxx;
    }
    std::vector<double> real_ref(nspin * nrxx);
    std::vector<double> real_save_ref(nspin * nrxx);
    std::vector<std::complex<double>> recip_ref(nspin * npw);
    std::vector<std::complex<double>> recip_save_ref(nspin * npw);
    for (int i = 0; i < nspin * npw; ++i)
    {
        recip_ref[i] = std::complex<double>(double(i), 1.0);
        recip_save_ref[i] = std::complex<double>(double(i), 0.0);
    }
    for (int i = 0; i < nspin; ++i)
    {
        pw_dbasis.recip2real(recip_ref.data() + i * npw, real_ref.data() + i * nrxx);
        pw_dbasis.recip2real(recip_save_ref.data() + i * npw, real_save_ref.data() + i * nrxx);
    }
    //--------------------------------MAIN BODY--------------------------------
    // RECIPROCAL
    Charge_Mixing CMtest_recip;
    CMtest_recip.set_rhopw(&pw_basis, &pw_dbasis);
    GlobalV::SCF_THR_TYPE = 1;
    CMtest_recip.set_mixing(GlobalV::MIXING_MODE,
                            GlobalV::MIXING_BETA,
                            GlobalV::MIXING_NDIM,
                            GlobalV::MIXING_GG0,
                            GlobalV::MIXING_TAU,
                            GlobalV::MIXING_BETA_MAG,
                            GlobalV::MIXING_GG0_MAG,
                            GlobalV::MIXING_GG0_MIN,
                            GlobalV::MIXING_ANGLE,
                            GlobalV::MIXING_DMR);
    CMtest_recip.init_mixing();
    for (int i = 0; i < nspin * npw; ++i)
    {
        charge._space_rhog[i] = recip_ref[i];
        charge._space_rhog_save[i] = recip_save_ref[i];
    }
    for (int i = 0; i < nspin * nrxx; ++i)
    {
        charge._space_rho[i] = real_ref[i];
        charge._space_rho_save[i] = real_save_ref[i];
    }
    CMtest_recip.mix_rho(&charge);
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            EXPECT_NEAR(charge.rho_save[is][ir], real_ref[is * nrxx + ir], 1e-8);
        }
        for (int ig = 0; ig < npw; ++ig)
        {
            EXPECT_NEAR(charge.rhog[is][ig].real(), recip_save_ref[is * npw + ig].real(), 1e-8);
            EXPECT_NEAR(charge.rhog[is][ig].imag(), recip_save_ref[is * npw + ig].imag() + 0.7, 1e-8);
        }
    }

    //-------------------------------------------------------------------------
    delete[] charge._space_rho;
    delete[] charge._space_rho_save;
    delete[] charge._space_rhog;
    delete[] charge._space_rhog_save;
    delete[] charge._space_kin_r;
    delete[] charge._space_kin_r_save;
    delete[] charge.rho;
    delete[] charge.rhog;
    delete[] charge.rho_save;
    delete[] charge.rhog_save;
    delete[] charge.kin_r;
    delete[] charge.kin_r_save;
}

TEST_F(ChargeMixingTest, MixDivCombTest)
{
    // NSPIN = 1
    GlobalV::NSPIN = 1;
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_dbasis);
    std::vector<std::complex<double>> data(pw_dbasis.npw, 1.0);
    std::complex<double>*datas, *datahf;
    std::complex<double>*datas2, *datahf2;
    CMtest.divide_data(data.data(), datas, datahf);
    EXPECT_EQ(datas, data.data());
    EXPECT_EQ(datahf, data.data() + pw_basis.npw);
    CMtest.combine_data(data.data(), datas, datahf);
    EXPECT_EQ(datas, nullptr);
    EXPECT_EQ(datahf, nullptr);

    CMtest.divide_data(data.data(), datas2, datahf2);
    CMtest.clean_data(datas2, datahf2);
    EXPECT_EQ(datas2, nullptr);
    EXPECT_EQ(datahf2, nullptr);

    // NSPIN = 2
    GlobalV::NSPIN = 2;
    data.resize(pw_dbasis.npw * 2, 1.0);
    std::vector<std::complex<double>> dataout(pw_dbasis.npw * 2, 1.0);
    CMtest.divide_data(data.data(), datas, datahf);
    CMtest.combine_data(dataout.data(), datas, datahf);
    EXPECT_EQ(datas, nullptr);
    EXPECT_EQ(datahf, nullptr);
    for (int i = 0; i < pw_dbasis.npw * 2; ++i)
    {
        EXPECT_EQ(dataout[i], data[i]);
    }

    CMtest.divide_data(data.data(), datas2, datahf2);
    CMtest.clean_data(datas2, datahf2);
    EXPECT_EQ(datas2, nullptr);
    EXPECT_EQ(datahf2, nullptr);
}

#undef private
