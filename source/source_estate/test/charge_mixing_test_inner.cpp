// Inner product tests for charge_mixing
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "../module_charge/charge_mixing.h"
#include "source_base/module_mixing/broyden_mixing.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

class ChargeMixingInnerTest : public ::testing::Test
{
  public:
    UnitCell ucell;
    ChargeMixingInnerTest()
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

TEST_F(ChargeMixingInnerTest, InnerDotRealTest)
{
    Charge_Mixing CMtest;
    // non mixing angle case
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
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    PARAM.input.nspin = 4;

    // a simple sum for inner product
    std::vector<double> drho1(pw_basis.nrxx * PARAM.input.nspin);
    std::vector<double> drho2(pw_basis.nrxx * PARAM.input.nspin);
    for (int i = 0; i < pw_basis.nrxx * PARAM.input.nspin; ++i)
    {
        drho1[i] = 1.0;
        drho2[i] = double(i);
    }
    double inner = CMtest.inner_product_real(drho1.data(), drho2.data());
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * PARAM.input.nspin  * (pw_basis.nrxx * PARAM.input.nspin - 1), 1e-8);

    // mixing angle case
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
    PARAM.input.nspin = 4;

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

TEST_F(ChargeMixingInnerTest, InnerDotRecipSimpleTest)
{
    Charge_Mixing CMtest;
    // non mixing angle case
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
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    PARAM.input.nspin = 2;

    // a simple sum for inner product
    std::vector<std::complex<double>> drhog1(pw_basis.npw * PARAM.input.nspin);
    std::vector<std::complex<double>> drhog2(pw_basis.npw * PARAM.input.nspin);
    for (int i = 0; i < pw_basis.npw * PARAM.input.nspin; ++i)
    {
        drhog1[i] = 1.0;
        drhog2[i] = double(i);
    }
    double inner = CMtest.inner_product_recip_simple(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 0.5 * pw_basis.npw * PARAM.input.nspin * (pw_basis.npw * PARAM.input.nspin - 1), 1e-8);
}

TEST_F(ChargeMixingInnerTest, InnerDotRecipHartreeTest)
{
    // REAL
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    const int npw = pw_basis.npw;
    const int nrxx = pw_basis.nrxx;
    PARAM.input.nspin = 1;
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
    ucell.tpiba2 = 1.0;
    ucell.omega = 2.0;
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
    PARAM.input.nspin = 1;
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
    PARAM.input.nspin = 2;
    drhog1.resize(pw_basis.npw * PARAM.input.nspin);
    drhog2.resize(pw_basis.npw * PARAM.input.nspin);
    std::vector<std::complex<double>> drhog1_mag(pw_basis.npw * PARAM.input.nspin);
    std::vector<std::complex<double>> drhog2_mag(pw_basis.npw * PARAM.input.nspin);
    for (int i = 0; i < pw_basis.npw * PARAM.input.nspin; ++i)
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
    PARAM.sys.gamma_only_pw= false;
    inner = CMtest.inner_product_recip_hartree(drhog1_mag.data(), drhog2_mag.data());
    EXPECT_NEAR(inner, 236763.82650318215, 1e-8);
    PARAM.sys.gamma_only_pw= true;
    inner = CMtest.inner_product_recip_hartree(drhog1_mag.data(), drhog2_mag.data());
    EXPECT_NEAR(inner, 236763.82650318215 * 2, 1e-8);

    // RECIPROCAL NSPIN=4 without mixing_angle
    PARAM.input.nspin = 4;
    drhog1.resize(pw_basis.npw * PARAM.input.nspin);
    drhog2.resize(pw_basis.npw * PARAM.input.nspin);
    for (int i = 0; i < pw_basis.npw * PARAM.input.nspin; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }

    PARAM.sys.domag = false;
    PARAM.sys.domag_z = false;
    inner = CMtest.inner_product_recip_hartree(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 28260.091995611871, 1e-8);
    PARAM.sys.gamma_only_pw= true;
    PARAM.sys.domag = true;
    PARAM.sys.domag_z = true;
    inner = CMtest.inner_product_recip_hartree(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 110668.61166927818, 1e-8);

    // RECIPROCAL NSPIN=4 with mixing_angle
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
    drhog1.resize(pw_basis.npw * 2);
    drhog2.resize(pw_basis.npw * 2);
    for (int i = 0; i < pw_basis.npw * 2; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }
    PARAM.sys.gamma_only_pw= false;
    inner = CMtest.inner_product_recip_hartree(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 36548.881431837777, 1e-8);
    PARAM.sys.gamma_only_pw= true;
    inner = CMtest.inner_product_recip_hartree(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 44776.555369916401, 1e-8);
}

TEST_F(ChargeMixingInnerTest, InnerDotRecipRhoTest)
{
    // REAL
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    PARAM.input.nspin = 1;
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
    ucell.tpiba2 = 1.0;
    ucell.omega = 2.0;
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
    PARAM.input.nspin = 1;
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

    PARAM.input.nspin = 2;
    drhog1.resize(pw_basis.npw * PARAM.input.nspin);
    drhog2.resize(pw_basis.npw * PARAM.input.nspin);
    for (int i = 0; i < pw_basis.npw * PARAM.input.nspin; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }
    PARAM.sys.gamma_only_pw= false;
    inner = CMtest.inner_product_recip_rho(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 236763.82650318215, 1e-8);
    PARAM.sys.gamma_only_pw= true;
    inner = CMtest.inner_product_recip_rho(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 236763.82650318215 * 2, 1e-8);

    PARAM.input.nspin = 4;
    drhog1.resize(pw_basis.npw * PARAM.input.nspin);
    drhog2.resize(pw_basis.npw * PARAM.input.nspin);
    for (int i = 0; i < pw_basis.npw * PARAM.input.nspin; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }

    PARAM.sys.domag = false;
    PARAM.sys.domag_z = false;
    inner = CMtest.inner_product_recip_rho(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 28260.091995611871, 1e-8);
    PARAM.sys.gamma_only_pw= true;
    PARAM.sys.domag = true;
    PARAM.sys.domag_z = true;
    inner = CMtest.inner_product_recip_rho(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 110668.61166927818, 1e-8);
}
