// Kerker screen tests for charge_mixing
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "../module_charge/charge_mixing.h"
#include "source_base/module_mixing/broyden_mixing.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_io/module_parameter/parameter.h"

class ChargeMixingKerkerTest : public ::testing::Test
{
  public:
    UnitCell ucell;
    ChargeMixingKerkerTest()
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

TEST_F(ChargeMixingKerkerTest, KerkerScreenRecipTest)
{
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    ucell.tpiba = 1.0;
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
    // nspin = 1
    PARAM.input.nspin = 1;
    std::complex<double>* drhog = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    std::complex<double>* drhog_old = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    // no kerker
    CMtest.mixing_gg0 = 0.0;
    CMtest.Kerker_screen_recip(drhog);
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }
    // kerker
    CMtest.mixing_gg0 = 1.0;
    CMtest.Kerker_screen_recip(drhog);
    double gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg0), 0.1 / CMtest.mixing_beta);
        EXPECT_NEAR(drhog[i].real(), ref, 1e-10);
        EXPECT_NEAR(drhog[i].imag(), ref, 1e-10);
    }
    delete[] drhog;
    delete[] drhog_old;

    // nspin = 2
    PARAM.input.nspin = 2;
    CMtest.mixing_beta = 0.4;
    CMtest.mixing_beta_mag = 1.6;
    drhog = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    drhog_old = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    // mixing_gg0 = 0.0
    CMtest.mixing_gg0 = 0.0;
    CMtest.Kerker_screen_recip(drhog);
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }
    // mixing_gg0 = 1.0, mixing_gg0_mag = 0.0
    CMtest.mixing_gg0 = 1.0;
    CMtest.Kerker_screen_recip(drhog);
    gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg0), 0.1 / CMtest.mixing_beta);
        // rho
        EXPECT_NEAR(drhog[i].real(), ref, 1e-10);
        EXPECT_NEAR(drhog[i].imag(), ref, 1e-10);
        // mag
        EXPECT_NEAR(drhog[i+pw_basis.npw].real(), 1.0, 1e-10);
        EXPECT_NEAR(drhog[i+pw_basis.npw].imag(), 1.0, 1e-10);
    }
    delete[] drhog;
    delete[] drhog_old;

    // nspin = 4
    PARAM.input.nspin = 4;
    drhog = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    drhog_old = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    // mixing_gg0 = 0.0
    CMtest.mixing_gg0 = 0.0;
    CMtest.Kerker_screen_recip(drhog);
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }
    // mixing_gg0 = 1.0, mixing_gg0_mag = 0.0
    CMtest.mixing_gg0 = 1.0;
    CMtest.Kerker_screen_recip(drhog);
    gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg0), 0.1 / CMtest.mixing_beta);
        // rho
        EXPECT_NEAR(drhog[i].real(), ref, 1e-10);
        EXPECT_NEAR(drhog[i].imag(), ref, 1e-10);
    }
    for (int i = 0; i < 3*pw_basis.npw; ++i)
    {
        EXPECT_NEAR(drhog[i + pw_basis.npw].real(), 1.0, 1e-10);
        EXPECT_NEAR(drhog[i + pw_basis.npw].imag(), 1.0, 1e-10);
    }
    // mixing_gg0 = 1.0, mixing_gg0_mag = 2.0
    CMtest.mixing_gg0 = 1.0;
    CMtest.mixing_gg0_mag = 2.0;
    CMtest.Kerker_screen_recip(drhog);
    double gg1 = std::pow(1.0 * ModuleBase::BOHR_TO_A, 2);
    double gg2 = std::pow(2.0 * ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg1), 0.1 / CMtest.mixing_beta);
        // rho
        EXPECT_NEAR(drhog[i].real(), ref * ref, 1e-10);
        EXPECT_NEAR(drhog[i].imag(), ref * ref, 1e-10);
    }
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg2), 0.1 / CMtest.mixing_beta_mag);
        // rho
        for (int j = 1; j < PARAM.input.nspin; ++j)
        {
            EXPECT_NEAR(drhog[i + pw_basis.npw * j].real(), ref, 1e-10);
            EXPECT_NEAR(drhog[i + pw_basis.npw * j].imag(), ref, 1e-10);
        }
    }
    delete[] drhog;
    delete[] drhog_old;
}

TEST_F(ChargeMixingKerkerTest, KerkerScreenRealTest)
{
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    ucell.tpiba = 1.0;
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
    // nspin = 1
    PARAM.input.nspin = 1;
    double* drhor = new double[PARAM.input.nspin*pw_basis.nrxx];
    double* drhor_ref = new double[PARAM.input.nspin*pw_basis.nrxx];
    for (int i = 0; i < PARAM.input.nspin*pw_basis.nrxx; ++i)
    {
        drhor_ref[i] = drhor[i] = 1.0;
    }
    // no kerker
    CMtest.mixing_gg0 = 0.0;
    CMtest.Kerker_screen_real(drhor);
    for (int i = 0; i < PARAM.input.nspin*pw_basis.nrxx; ++i)
    {
        EXPECT_EQ(drhor[i], drhor_ref[i]);
    }
    delete[] drhor;
    delete[] drhor_ref;

    // nspin = 2
    PARAM.input.nspin = 2;
    CMtest.mixing_gg0 = 0.0;
    std::complex<double>* drhog = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    std::complex<double>* drhog_old = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    drhor = new double[PARAM.input.nspin*pw_basis.nrxx];
    drhor_ref = new double[PARAM.input.nspin*pw_basis.nrxx];
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    CMtest.Kerker_screen_recip(drhog); // no kerker
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }

    // RECIPROCAL
    CMtest.mixing_gg0 = 1.0;
    PARAM.input.mixing_gg0_mag = 0.0;
    CMtest.Kerker_screen_recip(drhog);
    const double gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
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
    PARAM.input.mixing_gg0_mag = 0.0;
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
