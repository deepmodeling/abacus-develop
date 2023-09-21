#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_base/module_mixing/broyden_mixing.h"
#define private public
#include "../module_charge/charge_mixing.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_hamilt_general/module_xc/xc_functional.h"

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
 * Charge_Mixing::set_mixing(mixing_mode_in,mixing_beta_in,mixing_ndim_in,mixing_gg0_in,mixing_tau_in)
 *                    Charge_Mixing::set_rhopw(rhopw_in)
 *                    Charge_Mixing::get_mixing_mode()
 *                    Charge_Mixing::get_mixing_beta()
 *                    Charge_Mixing::get_mixing_ndim()
 *                    Charge_Mixing::get_mixing_gg0()
 *      - set the basic parameters of class charge_mixing
 *   - AutoSetTest: Charge_Mixing::auto_set(bandgap_in, ucell_)
 *                  Charge_Mixing::need_auto_set()
 *      - auto-set the parameters of class charge_mixing
 *   - KerkerScreenTest: Charge_Mixing::Kerker_screen_recip(drhog)
 *                       Charge_Mixing::Kerker_screen_real(drhog)
 *      - screen drho with Kerker method
 *   - InnerDotTest: Charge_Mixing::inner_product_recip(rhog1, rhog2)
 *                   Charge_Mixing::rhog_dot_product(rhog1, rhog2)
 *                   Charge_Mixing::inner_product_real(rho1, rho2)
 *      - calculate the inner product of two vectors
 *   - MixRhoTest: Charge_Mixing::mix_rho(iter, chr)
 *                 Charge_Mixing::mix_rho_recip(iter, chr)
 *                 Charge_Mixing::mix_rho_real(iter, chr)
 *      - mix rho with different methods
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
    }
    ModulePW::PW_Basis pw_basis;
    Charge charge;
};

TEST_F(ChargeMixingTest, SetMixingTest)
{
    GlobalV::NSPIN = 1;
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis);
    double beta = 1.0;
    int dim = 1;
    double gg0 = 1;

    FUNC_TYPE = 1;
    bool mixingtau = false;
    GlobalV::SCF_THR_TYPE = 1;
    std::string mode = "broyden";
    CMtest.set_mixing(mode, beta, dim, gg0, mixingtau);
    EXPECT_EQ(CMtest.rho_mdata.length, pw_basis.npw);

    GlobalV::SCF_THR_TYPE = 2;
    mode = "broyden";
    CMtest.set_mixing(mode, beta, dim, gg0, mixingtau);
    EXPECT_EQ(CMtest.rho_mdata.length, pw_basis.nrxx);
    EXPECT_EQ(CMtest.get_mixing_mode(), "broyden");
    EXPECT_EQ(CMtest.get_mixing_beta(), 1.0);
    EXPECT_EQ(CMtest.get_mixing_ndim(), 1);
    EXPECT_EQ(CMtest.get_mixing_gg0(), 1.0);

    FUNC_TYPE = 3;
    mixingtau = true;
    mode = "plain";
    GlobalV::SCF_THR_TYPE = 1;
    CMtest.set_mixing(mode, beta, dim, gg0, mixingtau);
    EXPECT_EQ(CMtest.tau_mdata.length, pw_basis.npw);

    GlobalV::SCF_THR_TYPE = 2;
    CMtest.set_mixing(mode, beta, dim, gg0, mixingtau);
    EXPECT_EQ(CMtest.tau_mdata.length, pw_basis.nrxx);

    mode = "nothing";
    std::string output;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(mode, beta, dim, gg0, mixingtau);, ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("This Mixing mode is not implemended yet,coming soon."));
}

TEST_F(ChargeMixingTest, AutoSetTest)
{
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis);
    GlobalV::SCF_THR_TYPE = 1;
    GlobalV::NSPIN = 1;

    CMtest.set_mixing("broyden", 1.0, 1, 0.2, false);
    CMtest.auto_set(0.0, GlobalC::ucell);
    EXPECT_EQ(CMtest.mixing_beta, 1.0);
    EXPECT_EQ(CMtest.mixing_gg0, 0.2);

    CMtest.need_auto_set();
    CMtest.auto_set(0.0, GlobalC::ucell);
    EXPECT_EQ(CMtest.mixing_beta, 0.2);
    EXPECT_EQ(CMtest.mixing->mixing_beta, 0.2);
    EXPECT_EQ(CMtest.mixing_gg0, 0.0);

    CMtest.need_auto_set();
    CMtest.auto_set(1.0, GlobalC::ucell);
    EXPECT_EQ(CMtest.mixing_beta, 0.7);
    EXPECT_EQ(CMtest.mixing->mixing_beta, 0.7);
    EXPECT_EQ(CMtest.mixing_gg0, 0.0);

    GlobalC::ucell.atoms = new Atom[1];
    GlobalC::ucell.ntype = 1;
    GlobalC::ucell.atoms[0].ncpp.psd = "Sc";
    CMtest.need_auto_set();
    CMtest.auto_set(1.0, GlobalC::ucell);
    EXPECT_EQ(CMtest.mixing_gg0, 1.5);
}

TEST_F(ChargeMixingTest, KerkerScreenRecipTest)
{
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis);
    GlobalV::NSPIN = 1;
    GlobalC::ucell.tpiba = 1.0;

    CMtest.mixing_gg0 = 0.0;
    std::complex<double>* drhog = new std::complex<double>[pw_basis.npw];
    std::complex<double>* drhog_old = new std::complex<double>[pw_basis.npw];
    double* drhor = new double[pw_basis.nrxx];
    double* drhor_ref = new double[pw_basis.nrxx];
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    CMtest.Kerker_screen_recip(drhog);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }

    // RECIPROCAL
    CMtest.mixing_gg0 = 1.0;
    CMtest.Kerker_screen_recip(drhog);
    const double gg0 = std::pow(0.529177, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        std::complex<double> ration = drhog[i] / drhog_old[i];
        double gg = this->pw_basis.gg[i];
        double ration_ref = std::max(gg / (gg + gg0), 0.1);
        EXPECT_NEAR(ration.real(), ration_ref, 1e-10);
        EXPECT_NEAR(ration.imag(), 0, 1e-10);
    }

    // REAL
    pw_basis.recip2real(drhog, drhor_ref);
    pw_basis.recip2real(drhog_old, drhor);

    CMtest.mixing_gg0 = 0.0;
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

TEST_F(ChargeMixingTest, InnerDotTest)
{
    // REAL
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis);
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

    inner = CMtest.inner_product_recip(drhog1.data(), drhog2.data());
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
    inner = CMtest.inner_product_recip(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 236763.82650318215, 1e-8);
    GlobalV::GAMMA_ONLY_PW = true;
    inner = CMtest.inner_product_recip(drhog1.data(), drhog2.data());
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
    inner = CMtest.inner_product_recip(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 28260.091995611871, 1e-8);
    GlobalV::GAMMA_ONLY_PW = true;
    GlobalV::DOMAG = true;
    GlobalV::DOMAG_Z = true;
    inner = CMtest.inner_product_recip(drhog1.data(), drhog2.data());
    EXPECT_NEAR(inner, 110668.61166927818, 1e-8);
}

TEST_F(ChargeMixingTest, MixRhoTest)
{
    const int nspin = GlobalV::NSPIN = 1;
    GlobalV::DOMAG_Z = false;
    FUNC_TYPE = 3;
    double beta = 0.7;
    int dim = 1;
    double gg0 = 0.0;
    bool mixingtau = true;
    std::string mode = "plain";
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
    CMtest_recip.set_rhopw(&pw_basis);
    GlobalV::SCF_THR_TYPE = 1;
    CMtest_recip.set_mixing(mode, beta, dim, gg0, mixingtau);
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
    CMtest_recip.mix_rho(1, &charge);
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
    CMtest_real.set_rhopw(&pw_basis);
    GlobalV::SCF_THR_TYPE = 2;
    CMtest_real.set_mixing(mode, beta, dim, gg0, mixingtau);
    for(int i = 0 ; i < nspin * nrxx; ++i)
    {
        charge._space_rho[i] = real_ref[i];
        charge._space_rho_save[i] = real_save_ref[i];
    }
    CMtest_recip.mix_rho(1, &charge);
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
#undef private
