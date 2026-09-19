#include "gmock/gmock.h"
#include "gtest/gtest.h"
// TODO(governance): remove this access hack once chg_drho.cpp and Charge no
// longer read global PARAM/RAW state (Step 4 of the module_charge refactor).
// The test still has to drive Charge::_space_* and XC_Functional privates.
#define private public
#include "../chg_mix.h"
#include "../chg_drho.h"
#include "../chg_drho_detail.h"
#include "../chg_precond.h"
#include "../chg_uspp.h"
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

// mock class cell
/************************************************
 *  unit test of chg_mix.cpp
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
 *   - KerkerScreenTest: module_charge::kerker_screen_recip(cfg, rhopw, tpiba, drhog)
 *                       module_charge::kerker_screen_real(cfg, rhopw, tpiba, drhog)
 *      - screen drho with Kerker method
 *   - InnerDotTest: module_charge::inner_product_recip_hartree(rhog1, rhog2)
 *                   module_charge::detail::inner_product_recip_rho(rhog1, rhog2)
 *                   module_charge::inner_product_real(rho1, rho2)
 *      - calculate the inner product of two vectors
 *   - MixRhoTest: Charge_Mixing::mix_rho(chr)
 *                 Charge_Mixing::mix_rho_recip(chr)
 *                 Charge_Mixing::mix_rho_real(chr)
 *      - mix rho with different methods
 *   - MixDivCombTest: module_charge::split_dgrid
 *                     module_charge::merge_dgrid
 *    - divide and combine data on the USPP double grid
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

    // Build a MixingConfig from the PARAM.input values set in the ctor, so
    // set_mixing is driven by explicit config instead of a 12-arg call.
    MixingConfig make_cfg()
    {
        MixingConfig cfg;
        cfg.mixing_mode = PARAM.input.mixing_mode;
        cfg.mixing_beta = PARAM.input.mixing_beta;
        cfg.mixing_ndim = PARAM.input.mixing_ndim;
        cfg.mixing_gg0 = PARAM.input.mixing_gg0;
        // Mirror the esolver-side resolution: tau mixing requires a
        // kinetic-energy-density functional.
        cfg.mixing_tau = PARAM.input.mixing_tau && XC_Functional::get_ked_flag();
        cfg.mixing_beta_mag = PARAM.input.mixing_beta_mag;
        cfg.mixing_gg0_mag = PARAM.input.mixing_gg0_mag;
        cfg.mixing_gg0_min = PARAM.input.mixing_gg0_min;
        cfg.mixing_angle = PARAM.input.mixing_angle;
        cfg.mixing_dmr = PARAM.input.mixing_dmr;
        cfg.nspin = PARAM.input.nspin;
        cfg.scf_thr_type = PARAM.input.scf_thr_type;
        cfg.double_grid = PARAM.globalv.double_grid;
        cfg.gamma_only_pw = PARAM.globalv.gamma_only_pw;
        cfg.domag = PARAM.globalv.domag;
        cfg.domag_z = PARAM.globalv.domag_z;
        cfg.scf_nmax = PARAM.input.scf_nmax;
        return cfg;
    }

    // Re-sync the runtime globals (nspin/scf_thr_type/gamma_only_pw/domag/domag_z)
    // into an already-configured Charge_Mixing. Tests mutate PARAM.sys/PARAM.input
    // after set_mixing to steer the residual/inner-product branches; the object
    // now reads them from cfg_, so the test must push the new values in.
    void sync_cfg(Charge_Mixing& cm)
    {
        cm.cfg_.nspin = PARAM.input.nspin;
        cm.cfg_.scf_thr_type = PARAM.input.scf_thr_type;
        cm.cfg_.gamma_only_pw = PARAM.sys.gamma_only_pw;
        cm.cfg_.domag = PARAM.sys.domag;
        cm.cfg_.domag_z = PARAM.sys.domag_z;
    }
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

    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    EXPECT_EQ(CMtest.get_mixing_mode(), "broyden");
    EXPECT_EQ(CMtest.get_mixing_beta(), 1.0);
    EXPECT_EQ(CMtest.get_mixing_ndim(), 1);
    EXPECT_EQ(CMtest.get_mixing_gg0(), 1.0);
    EXPECT_EQ(CMtest.get_mixing_config().mixing_tau, false);
    EXPECT_EQ(CMtest.mixing_beta_mag, 1.6);
    EXPECT_EQ(CMtest.mixing_gg0_mag, 0.0);
    EXPECT_EQ(CMtest.mixing_gg0_min, 0.1);
    EXPECT_EQ(CMtest.mixing_angle, -10.0);
    EXPECT_EQ(CMtest.mixing_dmr, false);

    PARAM.input.mixing_tau = true;
    XC_Functional::ked_flag = true;
    PARAM.input.mixing_mode = "plain";
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    EXPECT_EQ(CMtest.mixing_mode, "plain");
    EXPECT_EQ(CMtest.get_mixing_config().mixing_tau, true);
    XC_Functional::ked_flag = false;

    PARAM.input.mixing_beta = 1.1;
    std::string output;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);, ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("You'd better set mixing_beta to [0.0, 1.0]!"));

    PARAM.input.mixing_beta = 0.7;
    PARAM.input.mixing_beta_mag = -0.1;
    PARAM.input.nspin = 2;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);, ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("You'd better set mixing_beta_mag >= 0.0!"));

    PARAM.input.nspin = 1;
    PARAM.input.mixing_beta = 0.7;
    PARAM.input.mixing_beta_mag = 1.6;
    PARAM.input.mixing_mode = "nothing";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);, ::testing::ExitedWithCode(1), "");
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

    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    
    PARAM.input.scf_thr_type= 1;
    sync_cfg(CMtest);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, pw_basis.npw);
    
    PARAM.input.scf_thr_type= 2;
    sync_cfg(CMtest);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, pw_basis.nrxx);

    PARAM.input.nspin = 4;
    sync_cfg(CMtest);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, 4 * pw_basis.nrxx);

    PARAM.input.nspin = 1;
    PARAM.input.mixing_tau = true;
    XC_Functional::func_type = 3;
    XC_Functional::ked_flag = true;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.tau_mdata.length, pw_basis.nrxx);

    PARAM.input.nspin = 4;
    PARAM.input.mixing_angle = 1.0;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.rho_mdata.length, 2 * pw_basis.nrxx);
}

TEST_F(ChargeMixingTest, InnerDotRealTest)
{
    Charge_Mixing CMtest;
    // non mixing angle case
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    PARAM.input.nspin = 4;
    sync_cfg(CMtest);

    // a simple sum for inner product
    std::vector<double> drho1(pw_basis.nrxx * PARAM.input.nspin);
    std::vector<double> drho2(pw_basis.nrxx * PARAM.input.nspin);
    for (int i = 0; i < pw_basis.nrxx * PARAM.input.nspin; ++i)
    {
        drho1[i] = 1.0;
        drho2[i] = double(i);
    }
    double inner = module_charge::inner_product_real(drho1.data(), drho2.data(), pw_basis, CMtest.cfg_);
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * PARAM.input.nspin  * (pw_basis.nrxx * PARAM.input.nspin - 1), 1e-8);

    // mixing angle case
    PARAM.input.mixing_angle = 1.0;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    PARAM.input.nspin = 4;

    // a simple sum for inner product
    drho1.resize(pw_basis.nrxx * 2);
    drho2.resize(pw_basis.nrxx * 2);
    for (int i = 0; i < pw_basis.nrxx * 2; ++i)
    {
        drho1[i] = 1.0;
        drho2[i] = double(i);
    }
    inner = module_charge::inner_product_real(drho1.data(), drho2.data(), pw_basis, CMtest.cfg_);
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * 2  * (pw_basis.nrxx * 2 - 1), 1e-8);
}

TEST_F(ChargeMixingTest, InnerDotRecipHartreeTest)
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
    double inner = module_charge::inner_product_real(drhor1.data(), drhor2.data(), pw_basis, CMtest.cfg_);
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * (pw_basis.nrxx - 1), 1e-8);

    // RECIPROCAL NSPIN=1
    ucell.tpiba2 = 1.0;
    ucell.omega = 2.0;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    PARAM.input.nspin = 1;
    sync_cfg(CMtest);
    std::vector<std::complex<double>> drhog1(pw_basis.npw);
    std::vector<std::complex<double>> drhog2(pw_basis.npw);
    for (int i = 0; i < pw_basis.nrxx; ++i)
    {
        drhor1[i] = 0.0;
    }
    drhor1[2] = 1.0;
    pw_basis.real2recip(drhor1.data(), drhog1.data());
    pw_basis.real2recip(drhor2.data(), drhog2.data());

    inner = module_charge::inner_product_recip_hartree(drhog1.data(), drhog2.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, -0.3 * ModuleBase::e2 * ModuleBase::FOUR_PI, 1e-8);

    // RECIPROCAL NSPIN=2
    PARAM.input.nspin = 2;
    sync_cfg(CMtest);
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
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1_mag.data(), drhog2_mag.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 236763.82650318215, 1e-8);
    PARAM.sys.gamma_only_pw= true;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1_mag.data(), drhog2_mag.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 236763.82650318215 * 2, 1e-8);

    // RECIPROCAL NSPIN=4 without mixing_angle
    PARAM.input.nspin = 4;
    sync_cfg(CMtest);
    drhog1.resize(pw_basis.npw * PARAM.input.nspin);
    drhog2.resize(pw_basis.npw * PARAM.input.nspin);
    for (int i = 0; i < pw_basis.npw * PARAM.input.nspin; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }

    PARAM.sys.domag = false;
    PARAM.sys.domag_z = false;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1.data(), drhog2.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 28260.091995611871, 1e-8);
    PARAM.sys.gamma_only_pw= true;
    PARAM.sys.domag = true;
    PARAM.sys.domag_z = true;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1.data(), drhog2.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 110668.61166927818, 1e-8);

    // RECIPROCAL NSPIN=4 with mixing_angle
    PARAM.input.nspin = 4;
    PARAM.input.mixing_angle = 1.0;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    drhog1.resize(pw_basis.npw * 2);
    drhog2.resize(pw_basis.npw * 2);
    for (int i = 0; i < pw_basis.npw * 2; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }
    PARAM.sys.gamma_only_pw= false;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1.data(), drhog2.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 36548.881431837777, 1e-8);
    PARAM.sys.gamma_only_pw= true;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1.data(), drhog2.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 44776.555369916401, 1e-8);
}

TEST_F(ChargeMixingTest, InnerDotRecipRhoTest)
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
    double inner = module_charge::inner_product_real(drhor1.data(), drhor2.data(), pw_basis, CMtest.cfg_);
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * (pw_basis.nrxx - 1), 1e-8);

    // RECIPROCAL
    ucell.tpiba2 = 1.0;
    ucell.omega = 2.0;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    PARAM.input.nspin = 1;
    sync_cfg(CMtest);
    std::vector<std::complex<double>> drhog1(pw_basis.npw);
    std::vector<std::complex<double>> drhog2(pw_basis.npw);
    for (int i = 0; i < pw_basis.nrxx; ++i)
    {
        drhor1[i] = 0.0;
    }
    drhor1[2] = 1.0;
    pw_basis.real2recip(drhor1.data(), drhog1.data());
    pw_basis.real2recip(drhor2.data(), drhog2.data());

    inner = module_charge::detail::inner_product_recip_rho(drhog1.data(), drhog2.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, -0.3 * ModuleBase::e2 * ModuleBase::FOUR_PI, 1e-8);

    PARAM.input.nspin = 2;
    sync_cfg(CMtest);
    drhog1.resize(pw_basis.npw * PARAM.input.nspin);
    drhog2.resize(pw_basis.npw * PARAM.input.nspin);
    for (int i = 0; i < pw_basis.npw * PARAM.input.nspin; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }
    PARAM.sys.gamma_only_pw= false;
    sync_cfg(CMtest);
    inner = module_charge::detail::inner_product_recip_rho(drhog1.data(), drhog2.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 236763.82650318215, 1e-8);
    PARAM.sys.gamma_only_pw= true;
    sync_cfg(CMtest);
    inner = module_charge::detail::inner_product_recip_rho(drhog1.data(), drhog2.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 236763.82650318215 * 2, 1e-8);

    PARAM.input.nspin = 4;
    sync_cfg(CMtest);
    drhog1.resize(pw_basis.npw * PARAM.input.nspin);
    drhog2.resize(pw_basis.npw * PARAM.input.nspin);
    for (int i = 0; i < pw_basis.npw * PARAM.input.nspin; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }

    PARAM.sys.domag = false;
    PARAM.sys.domag_z = false;
    sync_cfg(CMtest);
    inner = module_charge::detail::inner_product_recip_rho(drhog1.data(), drhog2.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 28260.091995611871, 1e-8);
    PARAM.sys.gamma_only_pw= true;
    PARAM.sys.domag = true;
    PARAM.sys.domag_z = true;
    sync_cfg(CMtest);
    inner = module_charge::detail::inner_product_recip_rho(drhog1.data(), drhog2.data(), pw_basis, CMtest.cfg_, ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 110668.61166927818, 1e-8);
}

TEST_F(ChargeMixingTest, KerkerScreenRecipTest)
{
    ucell.tpiba = 1.0;
    // nspin = 1
    PARAM.input.nspin = 1;
    MixingConfig cfg = make_cfg();
    std::complex<double>* drhog = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    std::complex<double>* drhog_old = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    // no kerker
    cfg.mixing_gg0 = 0.0;
    module_charge::kerker_screen_recip(cfg, &pw_basis, ucell.tpiba, drhog);
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }
    // kerker
    cfg.mixing_gg0 = 1.0;
    module_charge::kerker_screen_recip(cfg, &pw_basis, ucell.tpiba, drhog);
    double gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg0), 0.1 / cfg.mixing_beta);
        EXPECT_NEAR(drhog[i].real(), ref, 1e-10);
        EXPECT_NEAR(drhog[i].imag(), ref, 1e-10);
    }
    delete[] drhog;
    delete[] drhog_old;

    // nspin = 2
    PARAM.input.nspin = 2;
    cfg = make_cfg();
    cfg.mixing_beta = 0.4;
    cfg.mixing_beta_mag = 1.6;
    drhog = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    drhog_old = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    // mixing_gg0 = 0.0
    cfg.mixing_gg0 = 0.0;
    module_charge::kerker_screen_recip(cfg, &pw_basis, ucell.tpiba, drhog);
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }
    // mixing_gg0 = 1.0, mixing_gg0_mag = 0.0
    cfg.mixing_gg0 = 1.0;
    module_charge::kerker_screen_recip(cfg, &pw_basis, ucell.tpiba, drhog);
    gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg0), 0.1 / cfg.mixing_beta);
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
    cfg = make_cfg();
    drhog = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    drhog_old = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    // mixing_gg0 = 0.0
    cfg.mixing_gg0 = 0.0;
    module_charge::kerker_screen_recip(cfg, &pw_basis, ucell.tpiba, drhog);
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }
    // mixing_gg0 = 1.0, mixing_gg0_mag = 0.0
    cfg.mixing_gg0 = 1.0;
    module_charge::kerker_screen_recip(cfg, &pw_basis, ucell.tpiba, drhog);
    gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg0), 0.1 / cfg.mixing_beta);
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
    cfg.mixing_gg0 = 1.0;
    cfg.mixing_gg0_mag = 2.0;
    module_charge::kerker_screen_recip(cfg, &pw_basis, ucell.tpiba, drhog);
    double gg1 = std::pow(1.0 * ModuleBase::BOHR_TO_A, 2);
    double gg2 = std::pow(2.0 * ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg1), 0.1 / cfg.mixing_beta);
        // rho
        EXPECT_NEAR(drhog[i].real(), ref * ref, 1e-10);
        EXPECT_NEAR(drhog[i].imag(), ref * ref, 1e-10);
    }
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg2), 0.1 / cfg.mixing_beta_mag);
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

TEST_F(ChargeMixingTest, KerkerScreenRealTest)
{
    ucell.tpiba = 1.0;
    // nspin = 1
    PARAM.input.nspin = 1;
    MixingConfig cfg = make_cfg();
    double* drhor = new double[PARAM.input.nspin*pw_basis.nrxx];
    double* drhor_ref = new double[PARAM.input.nspin*pw_basis.nrxx];
    for (int i = 0; i < PARAM.input.nspin*pw_basis.nrxx; ++i)
    {
        drhor_ref[i] = drhor[i] = 1.0;
    }
    // no kerker
    cfg.mixing_gg0 = 0.0;
    module_charge::kerker_screen_real(cfg, &pw_basis, ucell.tpiba, drhor);
    for (int i = 0; i < PARAM.input.nspin*pw_basis.nrxx; ++i)
    {
        EXPECT_EQ(drhor[i], drhor_ref[i]);
    }
    delete[] drhor;
    delete[] drhor_ref;

    // nspin = 2
    PARAM.input.nspin = 2;
    cfg = make_cfg();
    cfg.mixing_gg0 = 0.0;
    std::complex<double>* drhog = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    std::complex<double>* drhog_old = new std::complex<double>[PARAM.input.nspin*pw_basis.npw];
    drhor = new double[PARAM.input.nspin*pw_basis.nrxx];
    drhor_ref = new double[PARAM.input.nspin*pw_basis.nrxx];
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    module_charge::kerker_screen_recip(cfg, &pw_basis, ucell.tpiba, drhog); // no kerker
    for (int i = 0; i < PARAM.input.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }

    // RECIPROCAL
    cfg.mixing_gg0 = 1.0;
    cfg.mixing_gg0_mag = 0.0;
    module_charge::kerker_screen_recip(cfg, &pw_basis, ucell.tpiba, drhog);
    const double gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        std::complex<double> ration = drhog[i] / drhog[i+pw_basis.npw];
        double gg = this->pw_basis.gg[i];
        double ration_ref = std::max(gg / (gg + gg0), 0.1 / cfg.mixing_beta);
        EXPECT_NEAR(ration.real(), ration_ref, 1e-10);
        EXPECT_NEAR(ration.imag(), 0, 1e-10);
    }

    // REAL
    pw_basis.recip2real(drhog, drhor_ref);
    pw_basis.recip2real(drhog_old, drhor);

    cfg.mixing_gg0 = 0.0;
    cfg.mixing_gg0_mag = 0.0;
    // nothing happens
    module_charge::kerker_screen_real(cfg, &pw_basis, ucell.tpiba, drhor);

    cfg.mixing_gg0 = 1.0;
    module_charge::kerker_screen_real(cfg, &pw_basis, ucell.tpiba, drhor);
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
     PARAM.sys.double_grid = false;
    charge.set_rhopw(&pw_basis);
    const int nspin = PARAM.input.nspin = 1;
    PARAM.sys.domag_z = false;
    XC_Functional::func_type = 3;
    XC_Functional::ked_flag = true;
    PARAM.input.mixing_beta = 0.7;
    PARAM.input.mixing_ndim = 1;
    PARAM.input.mixing_gg0 = 0.0;
    PARAM.input.mixing_tau = true;
    PARAM.input.mixing_mode = "plain";
    const int nrxx = pw_basis.nrxx;
    const int npw = pw_basis.npw;
    charge._space_rho.resize(nspin * nrxx);
    charge._space_rho_save.resize(nspin * nrxx);
    charge._space_rhog.resize(nspin * npw);
    charge._space_rhog_save.resize(nspin * npw);
    charge._space_kin_r.resize(nspin * nrxx);
    charge._space_kin_r_save.resize(nspin * nrxx);
    charge.rho = new double*[nspin];
    charge.rhog = new std::complex<double>*[nspin];
    charge.rho_save = new double*[nspin];
    charge.rhog_save = new std::complex<double>*[nspin];
    charge.kin_r = new double*[nspin];
    charge.kin_r_save = new double*[nspin];
    for (int is = 0; is < nspin; is++)
    {
        charge.rho[is] = charge._space_rho.data() + is * nrxx;
        charge.rhog[is] = charge._space_rhog.data() + is * npw;
        charge.rho_save[is] = charge._space_rho_save.data() + is * nrxx;
        charge.rhog_save[is] = charge._space_rhog_save.data() + is * npw;
        charge.kin_r[is] = charge._space_kin_r.data() + is * nrxx;
        charge.kin_r_save[is] = charge._space_kin_r_save.data() + is * nrxx;
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
    PARAM.input.scf_thr_type= 1;
    CMtest_recip.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
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
    PARAM.input.scf_thr_type= 2;
    CMtest_real.set_rhopw(&pw_basis, &pw_basis);
    CMtest_real.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
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
    delete[] charge.rho;
    delete[] charge.rhog;
    delete[] charge.rho_save;
    delete[] charge.rhog_save;
    delete[] charge.kin_r;
    delete[] charge.kin_r_save;
}

TEST_F(ChargeMixingTest, MixDoubleGridRhoTest)
{
     PARAM.sys.double_grid = true;
    charge.set_rhopw(&pw_dbasis);
    const int nspin = PARAM.input.nspin = 1;
    PARAM.sys.domag_z = false;
    XC_Functional::func_type = 3;
    XC_Functional::ked_flag = true;
    PARAM.input.mixing_beta = 0.7;
    PARAM.input.mixing_ndim = 1;
    PARAM.input.mixing_gg0 = 0.0;
    PARAM.input.mixing_tau = true;
    PARAM.input.mixing_mode = "plain";
    const int nrxx = pw_dbasis.nrxx;
    const int npw = pw_dbasis.npw;
    charge._space_rho.resize(nspin * nrxx);
    charge._space_rho_save.resize(nspin * nrxx);
    charge._space_rhog.resize(nspin * npw);
    charge._space_rhog_save.resize(nspin * npw);
    charge._space_kin_r.resize(nspin * nrxx);
    charge._space_kin_r_save.resize(nspin * nrxx);
    charge.rho = new double*[nspin];
    charge.rhog = new std::complex<double>*[nspin];
    charge.rho_save = new double*[nspin];
    charge.rhog_save = new std::complex<double>*[nspin];
    charge.kin_r = new double*[nspin];
    charge.kin_r_save = new double*[nspin];
    for (int is = 0; is < nspin; is++)
    {
        charge.rho[is] = charge._space_rho.data() + is * nrxx;
        charge.rhog[is] = charge._space_rhog.data() + is * npw;
        charge.rho_save[is] = charge._space_rho_save.data() + is * nrxx;
        charge.rhog_save[is] = charge._space_rhog_save.data() + is * npw;
        charge.kin_r[is] = charge._space_kin_r.data() + is * nrxx;
        charge.kin_r_save[is] = charge._space_kin_r_save.data() + is * nrxx;
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

    PARAM.input.scf_thr_type= 1;
    CMtest_recip.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);

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
    PARAM.input.nspin = 1;
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_dbasis);
    std::vector<std::complex<double>> data(pw_dbasis.npw, 1.0);
    const int npw_smooth = pw_basis.npw;
    const int npw_dense = pw_dbasis.npw;
    const int npw_hf = npw_dense - npw_smooth;

    // split: smooth + high-frequency together reconstruct the dense data
    std::vector<std::complex<double>> datas(npw_smooth);
    std::vector<std::complex<double>> datahf(npw_hf);
    module_charge::split_dgrid(data.data(), datas, datahf,
                                1, npw_smooth, npw_dense);
    for (int i = 0; i < npw_smooth; ++i)
    {
        EXPECT_EQ(datas[i], data[i]);
    }
    for (int i = 0; i < npw_hf; ++i)
    {
        EXPECT_EQ(datahf[i], data[npw_smooth + i]);
    }

    // merge: inverse of split; output must equal input
    std::vector<std::complex<double>> dataout(npw_dense, std::complex<double>(0, 0));
    module_charge::merge_dgrid(dataout.data(), datas, datahf,
                                1, npw_smooth, npw_dense);
    for (int i = 0; i < npw_dense; ++i)
    {
        EXPECT_EQ(dataout[i], data[i]);
    }

    // No explicit cleanup call needed: vectors manage their own storage.

    // NSPIN = 2
    PARAM.input.nspin = 2;
    data.resize(npw_dense * 2, 1.0);
    dataout.assign(npw_dense * 2, std::complex<double>(0, 0));
    std::vector<std::complex<double>> datas2(npw_smooth * 2);
    std::vector<std::complex<double>> datahf2(npw_hf * 2);
    module_charge::split_dgrid(data.data(), datas2, datahf2,
                                2, npw_smooth, npw_dense);
    module_charge::merge_dgrid(dataout.data(), datas2, datahf2,
                                2, npw_smooth, npw_dense);
    for (int i = 0; i < npw_dense * 2; ++i)
    {
        EXPECT_EQ(dataout[i], data[i]);
    }
}

TEST_F(ChargeMixingTest, SCFOscillationTest)
{
    Charge_Mixing CMtest;
    int scf_nmax = 20;
    int scf_os_ndim = 3;
    double scf_os_thr = -0.05;
    bool scf_oscillate = false;
    std::vector<double> drho(scf_nmax, 0.0);
    std::vector<bool> scf_oscillate_ref(scf_nmax, false);
    drho = {6.83639633652e-05,
            4.93523029235e-05,
            3.59230097735e-05,
            2.68356403913e-05,
            2.17490806464e-05,
            2.14231642508e-05,
            1.67507494811e-05,
            1.53575889539e-05,
            1.26504511554e-05,
            1.04762016224e-05,
            8.10000162918e-06,
            7.66427917682e-06,
            6.70112820094e-06,
            5.68594436664e-06,
            4.80120233733e-06,
            4.86519757184e-06,
            4.37855804356e-06,
            4.29922703412e-06,
            4.36398486331e-06,
            4.94224615955e-06};
    scf_oscillate_ref = {false,false,false,false,false,true,false,false,false,false,
                        false,false,true,false,false,true,true,true,true,true};
    for (int i = 1; i <= scf_nmax; ++i)
    {
        scf_oscillate = CMtest.if_scf_oscillate(i,drho[i-1],scf_os_ndim,scf_os_thr);
        EXPECT_EQ(scf_oscillate, scf_oscillate_ref[i-1]);
    } 
}
