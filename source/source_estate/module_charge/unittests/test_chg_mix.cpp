#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "../chg_mix.h"
#include "../chg_drho.h"
#include "../chg_drho_detail.h"
#include "../chg_precond.h"
#include "../chg_uspp.h"
#include "source_base/module_mixing/broyden_mixing.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_hamilt/module_xc/xc_functional.h"

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
 *                    Charge_Mixing::get_mixing_config()
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
 *   - CloseKerkerGg0DisablesScreenReal: Charge_Mixing::close_kerker_gg0()
 *      - regression test: close_kerker_gg0() must short-circuit the Kerker
 *        screening lambda in mix_rho_real so output matches cfg.mixing_gg0=0
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
        cfg.mixing_mode = "broyden";
        cfg.mixing_beta = 0.8;
        cfg.mixing_ndim = 8;
        cfg.mixing_gg0  = 1.0;
        cfg.mixing_tau  = false;
        cfg.mixing_beta_mag = 1.6;
        cfg.mixing_gg0_mag = 0.0;
        cfg.mixing_gg0_min = 0.1;
        cfg.mixing_angle = -10.0;
        cfg.mixing_dmr = false;
        // the remaining fields are runtime globals the tests steer per case;
        // seed them with the Input_para / System_para defaults
        cfg.nspin = 1;
        cfg.scf_thr_type = -1;
        cfg.double_grid = false;
        cfg.gamma_only_pw = false;
        cfg.domag = false;
        cfg.domag_z = false;
        cfg.scf_nmax = 100;
        ucell.omega = 1.0;
        ucell.tpiba = 1.0;
    }
    ModulePW::PW_Basis pw_basis;
    ModulePW::PW_Basis_Sup pw_dbasis;
    Charge charge;

    /// The mixing configuration under test. The production mixing code takes
    /// this snapshot explicitly rather than reading any global, so the fixture
    /// owns it outright instead of staging values in the parameter singleton;
    /// each case amends the fields it cares about and hands it to set_mixing().
    MixingConfig cfg;

    /// Backing storage for the Charge views under test. charge.cpp is not one
    /// of this target's sources, so Charge::allocate() is unavailable here and
    /// the fixture supplies the buffers itself, pointing the public
    /// rho/rhog/kin_r views at them with the same is*nrxx (resp. is*npw)
    /// stride allocate() would lay out. Owning them as vectors also means the
    /// fixture cleans up after itself.
    std::vector<double> buf_rho, buf_rho_save, buf_kin_r, buf_kin_r_save;
    std::vector<std::complex<double>> buf_rhog, buf_rhog_save;
    std::vector<double*> ptr_rho, ptr_rho_save, ptr_kin_r, ptr_kin_r_save;
    std::vector<std::complex<double>*> ptr_rhog, ptr_rhog_save;

    void wire_charge(const int nspin, const int nrxx, const int npw, const bool kin_den)
    {
        buf_rho.assign(nspin * nrxx, 0.0);
        buf_rho_save.assign(nspin * nrxx, 0.0);
        buf_rhog.assign(nspin * npw, std::complex<double>(0.0, 0.0));
        buf_rhog_save.assign(nspin * npw, std::complex<double>(0.0, 0.0));
        ptr_rho.resize(nspin);
        ptr_rho_save.resize(nspin);
        ptr_rhog.resize(nspin);
        ptr_rhog_save.resize(nspin);
        for (int is = 0; is < nspin; is++)
        {
            ptr_rho[is] = buf_rho.data() + is * nrxx;
            ptr_rho_save[is] = buf_rho_save.data() + is * nrxx;
            ptr_rhog[is] = buf_rhog.data() + is * npw;
            ptr_rhog_save[is] = buf_rhog_save.data() + is * npw;
        }
        charge.rho = ptr_rho.data();
        charge.rho_save = ptr_rho_save.data();
        charge.rhog = ptr_rhog.data();
        charge.rhog_save = ptr_rhog_save.data();
        if (kin_den)
        {
            buf_kin_r.assign(nspin * nrxx, 0.0);
            buf_kin_r_save.assign(nspin * nrxx, 0.0);
            ptr_kin_r.resize(nspin);
            ptr_kin_r_save.resize(nspin);
            for (int is = 0; is < nspin; is++)
            {
                ptr_kin_r[is] = buf_kin_r.data() + is * nrxx;
                ptr_kin_r_save[is] = buf_kin_r_save.data() + is * nrxx;
            }
            charge.kin_r = ptr_kin_r.data();
            charge.kin_r_save = ptr_kin_r_save.data();
        }
    }

    // mixing_tau is only honoured when the functional actually carries a
    // kinetic energy density, which set_mixing's callers fold in. Cases flip
    // XC_Functional::set_ked_flag() around, so apply it at hand-off time.
    MixingConfig make_cfg()
    {
        MixingConfig out = cfg;
        out.mixing_tau = cfg.mixing_tau && XC_Functional::get_ked_flag();
        return out;
    }

    // Re-sync the runtime globals (nspin/scf_thr_type/gamma_only_pw/domag/domag_z)
    // into an already-configured Charge_Mixing. Cases mutate them after
    // set_mixing to steer the residual/inner-product branches, and the object
    // reads them from its own snapshot, so the test must push the new values in
    // without disturbing the mixing history that init_mixing() built.
    void sync_cfg(Charge_Mixing& cm)
    {
        MixingConfig cur = cm.get_mixing_config();
        cur.nspin = cfg.nspin;
        cur.scf_thr_type = cfg.scf_thr_type;
        cur.gamma_only_pw = cfg.gamma_only_pw;
        cur.domag = cfg.domag;
        cur.domag_z = cfg.domag_z;
        cm.set_mixing_config(cur);
    }
};

TEST_F(ChargeMixingTest, SetMixingTest)
{
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
    cfg.nspin = 1;
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    cfg.mixing_beta = 1.0;
    cfg.mixing_ndim = 1;
    cfg.mixing_gg0 = 1.0;

    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    EXPECT_EQ(CMtest.get_mixing_mode(), "broyden");
    EXPECT_EQ(CMtest.get_mixing_beta(), 1.0);
    EXPECT_EQ(CMtest.get_mixing_ndim(), 1);
    EXPECT_EQ(CMtest.get_mixing_config().mixing_gg0, 1.0);
    EXPECT_EQ(CMtest.get_mixing_config().mixing_tau, false);
    EXPECT_EQ(CMtest.get_mixing_config().mixing_beta_mag, 1.6);
    EXPECT_EQ(CMtest.get_mixing_config().mixing_gg0_mag, 0.0);
    EXPECT_EQ(CMtest.get_mixing_config().mixing_gg0_min, 0.1);
    EXPECT_EQ(CMtest.get_mixing_config().mixing_angle, -10.0);
    EXPECT_EQ(CMtest.get_mixing_config().mixing_dmr, false);

    cfg.mixing_tau = true;
    XC_Functional::set_ked_flag(true);
    cfg.mixing_mode = "plain";
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    EXPECT_EQ(CMtest.get_mixing_mode(), "plain");
    EXPECT_EQ(CMtest.get_mixing_config().mixing_tau, true);
    XC_Functional::set_ked_flag(false);

    cfg.mixing_beta = 1.1;
    std::string output;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);, ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("You'd better set mixing_beta to [0.0, 1.0]!"));

    cfg.mixing_beta = 0.7;
    cfg.mixing_beta_mag = -0.1;
    cfg.nspin = 2;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);, ::testing::ExitedWithCode(1), "");
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("You'd better set mixing_beta_mag >= 0.0!"));

    cfg.nspin = 1;
    cfg.mixing_beta = 0.7;
    cfg.mixing_beta_mag = 1.6;
    cfg.mixing_mode = "nothing";
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
    cfg.nspin = 1;
    XC_Functional::set_func_type(1);
    XC_Functional::set_ked_flag(false);
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);

    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    
    cfg.scf_thr_type= 1;
    sync_cfg(CMtest);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.get_rho_mdata().length, pw_basis.npw);
    
    cfg.scf_thr_type= 2;
    sync_cfg(CMtest);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.get_rho_mdata().length, pw_basis.nrxx);

    cfg.nspin = 4;
    sync_cfg(CMtest);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.get_rho_mdata().length, 4 * pw_basis.nrxx);

    cfg.nspin = 1;
    cfg.mixing_tau = true;
    XC_Functional::set_func_type(3);
    XC_Functional::set_ked_flag(true);
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.get_tau_mdata().length, pw_basis.nrxx);

    cfg.nspin = 4;
    cfg.mixing_angle = 1.0;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    CMtest.init_mixing();
    EXPECT_EQ(CMtest.get_rho_mdata().length, 2 * pw_basis.nrxx);
}

TEST_F(ChargeMixingTest, InnerDotRealTest)
{
    Charge_Mixing CMtest;
    // non mixing angle case
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    cfg.nspin = 4;
    sync_cfg(CMtest);

    // a simple sum for inner product
    std::vector<double> drho1(pw_basis.nrxx * cfg.nspin);
    std::vector<double> drho2(pw_basis.nrxx * cfg.nspin);
    for (int i = 0; i < pw_basis.nrxx * cfg.nspin; ++i)
    {
        drho1[i] = 1.0;
        drho2[i] = double(i);
    }
    double inner = module_charge::inner_product_real(drho1.data(), drho2.data(), pw_basis, CMtest.get_mixing_config());
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * cfg.nspin  * (pw_basis.nrxx * cfg.nspin - 1), 1e-8);

    // mixing angle case
    cfg.mixing_angle = 1.0;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    cfg.nspin = 4;

    // a simple sum for inner product
    drho1.resize(pw_basis.nrxx * 2);
    drho2.resize(pw_basis.nrxx * 2);
    for (int i = 0; i < pw_basis.nrxx * 2; ++i)
    {
        drho1[i] = 1.0;
        drho2[i] = double(i);
    }
    inner = module_charge::inner_product_real(drho1.data(), drho2.data(), pw_basis, CMtest.get_mixing_config());
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * 2  * (pw_basis.nrxx * 2 - 1), 1e-8);
}

TEST_F(ChargeMixingTest, InnerDotRecipHartreeTest)
{
    // REAL
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    const int npw = pw_basis.npw;
    const int nrxx = pw_basis.nrxx;
    cfg.nspin = 1;
    std::vector<double> drhor1(pw_basis.nrxx);
    std::vector<double> drhor2(pw_basis.nrxx);
    for (int i = 0; i < pw_basis.nrxx; ++i)
    {
        drhor1[i] = 1.0;
        drhor2[i] = double(i);
    }
    // Populate cfg_ before the first inner_product call: the function reads
    // nspin from cfg_, which is default-constructed (and thus invalid) until
    // set_mixing runs.
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    double inner = module_charge::inner_product_real(drhor1.data(), drhor2.data(), pw_basis, CMtest.get_mixing_config());
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * (pw_basis.nrxx - 1), 1e-8);

    // RECIPROCAL NSPIN=1
    ucell.tpiba2 = 1.0;
    ucell.omega = 2.0;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    cfg.nspin = 1;
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

    inner = module_charge::inner_product_recip_hartree(drhog1.data(), drhog2.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, -0.3 * ModuleBase::e2 * ModuleBase::FOUR_PI, 1e-8);

    // RECIPROCAL NSPIN=2
    cfg.nspin = 2;
    sync_cfg(CMtest);
    drhog1.resize(pw_basis.npw * cfg.nspin);
    drhog2.resize(pw_basis.npw * cfg.nspin);
    std::vector<std::complex<double>> drhog1_mag(pw_basis.npw * cfg.nspin);
    std::vector<std::complex<double>> drhog2_mag(pw_basis.npw * cfg.nspin);
    for (int i = 0; i < pw_basis.npw * cfg.nspin; ++i)
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
    cfg.gamma_only_pw= false;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1_mag.data(), drhog2_mag.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 236763.82650318215, 1e-8);
    cfg.gamma_only_pw= true;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1_mag.data(), drhog2_mag.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 236763.82650318215 * 2, 1e-8);

    // RECIPROCAL NSPIN=4 without mixing_angle
    cfg.nspin = 4;
    sync_cfg(CMtest);
    drhog1.resize(pw_basis.npw * cfg.nspin);
    drhog2.resize(pw_basis.npw * cfg.nspin);
    for (int i = 0; i < pw_basis.npw * cfg.nspin; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }

    cfg.domag = false;
    cfg.domag_z = false;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1.data(), drhog2.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 28260.091995611871, 1e-8);
    cfg.gamma_only_pw= true;
    cfg.domag = true;
    cfg.domag_z = true;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1.data(), drhog2.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 110668.61166927818, 1e-8);

    // RECIPROCAL NSPIN=4 with mixing_angle
    cfg.nspin = 4;
    cfg.mixing_angle = 1.0;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    drhog1.resize(pw_basis.npw * 2);
    drhog2.resize(pw_basis.npw * 2);
    for (int i = 0; i < pw_basis.npw * 2; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }
    cfg.gamma_only_pw= false;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1.data(), drhog2.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 36548.881431837777, 1e-8);
    cfg.gamma_only_pw= true;
    sync_cfg(CMtest);
    inner = module_charge::inner_product_recip_hartree(drhog1.data(), drhog2.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 44776.555369916401, 1e-8);
}

TEST_F(ChargeMixingTest, InnerDotRecipRhoTest)
{
    // REAL
    Charge_Mixing CMtest;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    cfg.nspin = 1;
    std::vector<double> drhor1(pw_basis.nrxx);
    std::vector<double> drhor2(pw_basis.nrxx);
    for (int i = 0; i < pw_basis.nrxx; ++i)
    {
        drhor1[i] = 1.0;
        drhor2[i] = double(i);
    }
    // Populate cfg_ before the first inner_product call (see the hartree test).
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    double inner = module_charge::inner_product_real(drhor1.data(), drhor2.data(), pw_basis, CMtest.get_mixing_config());
    EXPECT_NEAR(inner, 0.5 * pw_basis.nrxx * (pw_basis.nrxx - 1), 1e-8);

    // RECIPROCAL
    ucell.tpiba2 = 1.0;
    ucell.omega = 2.0;
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    cfg.nspin = 1;
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

    inner = module_charge::detail::inner_product_recip_rho(drhog1.data(), drhog2.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, -0.3 * ModuleBase::e2 * ModuleBase::FOUR_PI, 1e-8);

    cfg.nspin = 2;
    sync_cfg(CMtest);
    drhog1.resize(pw_basis.npw * cfg.nspin);
    drhog2.resize(pw_basis.npw * cfg.nspin);
    for (int i = 0; i < pw_basis.npw * cfg.nspin; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }
    cfg.gamma_only_pw= false;
    sync_cfg(CMtest);
    inner = module_charge::detail::inner_product_recip_rho(drhog1.data(), drhog2.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 236763.82650318215, 1e-8);
    cfg.gamma_only_pw= true;
    sync_cfg(CMtest);
    inner = module_charge::detail::inner_product_recip_rho(drhog1.data(), drhog2.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 236763.82650318215 * 2, 1e-8);

    cfg.nspin = 4;
    sync_cfg(CMtest);
    drhog1.resize(pw_basis.npw * cfg.nspin);
    drhog2.resize(pw_basis.npw * cfg.nspin);
    for (int i = 0; i < pw_basis.npw * cfg.nspin; ++i)
    {
        drhog1[i] = std::complex<double>(1.0, double(i));
        drhog2[i] = std::complex<double>(1.0, 1.0);
    }

    cfg.domag = false;
    cfg.domag_z = false;
    sync_cfg(CMtest);
    inner = module_charge::detail::inner_product_recip_rho(drhog1.data(), drhog2.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 28260.091995611871, 1e-8);
    cfg.gamma_only_pw= true;
    cfg.domag = true;
    cfg.domag_z = true;
    sync_cfg(CMtest);
    inner = module_charge::detail::inner_product_recip_rho(drhog1.data(), drhog2.data(), pw_basis, CMtest.get_mixing_config(), ucell.omega, ucell.tpiba);
    EXPECT_NEAR(inner, 110668.61166927818, 1e-8);
}

TEST_F(ChargeMixingTest, KerkerScreenRecipTest)
{
    ucell.tpiba = 1.0;
    // nspin = 1
    cfg.nspin = 1;
    // kcfg is the snapshot handed to the kernels; cfg (the fixture member) is
    // the source make_cfg() reads, so nspin has to be set on cfg, not on kcfg.
    MixingConfig kcfg = make_cfg();
    std::complex<double>* drhog = new std::complex<double>[cfg.nspin*pw_basis.npw];
    std::complex<double>* drhog_old = new std::complex<double>[cfg.nspin*pw_basis.npw];
    for (int i = 0; i < cfg.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    // no kerker
    kcfg.mixing_gg0 = 0.0;
    module_charge::kerker_screen_recip(kcfg, &pw_basis, ucell.tpiba, drhog);
    for (int i = 0; i < cfg.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }
    // kerker
    kcfg.mixing_gg0 = 1.0;
    module_charge::kerker_screen_recip(kcfg, &pw_basis, ucell.tpiba, drhog);
    double gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg0), 0.1 / kcfg.mixing_beta);
        EXPECT_NEAR(drhog[i].real(), ref, 1e-10);
        EXPECT_NEAR(drhog[i].imag(), ref, 1e-10);
    }
    delete[] drhog;
    delete[] drhog_old;

    // nspin = 2
    cfg.nspin = 2;
    kcfg = make_cfg();
    kcfg.mixing_beta = 0.4;
    kcfg.mixing_beta_mag = 1.6;
    drhog = new std::complex<double>[cfg.nspin*pw_basis.npw];
    drhog_old = new std::complex<double>[cfg.nspin*pw_basis.npw];
    for (int i = 0; i < cfg.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    // mixing_gg0 = 0.0
    kcfg.mixing_gg0 = 0.0;
    module_charge::kerker_screen_recip(kcfg, &pw_basis, ucell.tpiba, drhog);
    for (int i = 0; i < cfg.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }
    // mixing_gg0 = 1.0, mixing_gg0_mag = 0.0
    kcfg.mixing_gg0 = 1.0;
    module_charge::kerker_screen_recip(kcfg, &pw_basis, ucell.tpiba, drhog);
    gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg0), 0.1 / kcfg.mixing_beta);
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
    cfg.nspin = 4;
    kcfg = make_cfg();
    drhog = new std::complex<double>[cfg.nspin*pw_basis.npw];
    drhog_old = new std::complex<double>[cfg.nspin*pw_basis.npw];
    for (int i = 0; i < cfg.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    // mixing_gg0 = 0.0
    kcfg.mixing_gg0 = 0.0;
    module_charge::kerker_screen_recip(kcfg, &pw_basis, ucell.tpiba, drhog);
    for (int i = 0; i < cfg.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }
    // mixing_gg0 = 1.0, mixing_gg0_mag = 0.0
    kcfg.mixing_gg0 = 1.0;
    module_charge::kerker_screen_recip(kcfg, &pw_basis, ucell.tpiba, drhog);
    gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg0), 0.1 / kcfg.mixing_beta);
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
    kcfg.mixing_gg0 = 1.0;
    kcfg.mixing_gg0_mag = 2.0;
    module_charge::kerker_screen_recip(kcfg, &pw_basis, ucell.tpiba, drhog);
    double gg1 = std::pow(1.0 * ModuleBase::BOHR_TO_A, 2);
    double gg2 = std::pow(2.0 * ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg1), 0.1 / kcfg.mixing_beta);
        // rho
        EXPECT_NEAR(drhog[i].real(), ref * ref, 1e-10);
        EXPECT_NEAR(drhog[i].imag(), ref * ref, 1e-10);
    }
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        double gg = this->pw_basis.gg[i];
        double ref = std::max(gg / (gg + gg2), 0.1 / kcfg.mixing_beta_mag);
        // rho
        for (int j = 1; j < cfg.nspin; ++j)
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
    cfg.nspin = 1;
    // kcfg is the snapshot handed to the kernels; nspin is set on the fixture
    // member cfg, which is what make_cfg() reads.
    MixingConfig kcfg = make_cfg();
    double* drhor = new double[cfg.nspin*pw_basis.nrxx];
    double* drhor_ref = new double[cfg.nspin*pw_basis.nrxx];
    for (int i = 0; i < cfg.nspin*pw_basis.nrxx; ++i)
    {
        drhor_ref[i] = drhor[i] = 1.0;
    }
    // no kerker
    kcfg.mixing_gg0 = 0.0;
    module_charge::kerker_screen_real(kcfg, &pw_basis, ucell.tpiba, drhor);
    for (int i = 0; i < cfg.nspin*pw_basis.nrxx; ++i)
    {
        EXPECT_EQ(drhor[i], drhor_ref[i]);
    }
    delete[] drhor;
    delete[] drhor_ref;

    // nspin = 2
    cfg.nspin = 2;
    kcfg = make_cfg();
    kcfg.mixing_gg0 = 0.0;
    std::complex<double>* drhog = new std::complex<double>[cfg.nspin*pw_basis.npw];
    std::complex<double>* drhog_old = new std::complex<double>[cfg.nspin*pw_basis.npw];
    drhor = new double[cfg.nspin*pw_basis.nrxx];
    drhor_ref = new double[cfg.nspin*pw_basis.nrxx];
    for (int i = 0; i < cfg.nspin*pw_basis.npw; ++i)
    {
        drhog_old[i] = drhog[i] = std::complex<double>(1.0, 1.0);
    }
    module_charge::kerker_screen_recip(kcfg, &pw_basis, ucell.tpiba, drhog); // no kerker
    for (int i = 0; i < cfg.nspin*pw_basis.npw; ++i)
    {
        EXPECT_EQ(drhog[i], drhog_old[i]);
    }

    // RECIPROCAL
    kcfg.mixing_gg0 = 1.0;
    kcfg.mixing_gg0_mag = 0.0;
    module_charge::kerker_screen_recip(kcfg, &pw_basis, ucell.tpiba, drhog);
    const double gg0 = std::pow(ModuleBase::BOHR_TO_A, 2);
    for (int i = 0; i < pw_basis.npw; ++i)
    {
        std::complex<double> ration = drhog[i] / drhog[i+pw_basis.npw];
        double gg = this->pw_basis.gg[i];
        double ration_ref = std::max(gg / (gg + gg0), 0.1 / kcfg.mixing_beta);
        EXPECT_NEAR(ration.real(), ration_ref, 1e-10);
        EXPECT_NEAR(ration.imag(), 0, 1e-10);
    }

    // REAL
    pw_basis.recip2real(drhog, drhor_ref);
    pw_basis.recip2real(drhog_old, drhor);

    kcfg.mixing_gg0 = 0.0;
    kcfg.mixing_gg0_mag = 0.0;
    // nothing happens
    module_charge::kerker_screen_real(kcfg, &pw_basis, ucell.tpiba, drhor);

    kcfg.mixing_gg0 = 1.0;
    module_charge::kerker_screen_real(kcfg, &pw_basis, ucell.tpiba, drhor);
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
     cfg.double_grid = false;
    charge.set_rhopw(&pw_basis);
    const int nspin = cfg.nspin = 1;
    cfg.domag_z = false;
    XC_Functional::set_func_type(3);
    XC_Functional::set_ked_flag(true);
    cfg.mixing_beta = 0.7;
    cfg.mixing_ndim = 1;
    cfg.mixing_gg0 = 0.0;
    cfg.mixing_tau = true;
    cfg.mixing_mode = "plain";
    const int nrxx = pw_basis.nrxx;
    const int npw = pw_basis.npw;
    wire_charge(nspin, nrxx, npw, XC_Functional::get_ked_flag());
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
    cfg.scf_thr_type= 1;
    CMtest_recip.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    CMtest_recip.init_mixing();
    for(int i = 0 ; i < nspin * npw; ++i)
    {
        charge.rhog[i / npw][i % npw] = recip_ref[i];
        charge.rhog_save[i / npw][i % npw] = recip_save_ref[i];
    }
    for(int i = 0 ; i < nspin * nrxx; ++i)
    {
        charge.rho[i / nrxx][i % nrxx] = real_ref[i];
        charge.rho_save[i / nrxx][i % nrxx] = real_save_ref[i];
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
    cfg.scf_thr_type= 2;
    CMtest_real.set_rhopw(&pw_basis, &pw_basis);
    CMtest_real.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    CMtest_real.init_mixing();
    for(int i = 0 ; i < nspin * nrxx; ++i)
    {
        charge.rho[i / nrxx][i % nrxx] = real_ref[i];
        charge.rho_save[i / nrxx][i % nrxx] = real_save_ref[i];
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
    // the fixture owns the buffers behind charge.rho/rhog/kin_r, so there is
    // nothing to delete here
}

// Regression test: close_kerker_gg0() must short-circuit the Kerker screening
// lambda in mix_rho_real. Before the chg_precond refactor (commit 6d127d517)
// the kernels read this->mixing_gg0; after, they read cfg_ which is an
// immutable INPUT snapshot, so writing the dead member was a no-op and the
// non-separate-loop EXX path silently failed to disable Kerker. This test
// pins the fix: output after close_kerker_gg0() must match the cfg.mixing_gg0
// = 0 baseline.
TEST_F(ChargeMixingTest, CloseKerkerGg0DisablesScreenReal)
{
    cfg.double_grid = false;
    charge.set_rhopw(&pw_basis);
    const int nspin = cfg.nspin = 1;
    cfg.domag_z = false;
    XC_Functional::set_func_type(3);
    XC_Functional::set_ked_flag(false);
    cfg.mixing_beta = 0.7;
    cfg.mixing_ndim = 1;
    cfg.mixing_gg0 = 1.0; // Kerker active by default
    cfg.mixing_tau = false;
    cfg.mixing_mode = "plain";
    cfg.scf_thr_type = 2; // real-space path

    const int nrxx = pw_basis.nrxx;
    // Only rho/rho_save are exercised on this real-space path.
    wire_charge(nspin, nrxx, pw_basis.npw, XC_Functional::get_ked_flag());
    // Non-trivial real-space residual: linear ramp so Kerker (which damps
    // long wavelengths) actually changes the output vs the no-Kerker path.
    std::vector<double> real_ref(nspin * nrxx);
    std::vector<double> real_save_ref(nspin * nrxx);
    for (int i = 0; i < nspin * nrxx; ++i)
    {
        real_ref[i] = 0.3 + 0.01 * i;
        real_save_ref[i] = 0.1 + 0.005 * i;
    }

    // --- Run A: close_kerker_gg0() then mix_rho ---
    Charge_Mixing CM_disabled;
    CM_disabled.set_rhopw(&pw_basis, &pw_basis);
    CM_disabled.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    CM_disabled.init_mixing();
    CM_disabled.close_kerker_gg0();
    for (int i = 0; i < nspin * nrxx; ++i)
    {
        charge.rho[i / nrxx][i % nrxx] = real_ref[i];
        charge.rho_save[i / nrxx][i % nrxx] = real_save_ref[i];
    }
    CM_disabled.mix_rho(&charge);
    std::vector<double> rho_A(charge.rho[0], charge.rho[0] + nspin * nrxx);

    // --- Run B: cfg.mixing_gg0 = 0 baseline, no close_kerker_gg0 ---
    Charge_Mixing CM_baseline;
    CM_baseline.set_rhopw(&pw_basis, &pw_basis);
    MixingConfig cfg_off = make_cfg();
    cfg_off.mixing_gg0 = 0.0; // Kerker off at config level
    CM_baseline.set_mixing(cfg_off, ucell.omega, ucell.tpiba);
    CM_baseline.init_mixing();
    for (int i = 0; i < nspin * nrxx; ++i)
    {
        charge.rho[i / nrxx][i % nrxx] = real_ref[i];
        charge.rho_save[i / nrxx][i % nrxx] = real_save_ref[i];
    }
    CM_baseline.mix_rho(&charge);
    std::vector<double> rho_B(charge.rho[0], charge.rho[0] + nspin * nrxx);

    // close_kerker_gg0 path must match the Kerker-off baseline.
    for (int i = 0; i < nspin * nrxx; ++i)
    {
        EXPECT_NEAR(rho_A[i], rho_B[i], 1e-10)
            << "i=" << i << ": close_kerker_gg0 did not disable Kerker";
    }

    // --- Run C: Kerker active, no close_kerker_gg0. Output must differ from A
    // to prove the disable flag was load-bearing (not that Kerker was a no-op
    // for this input to begin with). ---
    Charge_Mixing CM_active;
    CM_active.set_rhopw(&pw_basis, &pw_basis);
    CM_active.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
    CM_active.init_mixing();
    for (int i = 0; i < nspin * nrxx; ++i)
    {
        charge.rho[i / nrxx][i % nrxx] = real_ref[i];
        charge.rho_save[i / nrxx][i % nrxx] = real_save_ref[i];
    }
    CM_active.mix_rho(&charge);
    std::vector<double> rho_C(charge.rho[0], charge.rho[0] + nspin * nrxx);

    bool any_diff = false;
    for (int i = 0; i < nspin * nrxx; ++i)
    {
        if (std::abs(rho_A[i] - rho_C[i]) > 1e-8)
        {
            any_diff = true;
            break;
        }
    }
    EXPECT_TRUE(any_diff)
        << "Kerker-active output equals Kerker-disabled output, so the "
           "close_kerker_gg0 test cannot prove the flag does anything";

    // the fixture owns the buffers behind charge.rho/rho_save
}

TEST_F(ChargeMixingTest, MixDoubleGridRhoTest)
{
     cfg.double_grid = true;
    charge.set_rhopw(&pw_dbasis);
    const int nspin = cfg.nspin = 1;
    cfg.domag_z = false;
    XC_Functional::set_func_type(3);
    XC_Functional::set_ked_flag(true);
    cfg.mixing_beta = 0.7;
    cfg.mixing_ndim = 1;
    cfg.mixing_gg0 = 0.0;
    cfg.mixing_tau = true;
    cfg.mixing_mode = "plain";
    const int nrxx = pw_dbasis.nrxx;
    const int npw = pw_dbasis.npw;
    // charge is pointed at pw_dbasis here, so size the buffers from the dense grid.
    wire_charge(nspin, nrxx, npw, XC_Functional::get_ked_flag());
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

    cfg.scf_thr_type= 1;
    CMtest_recip.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);

    CMtest_recip.init_mixing();
    for (int i = 0; i < nspin * npw; ++i)
    {
        charge.rhog[i / npw][i % npw] = recip_ref[i];
        charge.rhog_save[i / npw][i % npw] = recip_save_ref[i];
    }
    for (int i = 0; i < nspin * nrxx; ++i)
    {
        charge.rho[i / nrxx][i % nrxx] = real_ref[i];
        charge.rho_save[i / nrxx][i % nrxx] = real_save_ref[i];
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
    // the fixture owns the buffers behind charge.rho/rhog/kin_r, so there is
    // nothing to delete here
}

TEST_F(ChargeMixingTest, MixDivCombTest)
{
    // NSPIN = 1
    cfg.nspin = 1;
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
    cfg.nspin = 2;
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
    // if_scf_oscillate sizes _drho_history from cfg_.scf_nmax, so cfg_ must
    // be populated before the loop; a default-constructed cfg_ leaves it 0.
    cfg.scf_nmax = scf_nmax;
    CMtest.set_rhopw(&pw_basis, &pw_basis);
    CMtest.set_mixing(make_cfg(), ucell.omega, ucell.tpiba);
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
