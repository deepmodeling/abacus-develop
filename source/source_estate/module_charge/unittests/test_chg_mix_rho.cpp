#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_cell/magnetism.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_mix.h"
#include "source_estate/module_charge/chg_mix_cfg.h"
#include "source_io/module_parameter/parameter.h"

#include <vector>

// charge.cpp references Magnetism; provide a lightweight stub.
Magnetism::Magnetism()
{
    this->tot_mag = 0.0;
    this->abs_mag = 0.0;
}
Magnetism::~Magnetism()
{
}

/************************************************
 *  unit test of module_charge/chg_mix_rho.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Charge_Mixing::mix_rho: dispatches to mix_rho_recip (scf_thr_type==1)
 *     or mix_rho_real (scf_thr_type==2), then copies rho->rho_save.
 *     - abort on null chr / null chr->rhopw
 *     - abort when the grid was not set via set_mixing
 *     - abort when double_grid is on but rhodpw is null
 *     - real-space plain mixing: rho = rho_save + beta * (rho_new - rho_save)
 */

namespace
{

MixingConfig make_cfg(int nspin, int scf_thr_type, bool double_grid, bool mixing_tau)
{
    MixingConfig cfg{
        "plain",       // mixing_mode
        0.7,           // mixing_beta
        1,             // mixing_ndim
        0.0,           // mixing_gg0
        mixing_tau,    // mixing_tau
        1.6,           // mixing_beta_mag
        0.0,           // mixing_gg0_mag
        0.1,           // mixing_gg0_min
        -10.0,         // mixing_angle
        false,         // mixing_dmr
        nspin,         // nspin
        scf_thr_type,  // scf_thr_type
        double_grid,   // double_grid
        false,         // gamma_only_pw
        false,         // domag
        false,         // domag_z
        100            // scf_nmax
    };
    return cfg;
}

} // namespace

class ChargeMixRhoTest : public ::testing::Test
{
  public:
    ModulePW::PW_Basis pw_basis;
    ModulePW::PW_Basis_Sup pw_dbasis;
    Charge charge;
    double omega = 1.0;
    double tpiba = 1.0;

    ChargeMixRhoTest()
    {
        pw_basis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 20);
        pw_basis.initparameters(false, 20);
        pw_basis.setuptransform();
        pw_basis.collect_local_pw();
        pw_dbasis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 40);
        pw_dbasis.initparameters(false, 40);
        pw_dbasis.setuptransform(&pw_basis);
        pw_dbasis.collect_local_pw();
    }

    /// Configure a Charge_Mixing for plain mixing on the smooth grid.
    void setup_mixing(Charge_Mixing& cm, int nspin, int scf_thr_type, bool double_grid)
    {
        MixingConfig cfg = make_cfg(nspin, scf_thr_type, double_grid, false);
        if (double_grid)
        {
            cm.set_mixing(cfg, &pw_basis, &pw_dbasis, omega, tpiba);
        }
        else
        {
            cm.set_mixing(cfg, &pw_basis, &pw_basis, omega, tpiba);
        }
        cm.init_mixing();
    }

    /// Allocate Charge buffers (rho, rho_save) for nspin without kinetic density.
    void setup_charge(int nspin)
    {
        charge.set_rhopw(&pw_basis);
        const bool kin_den = false;
        const bool meta_gga = false;
        charge.allocate(nspin, kin_den, meta_gga, 0);
    }
};

// ---------------------------------------------------------------------------
// abort paths
// ---------------------------------------------------------------------------

TEST_F(ChargeMixRhoTest, MixRhoNullChrAborts)
{
    Charge_Mixing cm;
    MixingConfig cfg = make_cfg(1, 2, false, false);
    cm.set_mixing(cfg, &pw_basis, &pw_basis, omega, tpiba);
    cm.init_mixing();
    EXPECT_DEATH(cm.mix_rho(nullptr), "");
}

TEST_F(ChargeMixRhoTest, MixRhoNullChrRhopwAborts)
{
    Charge_Mixing cm;
    setup_mixing(cm, 1, 2, false);
    Charge empty_charge;
    EXPECT_DEATH(cm.mix_rho(&empty_charge), "");
}

TEST_F(ChargeMixRhoTest, MixRhoUnsetRhopwAborts)
{
    Charge_Mixing cm;
    MixingConfig cfg = make_cfg(1, 2, false, false);
    // Pass rhopw == nullptr to set_mixing to simulate a skipped grid setup.
    // Do NOT call init_mixing() here: init_mixing already WARNING_QUITs when
    // the grid is unset, which would kill the death-test parent process
    // before EXPECT_DEATH runs. The guard under test lives in mix_rho itself
    // and only checks this->rhopw == nullptr, independent of init_mixing.
    cm.set_mixing(cfg, nullptr, nullptr, omega, tpiba);
    setup_charge(1);
    EXPECT_DEATH(cm.mix_rho(&charge), "");
}

TEST_F(ChargeMixRhoTest, MixRhoDoubleGridWithoutRhodpwAborts)
{
    Charge_Mixing cm;
    MixingConfig cfg = make_cfg(1, 2, true, false);
    // set_mixing with rhodpw == nullptr while double_grid is on
    cm.set_mixing(cfg, &pw_basis, nullptr, omega, tpiba);
    cm.init_mixing();
    setup_charge(1);
    EXPECT_DEATH(cm.mix_rho(&charge), "");
}

// ---------------------------------------------------------------------------
// real-space plain mixing: rho = rho_save + beta * (rho - rho_save)
// ---------------------------------------------------------------------------

TEST_F(ChargeMixRhoTest, MixRhoRealPlainNspin1)
{
    Charge_Mixing cm;
    const int nspin = 1;
    setup_mixing(cm, nspin, 2, false);
    setup_charge(nspin);

    const double rho_save_val = 1.0;
    const double rho_new_val = 3.0;
    const double beta = 0.7;
    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        charge.rho[0][ir] = rho_new_val;
        charge.rho_save[0][ir] = rho_save_val;
    }

    cm.mix_rho(&charge);

    // mixed rho = rho_save + beta * (rho_new - rho_save) = 1 + 0.7 * 2 = 2.4
    const double expected = rho_save_val + beta * (rho_new_val - rho_save_val);
    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        EXPECT_NEAR(charge.rho[0][ir], expected, 1e-8);
    }
    // rho_save holds the pre-mixing rho (rho_new_val)
    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        EXPECT_NEAR(charge.rho_save[0][ir], rho_new_val, 1e-8);
    }
}
