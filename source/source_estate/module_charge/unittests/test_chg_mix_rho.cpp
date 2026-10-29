#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis_sup.h"
#include "source_cell/magnetism.h"
#include "source_estate/module_charge/charge.h"
#include "source_estate/module_charge/chg_mix.h"
#include "source_estate/module_charge/chg_mix_cfg.h"
#include "source_io/module_parameter/parameter.h"

#include <algorithm>
#include <cmath>
#include <utility>
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


namespace
{
void initialize_dense_charge(Charge& charge, const int nspin, const int nrxx, const bool kinetic)
{
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            const double wave = std::sin(0.11 * ir + 0.7 * is);
            charge.rho_save[is][ir] = 1.0 + 0.2 * is + 0.05 * wave;
            charge.rho[is][ir] = 1.1 + 0.3 * is + 0.07 * wave;
            if (kinetic)
            {
                charge.kin_r_save[is][ir] = 0.5 + 0.1 * is + 0.03 * wave;
                charge.kin_r[is][ir] = 0.8 + 0.15 * is + 0.06 * wave;
            }
        }
    }
}
}

class ChargeMixDensePlainTest : public ChargeMixRhoTest,
                                public ::testing::WithParamInterface<int>
{
};

TEST_P(ChargeMixDensePlainTest, MixesAllDensePointsAndSpinChannels)
{
    const int nspin = GetParam();
    const int nrxx = pw_dbasis.nrxx;
    ASSERT_GT(nrxx, pw_basis.nrxx);
    MixingConfig cfg = make_cfg(nspin, 2, true, true);
    cfg.mixing_beta_mag = 0.3;
    Charge_Mixing cm;
    cm.set_mixing(cfg, &pw_basis, &pw_dbasis, omega, tpiba);
    cm.init_mixing();
    charge.set_rhopw(&pw_dbasis);
    charge.allocate(nspin, true, true, 0);
    initialize_dense_charge(charge, nspin, nrxx, true);
    std::vector<double> original(nspin * nrxx);
    std::vector<double> original_tau(nspin * nrxx);
    std::vector<double> expected(nspin * nrxx);
    std::vector<double> expected_tau(nspin * nrxx);
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            const int index = is * nrxx + ir;
            original[index] = charge.rho[is][ir];
            original_tau[index] = charge.kin_r[is][ir];
            double beta = cfg.mixing_beta;
            if (nspin == 4 && is > 0)
            {
                beta = cfg.mixing_beta_mag;
            }
            expected[index] = charge.rho_save[is][ir]
                            + beta * (charge.rho[is][ir] - charge.rho_save[is][ir]);
            expected_tau[index] = charge.kin_r_save[is][ir]
                                + cfg.mixing_beta * (charge.kin_r[is][ir] - charge.kin_r_save[is][ir]);
        }
    }
    if (nspin == 2)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            const double sum_in = charge.rho_save[0][ir] + charge.rho_save[1][ir];
            const double sum_out = charge.rho[0][ir] + charge.rho[1][ir];
            const double mag_in = charge.rho_save[0][ir] - charge.rho_save[1][ir];
            const double mag_out = charge.rho[0][ir] - charge.rho[1][ir];
            const double sum = sum_in + cfg.mixing_beta * (sum_out - sum_in);
            const double mag = mag_in + cfg.mixing_beta_mag * (mag_out - mag_in);
            expected[ir] = 0.5 * (sum + mag);
            expected[nrxx + ir] = 0.5 * (sum - mag);
        }
    }
    cm.mix_rho(&charge);
    double max_error = 0.0;
    bool all_finite = true;
    for (int is = 0; is < nspin; ++is)
    {
        for (int ir = 0; ir < nrxx; ++ir)
        {
            const int index = is * nrxx + ir;
            const double rho_error = std::abs(charge.rho[is][ir] - expected[index]);
            const double tau_error = std::abs(charge.kin_r[is][ir] - expected_tau[index]);
            const double rho_save_error = std::abs(charge.rho_save[is][ir] - original[index]);
            const double tau_save_error = std::abs(charge.kin_r_save[is][ir] - original_tau[index]);
            all_finite = all_finite && std::isfinite(rho_error) && std::isfinite(tau_error);
            max_error = std::max(max_error, rho_error);
            max_error = std::max(max_error, tau_error);
            all_finite = all_finite && std::isfinite(rho_save_error) && std::isfinite(tau_save_error);
            max_error = std::max(max_error, rho_save_error);
            max_error = std::max(max_error, tau_save_error);
        }
    }
    EXPECT_TRUE(all_finite);
    EXPECT_LT(max_error, 1.0e-12);
}

INSTANTIATE_TEST_SUITE_P(SpinChannels, ChargeMixDensePlainTest, ::testing::Values(1, 2, 4));

class ChargeMixDenseHistoryTest : public ChargeMixRhoTest,
                                  public ::testing::WithParamInterface<std::pair<int, bool>>
{
};

TEST_P(ChargeMixDenseHistoryTest, DoubleGridRealKerkerAndHistoryMatchDenseGridOnly)
{
    const int nspin = GetParam().first;
    const int nrxx = pw_dbasis.nrxx;
    MixingConfig cfg = make_cfg(nspin, 2, true, true);
    cfg.mixing_mode = "broyden";
    cfg.mixing_ndim = 3;
    cfg.mixing_beta_mag = 0.3;
    cfg.mixing_gg0 = 1.0;
    cfg.mixing_gg0_mag = 0.8;
    if (GetParam().second)
    {
        cfg.mixing_angle = 0.5;
    }
    Charge_Mixing actual;
    actual.set_mixing(cfg, &pw_basis, &pw_dbasis, omega, tpiba);
    actual.init_mixing();
    cfg.double_grid = false;
    Charge_Mixing control;
    control.set_mixing(cfg, &pw_dbasis, &pw_dbasis, omega, tpiba);
    control.init_mixing();
    charge.set_rhopw(&pw_dbasis);
    charge.allocate(nspin, true, true, 0);
    Charge dense_control;
    dense_control.set_rhopw(&pw_dbasis);
    dense_control.allocate(nspin, true, true, 0);
    initialize_dense_charge(charge, nspin, nrxx, true);
    initialize_dense_charge(dense_control, nspin, nrxx, true);
    double max_error = 0.0;
    bool all_finite = true;
    for (int step = 0; step < 3; ++step)
    {
        if (step > 0)
        {
            for (int is = 0; is < nspin; ++is)
            {
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    const double perturbation = 0.03 * std::cos(0.07 * ir + 0.4 * is + 0.9 * step);
                    charge.rho[is][ir] += perturbation;
                    dense_control.rho[is][ir] += perturbation;
                    charge.kin_r[is][ir] += perturbation;
                    dense_control.kin_r[is][ir] += perturbation;
                }
            }
        }
        actual.mix_rho(&charge);
        control.mix_rho(&dense_control);
        for (int is = 0; is < nspin; ++is)
        {
            for (int ir = 0; ir < nrxx; ++ir)
            {
                const double rho_error = std::abs(charge.rho[is][ir] - dense_control.rho[is][ir]);
                const double tau_error = std::abs(charge.kin_r[is][ir] - dense_control.kin_r[is][ir]);
                all_finite = all_finite && std::isfinite(rho_error) && std::isfinite(tau_error);
                max_error = std::max(max_error, rho_error);
                max_error = std::max(max_error, tau_error);
            }
        }
    }
    EXPECT_TRUE(all_finite);
    EXPECT_LT(max_error, 1.0e-12);
}


namespace
{
const std::pair<int, bool> collinear_history(2, false);
const std::pair<int, bool> spinor_history(4, false);
const std::pair<int, bool> angle_history(4, true);
}
INSTANTIATE_TEST_SUITE_P(SpinChannels,
                         ChargeMixDenseHistoryTest,
                         ::testing::Values(collinear_history, spinor_history, angle_history));

TEST_F(ChargeMixRhoTest, RealDoubleGridInitWithoutDenseGridAborts)
{
    Charge_Mixing cm;
    MixingConfig cfg = make_cfg(1, 2, true, false);
    cm.set_mixing(cfg, &pw_basis, nullptr, omega, tpiba);
    EXPECT_DEATH(cm.init_mixing(), "");
}
