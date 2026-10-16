#include "gtest/gtest.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_basis/module_pw/pw_basis.h"
#include "source_estate/module_charge/chg_mix_cfg.h"
#include "source_estate/module_charge/chg_precond.h"

#include <algorithm>
#include <complex>
#include <vector>

/************************************************
 *  unit test of module_charge/chg_precond.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - kerker_screen_recip: multiply drhog[is*npw+ig] by
 *     max(gg/(gg+gg0), gg0_min/amin) per spin channel
 *     - early return when mixing_gg0 <= 0 or mixing_beta <= 0.1
 *     - nspin == 1: density channel only
 *     - nspin == 2: density + magnetization (mag skipped when gg0_mag/beta_mag too small)
 *     - nspin == 4: density + magnetization channels (resize_tmp with mixing_angle)
 *   - kerker_screen_real: FFT to reciprocal, apply (1 - filter_g), FFT back, subtract
 *     - early return when mixing_gg0 <= 0.0001 or mixing_beta <= 0.1
 *     - nspin == 1 real-space filtering matches reciprocal-space result
 */

namespace
{

MixingConfig make_cfg()
{
    MixingConfig cfg{
        "broyden",  // mixing_mode
        0.8,        // mixing_beta
        8,          // mixing_ndim
        1.0,        // mixing_gg0
        false,      // mixing_tau
        1.6,        // mixing_beta_mag
        0.0,        // mixing_gg0_mag
        0.1,        // mixing_gg0_min
        -10.0,      // mixing_angle
        false,      // mixing_dmr
        1,          // nspin
        2,          // scf_thr_type
        false,      // double_grid
        false,      // gamma_only_pw
        false,      // domag
        false,      // domag_z
        100         // scf_nmax
    };
    return cfg;
}

} // namespace

class ChgPrecondTest : public ::testing::Test
{
  protected:
    ModulePW::PW_Basis pw_basis;
    const double tpiba = 1.0;

    void SetUp() override
    {
        pw_basis.initgrids(4, ModuleBase::Matrix3(1, 0, 0, 0, 1, 0, 0, 0, 1), 20);
        pw_basis.initparameters(false, 20);
        pw_basis.setuptransform();
        pw_basis.collect_local_pw();
    }
};

// ---------------------------------------------------------------------------
// kerker_screen_recip
// ---------------------------------------------------------------------------

TEST_F(ChgPrecondTest, KerkerScreenRecipEarlyReturnGg0Zero)
{
    MixingConfig cfg = make_cfg();
    cfg.nspin = 1;
    cfg.mixing_gg0 = 0.0;

    std::vector<std::complex<double>> drhog(pw_basis.npw, std::complex<double>(1.0, 1.0));
    std::vector<std::complex<double>> drhog_old = drhog;

    module_charge::kerker_screen_recip(cfg, &pw_basis, tpiba, drhog.data());

    for (int ig = 0; ig < pw_basis.npw; ++ig)
    {
        EXPECT_EQ(drhog[ig], drhog_old[ig]);
    }
}

TEST_F(ChgPrecondTest, KerkerScreenRecipEarlyReturnBetaTooSmall)
{
    MixingConfig cfg = make_cfg();
    cfg.nspin = 1;
    cfg.mixing_gg0 = 1.0;
    cfg.mixing_beta = 0.1; // <= 0.1 triggers early return

    std::vector<std::complex<double>> drhog(pw_basis.npw, std::complex<double>(1.0, 1.0));
    std::vector<std::complex<double>> drhog_old = drhog;

    module_charge::kerker_screen_recip(cfg, &pw_basis, tpiba, drhog.data());

    for (int ig = 0; ig < pw_basis.npw; ++ig)
    {
        EXPECT_EQ(drhog[ig], drhog_old[ig]);
    }
}

TEST_F(ChgPrecondTest, KerkerScreenRecipNspin1Filter)
{
    MixingConfig cfg = make_cfg();
    cfg.nspin = 1;
    cfg.mixing_gg0 = 1.0;
    cfg.mixing_beta = 0.8;
    cfg.mixing_gg0_min = 0.1;

    std::vector<std::complex<double>> drhog(pw_basis.npw, std::complex<double>(1.0, 1.0));

    module_charge::kerker_screen_recip(cfg, &pw_basis, tpiba, drhog.data());

    const double gg0 = std::pow(cfg.mixing_gg0 * ModuleBase::BOHR_TO_A / tpiba, 2);
    const double gg0_amin = cfg.mixing_gg0_min / cfg.mixing_beta;
    for (int ig = 0; ig < pw_basis.npw; ++ig)
    {
        const double gg = pw_basis.gg[ig];
        const double ref = std::max(gg / (gg + gg0), gg0_amin);
        EXPECT_NEAR(drhog[ig].real(), ref, 1e-10);
        EXPECT_NEAR(drhog[ig].imag(), ref, 1e-10);
    }
}

TEST_F(ChgPrecondTest, KerkerScreenRecipNspin2MagSkippedWhenGg0MagZero)
{
    MixingConfig cfg = make_cfg();
    cfg.nspin = 2;
    cfg.mixing_gg0 = 1.0;
    cfg.mixing_beta = 0.8;
    cfg.mixing_gg0_mag = 0.0; // magnetization channel is skipped (break)
    cfg.mixing_gg0_min = 0.1;

    std::vector<std::complex<double>> drhog(2 * pw_basis.npw, std::complex<double>(1.0, 1.0));

    module_charge::kerker_screen_recip(cfg, &pw_basis, tpiba, drhog.data());

    const double gg0 = std::pow(cfg.mixing_gg0 * ModuleBase::BOHR_TO_A / tpiba, 2);
    const double gg0_amin = cfg.mixing_gg0_min / cfg.mixing_beta;
    for (int ig = 0; ig < pw_basis.npw; ++ig)
    {
        const double gg = pw_basis.gg[ig];
        const double ref = std::max(gg / (gg + gg0), gg0_amin);
        // density channel is filtered
        EXPECT_NEAR(drhog[ig].real(), ref, 1e-10);
        EXPECT_NEAR(drhog[ig].imag(), ref, 1e-10);
        // magnetization channel is untouched (break before processing is=1)
        EXPECT_NEAR(drhog[pw_basis.npw + ig].real(), 1.0, 1e-10);
        EXPECT_NEAR(drhog[pw_basis.npw + ig].imag(), 1.0, 1e-10);
    }
}

TEST_F(ChgPrecondTest, KerkerScreenRecipNspin2MagChannelFiltered)
{
    MixingConfig cfg = make_cfg();
    cfg.nspin = 2;
    cfg.mixing_gg0 = 1.0;
    cfg.mixing_beta = 0.8;
    cfg.mixing_gg0_mag = 2.0;
    cfg.mixing_beta_mag = 1.6;
    cfg.mixing_gg0_min = 0.1;

    std::vector<std::complex<double>> drhog(2 * pw_basis.npw, std::complex<double>(1.0, 1.0));

    module_charge::kerker_screen_recip(cfg, &pw_basis, tpiba, drhog.data());

    const double gg0_rho = std::pow(cfg.mixing_gg0 * ModuleBase::BOHR_TO_A / tpiba, 2);
    const double gg0_mag = std::pow(cfg.mixing_gg0_mag * ModuleBase::BOHR_TO_A / tpiba, 2);
    const double gg0_amin_rho = cfg.mixing_gg0_min / cfg.mixing_beta;
    const double gg0_amin_mag = cfg.mixing_gg0_min / cfg.mixing_beta_mag;
    for (int ig = 0; ig < pw_basis.npw; ++ig)
    {
        const double gg = pw_basis.gg[ig];
        const double ref_rho = std::max(gg / (gg + gg0_rho), gg0_amin_rho);
        const double ref_mag = std::max(gg / (gg + gg0_mag), gg0_amin_mag);
        EXPECT_NEAR(drhog[ig].real(), ref_rho, 1e-10);
        EXPECT_NEAR(drhog[pw_basis.npw + ig].real(), ref_mag, 1e-10);
    }
}

TEST_F(ChgPrecondTest, KerkerScreenRecipNspin4WithAngle)
{
    MixingConfig cfg = make_cfg();
    cfg.nspin = 4;
    cfg.mixing_gg0 = 1.0;
    cfg.mixing_beta = 0.8;
    cfg.mixing_gg0_mag = 2.0;
    cfg.mixing_beta_mag = 1.6;
    cfg.mixing_gg0_min = 0.1;
    cfg.mixing_angle = 1.0; // > 0 => resize_tmp = 2

    std::vector<std::complex<double>> drhog(4 * pw_basis.npw, std::complex<double>(1.0, 1.0));

    module_charge::kerker_screen_recip(cfg, &pw_basis, tpiba, drhog.data());

    // resize_tmp == 2 means only 4/2 = 2 channels are processed:
    // is=0 (density) and is=1 (magnetization).
    const double gg0_rho = std::pow(cfg.mixing_gg0 * ModuleBase::BOHR_TO_A / tpiba, 2);
    const double gg0_mag = std::pow(cfg.mixing_gg0_mag * ModuleBase::BOHR_TO_A / tpiba, 2);
    const double gg0_amin_rho = cfg.mixing_gg0_min / cfg.mixing_beta;
    const double gg0_amin_mag = cfg.mixing_gg0_min / cfg.mixing_beta_mag;
    for (int ig = 0; ig < pw_basis.npw; ++ig)
    {
        const double gg = pw_basis.gg[ig];
        const double ref_rho = std::max(gg / (gg + gg0_rho), gg0_amin_rho);
        const double ref_mag = std::max(gg / (gg + gg0_mag), gg0_amin_mag);
        // is=0 density channel
        EXPECT_NEAR(drhog[ig].real(), ref_rho, 1e-10);
        // is=1 magnetization channel
        EXPECT_NEAR(drhog[pw_basis.npw + ig].real(), ref_mag, 1e-10);
    }
}

// ---------------------------------------------------------------------------
// kerker_screen_real
// ---------------------------------------------------------------------------

TEST_F(ChgPrecondTest, KerkerScreenRealEarlyReturnGg0Zero)
{
    MixingConfig cfg = make_cfg();
    cfg.nspin = 1;
    cfg.mixing_gg0 = 0.0;

    std::vector<double> drhor(pw_basis.nrxx, 1.0);
    std::vector<double> drhor_old = drhor;

    module_charge::kerker_screen_real(cfg, &pw_basis, tpiba, drhor.data());

    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        EXPECT_EQ(drhor[ir], drhor_old[ir]);
    }
}

TEST_F(ChgPrecondTest, KerkerScreenRealNspin1MatchesRecip)
{
    MixingConfig cfg = make_cfg();
    cfg.nspin = 1;
    cfg.mixing_gg0 = 1.0;
    cfg.mixing_beta = 0.8;
    cfg.mixing_gg0_min = 0.1;

    // start from a uniform real-space field; its reciprocal image is nonzero
    // only at G=0, which lets us verify the (1 - filter_g) subtraction.
    std::vector<double> drhor(pw_basis.nrxx, 1.0);

    // reference: apply kerker_screen_recip to the FFT of drhor, then FFT back
    std::vector<std::complex<double>> drhog(pw_basis.npw);
    pw_basis.real2recip(drhor.data(), drhog.data());
    module_charge::kerker_screen_recip(cfg, &pw_basis, tpiba, drhog.data());
    std::vector<double> drhor_ref(pw_basis.nrxx);
    pw_basis.recip2real(drhog.data(), drhor_ref.data());

    module_charge::kerker_screen_real(cfg, &pw_basis, tpiba, drhor.data());

    for (int ir = 0; ir < pw_basis.nrxx; ++ir)
    {
        EXPECT_NEAR(drhor[ir], drhor_ref[ir], 1e-8);
    }
}
