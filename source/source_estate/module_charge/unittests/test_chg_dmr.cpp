#include "gtest/gtest.h"

#include "source_base/module_mixing/plain_mixing.h"
#include "source_estate/module_charge/chg_dmr.h"
#include "source_estate/module_charge/chg_mix_cfg.h"

#include <vector>

/************************************************
 *  unit test of module_charge/chg_dmr.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - init_mixing_dmr: allocate the DMR mixing buffer and reset history
 *     - scf_thr_type == 2: mdata.length == nnr * dmr_nspin, counters reset
 *     - scf_thr_type == 1: abort (PW basis not supported)
 *     - abort on mixing == nullptr or nnr < 0
 *   - mix_dmr: mix the real-space density matrix
 *     - nspin == 1: out = in + beta * (out_new - in)
 *     - nspin == 2: up/down -> charge/mag channels, two betas, back
 *     - nspin == 4: treated like nspin == 1 (single channel)
 *     - nnr == 0: null buffers are allowed (empty partition)
 *     - abort on invalid inputs (null mixing, bad nspin, null buffer with nnr>0)
 */

namespace
{

/// Build a MixingConfig with all fields explicitly initialized.
MixingConfig make_cfg(int nspin, double beta, double beta_mag, int scf_thr_type)
{
    MixingConfig cfg{
        "plain",     // mixing_mode
        beta,        // mixing_beta
        4,           // mixing_ndim
        0.0,         // mixing_gg0
        false,       // mixing_tau
        beta_mag,    // mixing_beta_mag
        0.0,         // mixing_gg0_mag
        0.1,         // mixing_gg0_min
        -10.0,       // mixing_angle
        true,        // mixing_dmr
        nspin,       // nspin
        scf_thr_type,// scf_thr_type
        false,       // double_grid
        false,       // gamma_only_pw
        false,       // domag
        false,       // domag_z
        100          // scf_nmax
    };
    return cfg;
}

} // namespace

class ChgDmrTest : public ::testing::Test
{
  protected:
    Base_Mixing::Plain_Mixing mixing;
    Base_Mixing::Mixing_Data mdata;
};

// ---------------------------------------------------------------------------
// init_mixing_dmr
// ---------------------------------------------------------------------------

TEST_F(ChgDmrTest, InitMixingDmrNspin1AllocatesAndResets)
{
    const int nnr = 10;
    MixingConfig cfg = make_cfg(1, 0.5, 0.5, 2);
    module_charge::init_mixing_dmr(&mixing, mdata, nnr, cfg);

    // dmr_nspin = 1 for nspin == 1
    EXPECT_EQ(mdata.length, static_cast<size_t>(nnr * 1));
    EXPECT_EQ(mdata.ndim_use, 0);
    EXPECT_EQ(mdata.ndim_history, 0);
    EXPECT_EQ(mdata.start, -1);
}

TEST_F(ChgDmrTest, InitMixingDmrNspin2AllocatesTwoChannels)
{
    const int nnr = 7;
    MixingConfig cfg = make_cfg(2, 0.5, 0.8, 2);
    module_charge::init_mixing_dmr(&mixing, mdata, nnr, cfg);

    // dmr_nspin = 2 for nspin == 2
    EXPECT_EQ(mdata.length, static_cast<size_t>(nnr * 2));
    EXPECT_EQ(mdata.ndim_use, 0);
}

TEST_F(ChgDmrTest, InitMixingDmrPwThresholdAborts)
{
    const int nnr = 5;
    MixingConfig cfg = make_cfg(1, 0.5, 0.5, 1); // scf_thr_type == 1
    EXPECT_DEATH(module_charge::init_mixing_dmr(&mixing, mdata, nnr, cfg), "");
}

TEST_F(ChgDmrTest, InitMixingDmrNullMixingAborts)
{
    const int nnr = 5;
    MixingConfig cfg = make_cfg(1, 0.5, 0.5, 2);
    EXPECT_DEATH(module_charge::init_mixing_dmr(nullptr, mdata, nnr, cfg), "");
}

TEST_F(ChgDmrTest, InitMixingDmrNegativeNnrAborts)
{
    MixingConfig cfg = make_cfg(1, 0.5, 0.5, 2);
    EXPECT_DEATH(module_charge::init_mixing_dmr(&mixing, mdata, -1, cfg), "");
}

// ---------------------------------------------------------------------------
// mix_dmr nspin == 1
// ---------------------------------------------------------------------------

TEST_F(ChgDmrTest, MixDmrNspin1PlainStep)
{
    const int nnr = 4;
    MixingConfig cfg = make_cfg(1, 0.5, 0.5, 2);
    mixing.mixing_beta = 0.5;

    module_charge::init_mixing_dmr(&mixing, mdata, nnr, cfg);

    std::vector<double> dmr_in(nnr, 1.0);   // saved (previous) density matrix
    std::vector<double> dmr_out(nnr, 3.0);  // new density matrix from this step
    std::vector<double*> out_ptrs = {dmr_out.data()};
    std::vector<const double*> in_ptrs = {dmr_in.data()};

    module_charge::mix_dmr(out_ptrs, in_ptrs, nnr, &mixing, mdata, cfg);

    // plain mixing: out = in + beta * (out_new - in) = 1 + 0.5 * (3 - 1) = 2
    for (int i = 0; i < nnr; ++i)
    {
        EXPECT_NEAR(dmr_out[i], 2.0, 1e-12);
    }
}

// ---------------------------------------------------------------------------
// mix_dmr nspin == 2 (charge / magnetization channels with two betas)
// ---------------------------------------------------------------------------

TEST_F(ChgDmrTest, MixDmrNspin2ChargeConservation)
{
    const int nnr = 3;
    const double beta = 0.5;
    const double beta_mag = 0.8;
    MixingConfig cfg = make_cfg(2, beta, beta_mag, 2);
    mixing.mixing_beta = beta;

    module_charge::init_mixing_dmr(&mixing, mdata, nnr, cfg);

    // up/down saved and new
    std::vector<double> up_in(nnr, 1.0);
    std::vector<double> dn_in(nnr, 2.0);
    std::vector<double> up_out(nnr, 3.0);
    std::vector<double> dn_out(nnr, 4.0);

    std::vector<double*> out_ptrs = {up_out.data(), dn_out.data()};
    std::vector<const double*> in_ptrs = {up_in.data(), dn_in.data()};

    module_charge::mix_dmr(out_ptrs, in_ptrs, nnr, &mixing, mdata, cfg);

    // charge channel: c_save = 1+2 = 3, c_new = 3+4 = 7
    // c_mix = 3 + 0.5 * (7 - 3) = 5
    // mag channel: m_save = 1-2 = -1, m_new = 3-4 = -1
    // m_mix = -1 + 0.8 * (-1 - (-1)) = -1
    // up = 0.5 * (5 + (-1)) = 2
    // dn = 0.5 * (5 - (-1)) = 3
    for (int i = 0; i < nnr; ++i)
    {
        EXPECT_NEAR(up_out[i], 2.0, 1e-12);
        EXPECT_NEAR(dn_out[i], 3.0, 1e-12);
    }
}

// ---------------------------------------------------------------------------
// mix_dmr nspin == 4 (treated as single channel like nspin == 1)
// ---------------------------------------------------------------------------

TEST_F(ChgDmrTest, MixDmrNspin4SingleChannel)
{
    const int nnr = 3;
    MixingConfig cfg = make_cfg(4, 0.5, 0.5, 2);
    mixing.mixing_beta = 0.5;

    module_charge::init_mixing_dmr(&mixing, mdata, nnr, cfg);

    std::vector<double> dmr_in(nnr, 2.0);
    std::vector<double> dmr_out(nnr, 6.0);
    std::vector<double*> out_ptrs = {dmr_out.data()};
    std::vector<const double*> in_ptrs = {dmr_in.data()};

    module_charge::mix_dmr(out_ptrs, in_ptrs, nnr, &mixing, mdata, cfg);

    // out = 2 + 0.5 * (6 - 2) = 4
    for (int i = 0; i < nnr; ++i)
    {
        EXPECT_NEAR(dmr_out[i], 4.0, 1e-12);
    }
}

// ---------------------------------------------------------------------------
// boundary: empty partition (nnr == 0) allows null buffers
// ---------------------------------------------------------------------------

TEST_F(ChgDmrTest, MixDmrZeroNnrAllowsNullBuffers)
{
    const int nnr = 0;
    MixingConfig cfg = make_cfg(1, 0.5, 0.5, 2);
    module_charge::init_mixing_dmr(&mixing, mdata, nnr, cfg);

    // null buffers are legitimate when the rank owns no DMR elements
    std::vector<double*> out_ptrs = {nullptr};
    std::vector<const double*> in_ptrs = {nullptr};
    EXPECT_NO_THROW(module_charge::mix_dmr(out_ptrs, in_ptrs, nnr, &mixing, mdata, cfg));
}

// ---------------------------------------------------------------------------
// abort cases
// ---------------------------------------------------------------------------

TEST_F(ChgDmrTest, MixDmrNullMixingAborts)
{
    const int nnr = 4;
    MixingConfig cfg = make_cfg(1, 0.5, 0.5, 2);
    std::vector<double> buf(nnr, 0.0);
    std::vector<double*> out_ptrs = {buf.data()};
    std::vector<const double*> in_ptrs = {buf.data()};
    EXPECT_DEATH(module_charge::mix_dmr(out_ptrs, in_ptrs, nnr, nullptr, mdata, cfg), "");
}

TEST_F(ChgDmrTest, MixDmrBadNspinAborts)
{
    const int nnr = 4;
    MixingConfig cfg = make_cfg(3, 0.5, 0.5, 2); // nspin == 3 not supported
    std::vector<double> buf(nnr, 0.0);
    std::vector<double*> out_ptrs = {buf.data()};
    std::vector<const double*> in_ptrs = {buf.data()};
    EXPECT_DEATH(module_charge::mix_dmr(out_ptrs, in_ptrs, nnr, &mixing, mdata, cfg), "");
}

TEST_F(ChgDmrTest, MixDmrNullBufferWithPositiveNnrAborts)
{
    const int nnr = 4;
    MixingConfig cfg = make_cfg(1, 0.5, 0.5, 2);
    std::vector<double> buf(nnr, 0.0);
    // dmr_out is null while nnr > 0
    std::vector<double*> out_ptrs = {nullptr};
    std::vector<const double*> in_ptrs = {buf.data()};
    EXPECT_DEATH(module_charge::mix_dmr(out_ptrs, in_ptrs, nnr, &mixing, mdata, cfg), "");
}
