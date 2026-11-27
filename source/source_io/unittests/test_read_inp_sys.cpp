#include "gtest/gtest.h"

#include "source_io/module_parameter/read_input_tool.h"

/************************************************
 *  unit tests of dsp_kpar_layout_valid
 *
 *  The helper encodes the DSP vendor library (mtblas/C_PACK) constraint
 *  checked by the kpar input validation under __DSP: the effective k-point
 *  pool count (kpar for PW, kpar_lcao for LCAO, the other being 1) times
 *  bndpar must equal the number of MPI processes. See issue #8228.
 ***********************************************/

TEST(DspKparLayoutValidTest, PWKparPools)
{
    // PW: kpar keeps the user value, kpar_lcao stays 1
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, "pw", 4, 1, 1, "ksdft"));
    // default kpar=1 on 4 processes: the silently diverging case of #8228
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, "pw", 1, 1, 1, "ksdft"));
    // serial run is always fine
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(1, "pw", 1, 1, 1, "ksdft"));
    // more processes than pools
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(8, "pw", 4, 1, 1, "ksdft"));
}

TEST(DspKparLayoutValidTest, LCAOKparReset)
{
    // LCAO moves the user kpar into kpar_lcao and forces kpar to 1
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, "lcao", 1, 4, 1, "ksdft"));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, "lcao", 1, 2, 1, "ksdft"));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, "lcao", 1, 1, 1, "ksdft"));
}

TEST(DspKparLayoutValidTest, EffectiveKparSelection)
{
    // Regression guard for the removed std::max(kpar, kpar_lcao):
    // the field NOT used by the current basis_type must be ignored,
    // even when it is larger than the effective one.

    // LCAO with a stale kpar=8: old std::max took 8 and rejected a
    // perfectly valid LCAO layout; the explicit selection accepts it.
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, "lcao", 8, 4, 1, "ksdft"));

    // PW with a stale kpar_lcao=8: old std::max took 8 and wrongly
    // accepted an invalid PW layout (4 pools on 8 processes).
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(8, "pw", 4, 8, 1, "ksdft"));

    // unknown basis_type falls back to the PW field
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, "lcao_in_pw", 4, 1, 1, "ksdft"));
}

TEST(DspKparLayoutValidTest, BndparExemption)
{
    // the bndpar+kpar layout that #7420 intended to allow
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, "pw", 2, 1, 2, "ksdft"));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, "pw", 1, 1, 2, "ksdft"));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(8, "pw", 2, 1, 2, "ksdft"));
    // LCAO variant: kpar_lcao * bndpar == nproc
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, "lcao", 1, 2, 2, "ksdft"));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, "lcao", 1, 1, 2, "ksdft"));
}

TEST(DspKparLayoutValidTest, EsolverExemption)
{
    // dp / nep never run an eigensolver and never touch mtblas/C_PACK:
    // any layout must be accepted, however mismatched
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, "pw",   2, 1, 3, "dp"));
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(7, "lcao", 1, 2, 3, "nep"));
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(0, "pw",   0, 0, 0, "dp"));
    // case-insensitive exemption
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, "pw", 2, 1, 3, "DP"));
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, "pw", 2, 1, 3, "Nep"));
    // non-ML esolvers are still checked
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, "pw", 1, 1, 2, "ksdft"));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, "pw", 1, 1, 2, "sdft"));
}

TEST(DspKparLayoutValidTest, DegenerateValues)
{
    // non-positive factors are invalid layouts, rejected by plain
    // comparisons without any division so no UB can occur
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, "pw",   0, 0, 1, "ksdft"));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, "pw",   1, 1, 0, "ksdft"));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, "pw",  -1, -1, 1, "ksdft"));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(0, "pw",   1, 1, 1, "ksdft"));
    // the 0*1 == 0 hole: nproc=0 with all-zero factors must NOT pass
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(0, "pw",   0, 0, 1, "ksdft"));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(0, "lcao", 0, 0, 0, "ksdft"));
    // negative product accidentally matching a negative nproc must not pass
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(-4, "pw", -2, 1, 2, "ksdft"));
}

//---------------------------------------------------------------
// sub-helpers: unit-level coverage of the two extracted decisions
//---------------------------------------------------------------

TEST(DspLayoutCheckNeededTest, EsolverType)
{
    EXPECT_FALSE(ModuleIO::dsp_layout_check_needed("dp"));
    EXPECT_FALSE(ModuleIO::dsp_layout_check_needed("nep"));
    EXPECT_FALSE(ModuleIO::dsp_layout_check_needed("DP"));
    EXPECT_FALSE(ModuleIO::dsp_layout_check_needed("NEP"));
    EXPECT_TRUE(ModuleIO::dsp_layout_check_needed("ksdft"));
    EXPECT_TRUE(ModuleIO::dsp_layout_check_needed("sdft"));
    EXPECT_TRUE(ModuleIO::dsp_layout_check_needed(""));
    // trailing whitespace is NOT normalized: input parsing lowercases
    // esolver_type upstream; an unnormalized value is checked strictly
    EXPECT_TRUE(ModuleIO::dsp_layout_check_needed("dp "));
}

TEST(DspEffectiveKparTest, BasisTypeSelection)
{
    EXPECT_EQ(ModuleIO::dsp_effective_kpar("pw",   4, 1), 4);
    EXPECT_EQ(ModuleIO::dsp_effective_kpar("lcao", 1, 4), 4);
    // invariant after reset_value: exactly one of the two is 1
    EXPECT_EQ(ModuleIO::dsp_effective_kpar("lcao", 1, 1), 1);
    // unknown basis types take the PW field, matching reset_value,
    // which only reroutes kpar when basis_type == "lcao"
    EXPECT_EQ(ModuleIO::dsp_effective_kpar("", 3, 1), 3);
}
