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
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, 4, 1, 1));
    // default kpar=1 on 4 processes: the silently diverging case of #8228
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, 1, 1, 1));
    // serial run is always fine
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(1, 1, 1, 1));
    // more processes than pools
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(8, 4, 1, 1));
}

TEST(DspKparLayoutValidTest, LCAOKparReset)
{
    // LCAO moves the user kpar into kpar_lcao and forces kpar to 1
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, 1, 4, 1));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, 1, 2, 1));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, 1, 1, 1));
}

TEST(DspKparLayoutValidTest, BndparExemption)
{
    // the bndpar+kpar layout that #7420 intended to allow
    EXPECT_TRUE(ModuleIO::dsp_kpar_layout_valid(4, 2, 1, 2));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, 1, 1, 2));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(8, 2, 1, 2));
}

TEST(DspKparLayoutValidTest, DegenerateValues)
{
    // non-positive factors are invalid layouts, and must be rejected by plain
    // arithmetic without any division so no undefined behavior can occur
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, 0, 0, 1));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, 1, 1, 0));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(4, -1, -1, 1));
    EXPECT_FALSE(ModuleIO::dsp_kpar_layout_valid(0, 1, 1, 1));
}
