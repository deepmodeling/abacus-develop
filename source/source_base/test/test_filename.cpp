#ifndef TEST_FILENAME_CPP
#define TEST_FILENAME_CPP

#include <string>
#include <vector>
#include <gtest/gtest.h>
#include "source_base/module_out/filename.h"

// Keep local names short; do not pollute.
using ModuleIO::filename_output;

// ---------------------------------------------------------------------------
// Fixture: minimal ik2iktot mapping for a 2-k-point pool.
// ---------------------------------------------------------------------------
class FilenameOutputTest : public ::testing::Test
{
protected:
    const std::vector<int> ik2iktot = {0, 1};
    const int ik_local = 0;
    const int nkstot = 2;
    const int nspin = 2;
    const std::string dir = "OUT.test/";
};

// ---------------------------------------------------------------------------
// mode 0 = no append -> per-step g{N} suffix included.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, Mode0IstepAddsG)
{
    // PW wavefunction, binary, k1, spin up (ik0=0 < half_k=1 -> s1).
    const std::string fn = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/2, /*out_app_flag=*/0, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/-1);
    EXPECT_EQ(fn, "OUT.test/wfk1s1g2_pw.dat");
}

// ---------------------------------------------------------------------------
// mode 1 = append -> fixed name, no g suffix even when istep >= 0.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, Mode1FixedNameNoG)
{
    const std::string fn = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/2, /*out_app_flag=*/1, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/-1);
    EXPECT_EQ(fn, "OUT.test/wfk1s1_pw.dat");
}

// ---------------------------------------------------------------------------
// mode 2 = write to subdirectory -> same naming as mode 0 (g suffix kept).
// The directory routing is handled by the caller (set_global_dir), so
// filename_output itself treats modes 0 and 2 identically.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, Mode2SameNamingAsMode0)
{
    const std::string fn = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/2, /*out_app_flag=*/2, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/-1);
    EXPECT_EQ(fn, "OUT.test/wfk1s1g2_pw.dat");
}

// ---------------------------------------------------------------------------
// istep < 0 (SCF single-point) -> no g suffix in any mode.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, SinglePointNoG)
{
    const std::string fn = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/2, /*out_app_flag=*/0, /*gamma_only=*/false,
        /*istep=*/-1, /*iter=*/-1);
    EXPECT_EQ(fn, "OUT.test/wfk1s1_pw.dat");
}

// ---------------------------------------------------------------------------
// iter >= 0 adds e{N} in modes 0/2, but is suppressed in mode 1 too.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, IterIndexRespectsAppendFlag)
{
    // mode 0: both g and e present.
    const std::string fn_mode0 = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/1, /*out_app_flag=*/0, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/3);
    EXPECT_EQ(fn_mode0, "OUT.test/wfk1s1g2e3_pw.txt");

    // mode 1: neither g nor e.
    const std::string fn_mode1 = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/1, /*out_app_flag=*/1, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/3);
    EXPECT_EQ(fn_mode1, "OUT.test/wfk1s1_pw.txt");
}

// ---------------------------------------------------------------------------
// gamma_only suppresses the k-point label for non-hk/sk/tk properties.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, GammaOnlyNoKLabel)
{
    const std::string fn = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/2, /*out_app_flag=*/1, /*gamma_only=*/true,
        /*istep=*/-1, /*iter=*/-1);
    EXPECT_EQ(fn, "OUT.test/wfs1_pw.dat");
}

// ---------------------------------------------------------------------------
// nspin=4 (noncollinear SOC) -> spin block is s4, no half-k split.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, NoncollinearSpin4)
{
    const std::string fn = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, /*nspin=*/4, nkstot,
        /*out_type=*/2, /*out_app_flag=*/2, /*gamma_only=*/false,
        /*istep=*/0, /*iter=*/-1);
    EXPECT_EQ(fn, "OUT.test/wfk1s4g1_pw.dat");
}

#endif // TEST_FILENAME_CPP
