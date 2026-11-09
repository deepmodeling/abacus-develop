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
// mode 0 = overwrite a fixed filename -> no g{N} suffix even when istep >= 0.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, Mode0FixedNameIgnoresIstep)
{
    // PW wavefunction, binary, k1, spin up (ik0=0 < half_k=1 -> s1).
    const std::string fn = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/2, /*out_app_flag=*/0, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/-1);
    EXPECT_EQ(fn, "OUT.test/wfk1s1_pw.dat");
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
// mode 2 = per-step file -> the only mode that embeds the g{N} suffix.
// The directory routing is handled by the caller (set_global_dir).
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, Mode2EmbedsStep)
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
// iter >= 0 adds e{N} only in mode 2; modes 0 and 1 keep the fixed filename.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, IterIndexOnlyInMode2)
{
    // mode 0: neither g nor e.
    const std::string fn_mode0 = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/1, /*out_app_flag=*/0, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/3);
    EXPECT_EQ(fn_mode0, "OUT.test/wfk1s1_pw.txt");

    // mode 1: neither g nor e.
    const std::string fn_mode1 = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/1, /*out_app_flag=*/1, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/3);
    EXPECT_EQ(fn_mode1, "OUT.test/wfk1s1_pw.txt");

    // mode 2: both g and e present.
    const std::string fn_mode2 = filename_output(
        dir, "wf", "pw", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/1, /*out_app_flag=*/2, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/3);
    EXPECT_EQ(fn_mode2, "OUT.test/wfk1s1g2e3_pw.txt");
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
// The kinetic matrix T(k) is spin-independent, so its caller always passes
// nspin = 1 together with the property "t".  Gamma-only must not emit the
// k-point label, and multi-k must not double the "k" (tk1, not tkk1).
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, KineticPropertyUsesTWithoutDoubleK)
{
    const std::string fn_gamma = filename_output(
        dir, "t", "nao", ik_local, ik2iktot, /*nspin=*/1, nkstot,
        /*out_type=*/1, /*out_app_flag=*/1, /*gamma_only=*/true,
        /*istep=*/-1, /*iter=*/-1);
    EXPECT_EQ(fn_gamma, "OUT.test/t_nao.txt");

    const std::string fn_multi = filename_output(
        dir, "t", "nao", ik_local, ik2iktot, /*nspin=*/1, nkstot,
        /*out_type=*/1, /*out_app_flag=*/1, /*gamma_only=*/false,
        /*istep=*/-1, /*iter=*/-1);
    EXPECT_EQ(fn_multi, "OUT.test/tk1_nao.txt");
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

// ---------------------------------------------------------------------------
// The MO-basis Vxc matrix (property "vxcmo") and the NAO-basis one
// (property "vxck") must never resolve to the same filename.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, VxcMoDoesNotCollideWithVxck)
{
    const std::string fn_mo = filename_output(
        dir, "vxcmo", "nao", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/1, /*out_app_flag=*/0, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/-1);
    const std::string fn_nao = filename_output(
        dir, "vxck", "nao", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/1, /*out_app_flag=*/0, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/-1);

    EXPECT_EQ(fn_mo, "OUT.test/vxcmok1s1_nao.txt");
    EXPECT_EQ(fn_nao, "OUT.test/vxck1s1_nao.txt");
    EXPECT_NE(fn_mo, fn_nao);
}

// ---------------------------------------------------------------------------
// In mode 2 the MO-basis Vxc matrix carries the g{N} suffix like any other
// per-step matrix.
// ---------------------------------------------------------------------------
TEST_F(FilenameOutputTest, VxcMoMode2EmbedsStep)
{
    const std::string fn = filename_output(
        dir, "vxcmo", "nao", ik_local, ik2iktot, nspin, nkstot,
        /*out_type=*/1, /*out_app_flag=*/2, /*gamma_only=*/false,
        /*istep=*/1, /*iter=*/-1);
    EXPECT_EQ(fn, "OUT.test/vxcmok1s1g2_nao.txt");
}

#endif // TEST_FILENAME_CPP
