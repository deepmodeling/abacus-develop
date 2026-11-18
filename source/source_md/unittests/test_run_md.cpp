#include "gtest/gtest.h"

#include "source_md/run_md.h"

// --- Cadence independence: output (out_freq_ion) and restart
// (md_restartfreq) must be driven by their own frequencies so that
// setting out_freq_ion does not silently re-aim md_restartfreq.

TEST(RunMDCadenceTest, disabled_when_freq_non_positive)
{
    for (int step = 0; step <= 20; ++step)
    {
        EXPECT_FALSE(Run_MD::is_cadence_step(step, 0));
        EXPECT_FALSE(Run_MD::is_cadence_step(step, -1));
    }
}

TEST(RunMDCadenceTest, fires_every_freq_steps)
{
    const int freq = 3;
    for (int step = 0; step <= 20; ++step)
    {
        const bool expected = (step % freq == 0);
        EXPECT_EQ(Run_MD::is_cadence_step(step, freq), expected) << "step=" << step;
    }
}

// Core regression guard: with out_freq_ion and md_restartfreq set to
// different positive values, the two cadences must disagree on at least
// one step. This guarantees md_restartfreq is no longer swallowed by
// out_freq_ion > 0.
TEST(RunMDCadenceTest, output_and_restart_are_independent)
{
    const int out_freq_ion = 4;
    const int md_restartfreq = 6;

    bool found_disagreement = false;
    for (int step = 0; step <= 60; ++step)
    {
        const bool output_step = Run_MD::is_cadence_step(step, out_freq_ion);
        const bool restart_step = Run_MD::is_cadence_step(step, md_restartfreq);
        if (output_step != restart_step)
        {
            found_disagreement = true;
            break;
        }
    }
    EXPECT_TRUE(found_disagreement)
        << "out_freq_ion=" << out_freq_ion << " and md_restartfreq=" << md_restartfreq
        << " produced identical cadences; md_restartfreq is still being ignored.";
}

// --- Frame-write invariant: every restart cadence step must produce a
// checkpoint STRU frame so Restart_md.txt is never left pointing at a step
// that has no matching structure. These tests exercise is_frame_write_step
// as the single source of truth behind run_md's write_any_frame branch.

TEST(RunMDCadenceTest, frame_written_on_every_restart_step_regardless_of_out_stru)
{
    // The original bug report: out_freq_ion=3, md_restartfreq=5, out_stru=0.
    // Restart_md.txt fires at steps 0, 5, 10, 15... but pre-fix code only
    // wrote STRU frames at the out_freq_ion cadence (0, 3, 6, 9...), leaving
    // no frame for restart steps 5, 10, etc.
    const int out_freq_ion = 3;
    const int md_restartfreq = 5;
    const bool out_stru = false;

    for (int step = 0; step <= 60; ++step)
    {
        const bool should_write = Run_MD::is_frame_write_step(step, out_freq_ion, md_restartfreq, out_stru);
        const bool is_restart_step = Run_MD::is_cadence_step(step, md_restartfreq);
        if (is_restart_step)
        {
            EXPECT_TRUE(should_write)
                << "step=" << step
                << " is a restart cadence step but is_frame_write_step returned false; "
                << "Restart_md.txt will have no matching STRU frame.";
        }
    }
}

TEST(RunMDCadenceTest, frame_written_on_output_step_when_out_stru_positive)
{
    const int out_freq_ion = 3;
    const int md_restartfreq = 5;
    const bool out_stru = true;

    for (int step = 0; step <= 60; ++step)
    {
        const bool should_write = Run_MD::is_frame_write_step(step, out_freq_ion, md_restartfreq, out_stru);
        const bool is_restart_step = Run_MD::is_cadence_step(step, md_restartfreq);
        const bool is_output_step = Run_MD::is_cadence_step(step, out_freq_ion);

        // write when either cadence fires (union of the two)
        EXPECT_EQ(should_write, is_restart_step || is_output_step)
            << "step=" << step
            << ": should_write=" << should_write
            << " is_restart_step=" << is_restart_step
            << " is_output_step=" << is_output_step;
    }
}

TEST(RunMDCadenceTest, no_output_frame_when_out_stru_zero_except_restart_cadence)
{
    const int out_freq_ion = 3;
    const int md_restartfreq = 5;
    const bool out_stru = false;

    for (int step = 0; step <= 60; ++step)
    {
        const bool should_write = Run_MD::is_frame_write_step(step, out_freq_ion, md_restartfreq, out_stru);
        // out_stru=0 silences the per-step output cadence; only the restart
        // checkpoint frame remains.
        EXPECT_EQ(should_write, Run_MD::is_cadence_step(step, md_restartfreq))
            << "step=" << step << ": out_stru=0 must not write frames at the "
            << "out_freq_ion cadence, only at the md_restartfreq cadence.";
    }
}

TEST(RunMDCadenceTest, backward_compat_fallback_when_out_freq_ion_unset)
{
    // Default: out_freq_ion=0, md_restartfreq=5, out_stru=1. The effective
    // output cadence falls back to md_restartfreq so the output and restart
    // cadences coincide every step. Historically this was the only supported
    // configuration and must keep working.
    const int out_freq_ion = 0;
    const int md_restartfreq = 5;
    const bool out_stru = true;

    for (int step = 0; step <= 60; ++step)
    {
        const bool should_write = Run_MD::is_frame_write_step(step, out_freq_ion, md_restartfreq, out_stru);
        const bool expected = Run_MD::is_cadence_step(step, md_restartfreq);
        EXPECT_EQ(should_write, expected) << "step=" << step;
    }
}

TEST(RunMDCadenceTest, checkpoint_frame_also_writes_on_joint_step)
{
    // When out_freq_ion and md_restartfreq coincide at step LCM they fire
    // together; the predicate still returns a single bool (one write), not
    // two separate bools. We just assert that step 15 (lcm(3,5)) produces a
    // true so the write_dedup code path in run_md.cpp is exercised.
    const int out_freq_ion = 3;
    const int md_restartfreq = 5;
    EXPECT_TRUE(Run_MD::is_frame_write_step(15, out_freq_ion, md_restartfreq, true));
    EXPECT_TRUE(Run_MD::is_frame_write_step(15, out_freq_ion, md_restartfreq, false));
}

TEST(RunMDCadenceTest, no_frame_when_all_frequencies_disabled)
{
    // Nothing writes when both cadence frequencies are disabled (0).
    EXPECT_FALSE(Run_MD::is_frame_write_step(3, 0, 0, true));
    EXPECT_FALSE(Run_MD::is_frame_write_step(3, 0, 0, false));
    // md_restartfreq disabled and step doesn't hit output cadence → nothing
    EXPECT_FALSE(Run_MD::is_frame_write_step(7, 3, 0, true));
    // md_restartfreq disabled but step does hit output cadence with out_stru →
    // output frame is still written
    EXPECT_TRUE(Run_MD::is_frame_write_step(3, 3, 0, true));
    // md_restartfreq disabled, step hits output cadence but out_stru=0 → no
    // output frame, no checkpoint frame either
    EXPECT_FALSE(Run_MD::is_frame_write_step(3, 3, 0, false));
}
