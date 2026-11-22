#include "source_md/run_md.h"

#include "gtest/gtest.h"

// --- Final step predicate tests (local additions to upstream's consolidate)
// ---
// These only exercise the inline helpers declared in run_md.h, so the target
// can stay free of the heavy MD/cell link dependencies.

namespace
{

/**
 * @brief The global step (local step + restart offset) md_nstep - 1 is the
 * only ionic step that counts as final, so the output layer can always
 * persist the final snapshot.
 */
TEST(RunMdFinalStepTest, IdentifiesTheFinalIonicStep)
{
    EXPECT_FALSE(Run_MD::is_last_step(0, 4));
    EXPECT_FALSE(Run_MD::is_last_step(2, 4));
    EXPECT_TRUE(Run_MD::is_last_step(3, 4));
    // Restart: continue from step 5 up to md_nstep=6; local step 1 is final.
    EXPECT_TRUE(Run_MD::is_last_step(5, 6));
    // Single-step run.
    EXPECT_TRUE(Run_MD::is_last_step(0, 1));
    // Invalid input never yields a final step.
    EXPECT_FALSE(Run_MD::is_last_step(0, 0));
}

/**
 * @brief The supplementary final frame must not duplicate a frame that the
 * MD loop already wrote on the output or restart cadence.
 */
TEST(RunMdFinalStepTest, SupplementaryFinalFrameSkipsStepsAlreadyOnCadence)
{
    // Final step off the output cadence: the loop skips it, a frame is needed.
    EXPECT_FALSE(Run_MD::is_frame_write_step(3, 2, 0, true));
    // Final step on the output cadence: already written by the loop.
    EXPECT_TRUE(Run_MD::is_frame_write_step(4, 2, 0, true));
    // Final step on the restart cadence: already written by the loop.
    EXPECT_TRUE(Run_MD::is_frame_write_step(4, 2, 4, true));
    // md_nstep=1: step 0 is final and is inherently a cadence step.
    EXPECT_TRUE(Run_MD::is_last_step(0, 1));
    EXPECT_TRUE(Run_MD::is_frame_write_step(0, 2, 0, true));
    // Both frequencies disabled: the loop writes nothing, a final frame is needed.
    EXPECT_FALSE(Run_MD::is_frame_write_step(0, 0, 0, true));
}

} // namespace
