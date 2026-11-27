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
 * @brief Ordinary structure cadence must be independent of restart output.
 */
TEST(RunMdFinalStepTest, StructureOutputRequiresItsOwnPositiveCadence)
{
    EXPECT_FALSE(Run_MD::is_structure_output_step(0, 0, true));
    EXPECT_FALSE(Run_MD::is_structure_output_step(0, 3, false));
    EXPECT_TRUE(Run_MD::is_structure_output_step(0, 3, true));
    EXPECT_FALSE(Run_MD::is_structure_output_step(2, 3, true));
    EXPECT_TRUE(Run_MD::is_structure_output_step(3, 3, true));
    EXPECT_TRUE(Run_MD::is_cadence_step(0, 5));
    EXPECT_FALSE(Run_MD::is_cadence_step(0, 0));
    EXPECT_FALSE(Run_MD::is_cadence_step(3, 5));
}

} // namespace
