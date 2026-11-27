#include "gtest/gtest.h"

#include "source_pw/module_pwdft/dftu_base_io.h"

#include <string>

#ifdef __MPI
#include "mpi.h"
#endif

// Focused tests for the occupation-matrix filename helpers and gating logic
// implemented in dftu_base_io.cpp.

TEST(DFTUBaseIoTest, IonStepGateZeroFreqDisablesNumberedFiles)
{
    const DFTU_BASE::OccmatOutputCfg cfg{0, 10, 10, true, 1};
    EXPECT_FALSE(DFTU_BASE::is_ion_step_output_step(0, cfg));
    EXPECT_FALSE(DFTU_BASE::is_ion_step_output_step(1, cfg));
    EXPECT_FALSE(DFTU_BASE::is_ion_step_output_step(2, cfg));
}

TEST(DFTUBaseIoTest, IonStepGateDivisibleIonStepsOnly)
{
    const DFTU_BASE::OccmatOutputCfg cfg{2, 10, 10, true, 1};
    EXPECT_TRUE(DFTU_BASE::is_ion_step_output_step(0, cfg));
    EXPECT_FALSE(DFTU_BASE::is_ion_step_output_step(1, cfg));
    EXPECT_TRUE(DFTU_BASE::is_ion_step_output_step(2, cfg));
    EXPECT_FALSE(DFTU_BASE::is_ion_step_output_step(3, cfg));
    EXPECT_TRUE(DFTU_BASE::is_ion_step_output_step(4, cfg));
}

TEST(DFTUBaseIoTest, ElecGatePeriodicTrigger)
{
    const DFTU_BASE::OccmatOutputCfg cfg{1, 2, 10, true, 1};
    EXPECT_FALSE(DFTU_BASE::is_elec_snapshot_trigger(1, false, cfg));
    EXPECT_TRUE(DFTU_BASE::is_elec_snapshot_trigger(2, false, cfg));
    EXPECT_FALSE(DFTU_BASE::is_elec_snapshot_trigger(3, false, cfg));
    EXPECT_TRUE(DFTU_BASE::is_elec_snapshot_trigger(4, false, cfg));
}

TEST(DFTUBaseIoTest, ElecGateScfNmaxTrigger)
{
    const DFTU_BASE::OccmatOutputCfg cfg{1, 10, 10, true, 1};
    EXPECT_FALSE(DFTU_BASE::is_elec_snapshot_trigger(1, false, cfg));
    EXPECT_FALSE(DFTU_BASE::is_elec_snapshot_trigger(9, false, cfg));
    EXPECT_TRUE(DFTU_BASE::is_elec_snapshot_trigger(10, false, cfg));
}

TEST(DFTUBaseIoTest, ElecGateConvergenceTrigger)
{
    const DFTU_BASE::OccmatOutputCfg cfg{1, 10, 10, true, 1};
    EXPECT_TRUE(DFTU_BASE::is_elec_snapshot_trigger(3, true, cfg));
    EXPECT_TRUE(DFTU_BASE::is_elec_snapshot_trigger(7, true, cfg));
}

TEST(DFTUBaseIoTest, IonStepFilenameStartsFromOne)
{
    const std::string fn0 = DFTU_BASE::gen_ion_step_occ_mat_filename("OUT.ABACUS/", 0);
    EXPECT_EQ(fn0, "OUT.ABACUS/occ_matg1.txt");

    const std::string fn4 = DFTU_BASE::gen_ion_step_occ_mat_filename("OUT.ABACUS/", 4);
    EXPECT_EQ(fn4, "OUT.ABACUS/occ_matg5.txt");
}

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
