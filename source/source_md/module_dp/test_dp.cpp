#include "source_md/module_dp/dp.h"

#include "source_md/test/setcell.h"

#include "gtest/gtest.h"

#include <cmath>
#include <fstream>
#include <string>

TEST(DPPotential, UnitCellCalculationOrMissingBackend)
{
    const std::string source_file = __FILE__;
    const std::string marker = "/source/source_md/";
    const std::size_t marker_pos = source_file.find(marker);
    ASSERT_NE(marker_pos, std::string::npos);
    const std::string model_file = source_file.substr(0, marker_pos) + "/tests/PP_ORB/Al-SCAN.pb";

    UnitCell cell;
    Setcell::setupcell(cell);
    cell.atoms[0].label = "Al";
    Input_para input;
    Setcell::parameters(input);
    const std::string output_dir = testing::TempDir();
    std::ofstream running_log;
    ModuleMD::DP potential(model_file, output_dir, running_log);
    potential.initialize(cell, input);

    double energy = 0.0;
    ModuleBase::matrix force;
    ModuleBase::matrix stress(3, 3);
#ifdef __DPMD
    potential.compute(cell, energy, force, stress);
    EXPECT_TRUE(std::isfinite(energy));
    for (int i = 0; i < cell.nat; ++i)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            EXPECT_TRUE(std::isfinite(force(i, direction)));
        }
    }
#else
    EXPECT_EXIT(potential.compute(cell, energy, force, stress), ::testing::ExitedWithCode(1), "");
#endif
}
