#include "source_lcao/module_lr/utils/lr_util_print.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <unistd.h>

#include <gtest/gtest.h>

TEST(LRUtilPrint, ExcitationEnergyFileRoundTripsDoublePrecision)
{
    constexpr int nstates = 2;
    // At 16 significant digits, the next double above 1.0 rounds to 1.0.
    const std::array<double, nstates> energies = {std::nextafter(1.0, 2.0), 1e-12};
    const std::string filename = testing::TempDir() + "lr_excitation_energies_" + std::to_string(getpid()) + ".dat";
    EXPECT_TRUE(LR_Util::write_excitation_energy_file(filename, energies.data(), nstates));

    std::array<double, nstates> restored{};
    std::ifstream input(filename);
    EXPECT_TRUE(static_cast<bool>(input >> restored[0] >> restored[1]));
    input.close();
    EXPECT_EQ(std::remove(filename.c_str()), 0);
    EXPECT_EQ(restored, energies);
}
