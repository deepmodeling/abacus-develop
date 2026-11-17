/**
 * @file test_nacv_tools_io.cpp
 * @brief Exercise the real NAC eigenstate writer, including Release builds.
 */
#include "../nacv_tools.h"

#include <gtest/gtest.h>

#include <complex>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace
{
template <typename T> void expect_file_values(const std::string &filename, const std::vector<T> &expected)
{
    std::ifstream input(filename);
    ASSERT_TRUE(input.is_open()) << filename;
    for (const T &reference : expected)
    {
        T actual{};
        input >> actual;
        ASSERT_FALSE(input.fail()) << filename;
        EXPECT_EQ(actual, reference);
    }
    input >> std::ws;
    EXPECT_TRUE(input.eof()) << filename;
}

/** @brief Give each writer test its own filenames in the gtest scratch directory. */
class NacvIoTest : public ::testing::Test
{
  protected:
    std::string directory;
    std::string label;
    std::string suffix;
    std::string energy_filename;
    std::vector<std::string> amplitude_filenames;

    void SetUp() override
    {
        directory = ::testing::TempDir();
        label = "nacv_";
        label += ::testing::UnitTest::GetInstance()->current_test_info()->name();
        suffix = "_g2";
        energy_filename = directory + "Excitation_Energy_" + label + suffix + ".dat";
        for (int rank = 0; rank < 3; ++rank)
        {
            const std::string filename =
                directory + "Excitation_Amplitude_" + label + suffix + "_" + std::to_string(rank + 1) + ".dat";
            amplitude_filenames.push_back(filename);
        }
    }

    void TearDown() override
    {
        std::remove(energy_filename.c_str());
        for (const std::string &filename : amplitude_filenames)
        {
            std::remove(filename.c_str());
        }
    }
};
} // namespace

TEST_F(NacvIoTest, RealStatesRoundTripAtDoublePrecision)
{
    const std::vector<double> energies{1.1234567890123457, 2.9876543210987654};
    const std::vector<double> amplitudes{0.12345678901234568,
        -0.9876543210987654,
        1.2345678901234568e-9,
        -0.3333333333333333,
        0.7777777777777778,
        0.0};
    const double *const energy_data = energies.data();
    const double *const amplitude_data = amplitudes.data();
    LR::write_nac_states(
        directory,
        label,
        suffix,
        0,
        energy_data,
        amplitude_data,
        3,
        2);
    expect_file_values(energy_filename, energies);
    expect_file_values(amplitude_filenames[0], amplitudes);
}

TEST_F(NacvIoTest, ComplexStatesKeepRealAndImaginaryComponents)
{
    const std::vector<double> energies{0.12345678901234568, 0.9876543210987654};
    const std::vector<std::complex<double>> amplitudes{{0.12345678901234568, -0.9876543210987654},
        {-0.3333333333333333, 0.7777777777777778}, {1.2345678901234568e-9, -2.345678901234568e-9}, {0.0, 0.0}};
    const double *const energy_data = energies.data();
    const std::complex<double> *const amplitude_data = amplitudes.data();
    LR::write_nac_states(
        directory,
        label,
        suffix,
        0,
        energy_data,
        amplitude_data,
        2,
        2);
    expect_file_values(energy_filename, energies);
    expect_file_values(amplitude_filenames[0], amplitudes);
}

TEST_F(NacvIoTest, NonwriterCreatesOnlyItsLocalAmplitudeFile)
{
    const std::vector<double> energies{1.5};
    const std::vector<double> amplitudes{0.25, -0.75};
    const double *const energy_data = energies.data();
    const double *const amplitude_data = amplitudes.data();
    LR::write_nac_states(
        directory,
        label,
        suffix,
        1,
        energy_data,
        amplitude_data,
        2,
        1);
    std::ifstream energy_file(energy_filename);
    EXPECT_FALSE(energy_file.is_open());
    expect_file_values(amplitude_filenames[1], amplitudes);
    std::ifstream other_rank(amplitude_filenames[0]);
    EXPECT_FALSE(other_rank.is_open());
}

TEST_F(NacvIoTest, NonwriterPreservesTheWriterEnergyFile)
{
    {
        std::ofstream output(energy_filename);
        ASSERT_TRUE(output.is_open());
        output << "0.125\n";
    }
    const std::vector<double> energies{1.5};
    const std::vector<double> amplitudes{-0.5};
    const double *const energy_data = energies.data();
    const double *const amplitude_data = amplitudes.data();
    LR::write_nac_states(
        directory,
        label,
        suffix,
        2,
        energy_data,
        amplitude_data,
        1,
        1);
    const std::vector<double> original_energies{0.125};
    expect_file_values(energy_filename, original_energies);
    expect_file_values(amplitude_filenames[2], amplitudes);
}

TEST_F(NacvIoTest, EmptyLocalAmplitudeWindowProducesAnEmptyRecord)
{
    const double energies[]{0.5, 1.0};
    const double unused_amplitude = 7.0;
    LR::write_nac_states(
        directory,
        label,
        suffix,
        1,
        energies,
        &unused_amplitude,
        0,
        2);
    const std::vector<double> empty;
    expect_file_values(amplitude_filenames[1], empty);
    std::ifstream energy_file(energy_filename);
    EXPECT_FALSE(energy_file.is_open());
}
