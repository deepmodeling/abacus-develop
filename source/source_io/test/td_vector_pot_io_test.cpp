#include "source_io/module_efield/td_vector_pot_io.h"

#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace
{

class TDVectorPotIOTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        clean_output_file();
    }

    void TearDown() override
    {
        clean_output_file();
    }

    void clean_output_file() const
    {
        std::remove(output_path_.c_str());
    }

    std::vector<std::string> read_lines() const
    {
        std::ifstream input(output_path_.c_str());
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(input, line))
        {
            lines.push_back(line);
        }
        return lines;
    }

    const std::string output_prefix_ = "td_vector_pot_io_test_";
    const std::string output_path_ = output_prefix_ + "vector_pot.txt";
};

TEST_F(TDVectorPotIOTest, FreshCalculationTruncatesAndWritesHeader)
{
    {
        std::ofstream output(output_path_.c_str());
        output << "old data\n";
    }

    ModuleIO::prepare_td_vector_pot_output(output_prefix_, false);
    ModuleIO::write_td_vector_pot(output_prefix_, 0, ModuleBase::Vector3<double>(1.0, 2.0, 3.0));
    ModuleIO::write_td_vector_pot(output_prefix_, 1, ModuleBase::Vector3<double>(4.0, 5.0, 6.0));

    const std::vector<std::string> lines = read_lines();
    ASSERT_EQ(lines.size(), 3U);
    EXPECT_EQ(lines[0].find("#istep"), 0U);
    EXPECT_EQ(lines[1].find("1"), 0U);
    EXPECT_EQ(lines[2].find("2"), 0U);

    const std::vector<ModuleBase::Vector3<double>> vector_potentials = ModuleIO::read_td_vector_pot(output_prefix_);
    ASSERT_EQ(vector_potentials.size(), 2U);
    EXPECT_DOUBLE_EQ(vector_potentials[0][0], 1.0);
    EXPECT_DOUBLE_EQ(vector_potentials[0][1], 2.0);
    EXPECT_DOUBLE_EQ(vector_potentials[0][2], 3.0);
    EXPECT_DOUBLE_EQ(vector_potentials[1][0], 4.0);
    EXPECT_DOUBLE_EQ(vector_potentials[1][1], 5.0);
    EXPECT_DOUBLE_EQ(vector_potentials[1][2], 6.0);
}

TEST_F(TDVectorPotIOTest, RestartPreservesExistingSamples)
{
    {
        std::ofstream output(output_path_.c_str());
        output << "#istep A_x A_y A_z\n";
        output << "7 1.0 2.0 3.0\n";
    }

    ModuleIO::prepare_td_vector_pot_output(output_prefix_, true);
    ModuleIO::write_td_vector_pot(output_prefix_, 7, ModuleBase::Vector3<double>(4.0, 5.0, 6.0));

    const std::vector<std::string> lines = read_lines();
    ASSERT_EQ(lines.size(), 3U);
    EXPECT_EQ(lines[0], "#istep A_x A_y A_z");
    EXPECT_EQ(lines[1], "7 1.0 2.0 3.0");
    EXPECT_EQ(lines[2].find("8"), 0U);
}

TEST_F(TDVectorPotIOTest, RestartCreatesHeaderForMissingOrEmptyFile)
{
    ModuleIO::prepare_td_vector_pot_output(output_prefix_, true);
    std::vector<std::string> lines = read_lines();
    ASSERT_EQ(lines.size(), 1U);
    EXPECT_EQ(lines[0].find("#istep"), 0U);

    {
        std::ofstream output(output_path_.c_str(), std::ofstream::out);
    }
    ModuleIO::prepare_td_vector_pot_output(output_prefix_, true);
    lines = read_lines();
    ASSERT_EQ(lines.size(), 1U);
    EXPECT_EQ(lines[0].find("#istep"), 0U);
}

TEST_F(TDVectorPotIOTest, ReaderSkipsBlankAndCommentLinesAndIgnoresLabels)
{
    {
        std::ofstream output(output_path_.c_str());
        output << "  # comment\n";
        output << "\n";
        output << "19 1.5 -2.5 3.5\n";
        output << "42 4.5 5.5 -6.5\n";
    }

    const std::vector<ModuleBase::Vector3<double>> vector_potentials = ModuleIO::read_td_vector_pot(output_prefix_);
    ASSERT_EQ(vector_potentials.size(), 2U);
    EXPECT_DOUBLE_EQ(vector_potentials[0][0], 1.5);
    EXPECT_DOUBLE_EQ(vector_potentials[0][1], -2.5);
    EXPECT_DOUBLE_EQ(vector_potentials[0][2], 3.5);
    EXPECT_DOUBLE_EQ(vector_potentials[1][0], 4.5);
    EXPECT_DOUBLE_EQ(vector_potentials[1][1], 5.5);
    EXPECT_DOUBLE_EQ(vector_potentials[1][2], -6.5);
}

TEST_F(TDVectorPotIOTest, ReaderRejectsEmptyAndMalformedFiles)
{
    {
        std::ofstream output(output_path_.c_str());
        output << "# no samples\n";
    }
    EXPECT_EXIT(ModuleIO::read_td_vector_pot(output_prefix_), testing::ExitedWithCode(1), "");

    {
        std::ofstream output(output_path_.c_str(), std::ofstream::out);
        output << "1 2.0 invalid 4.0\n";
    }
    EXPECT_EXIT(ModuleIO::read_td_vector_pot(output_prefix_), testing::ExitedWithCode(1), "");

    {
        std::ofstream output(output_path_.c_str(), std::ofstream::out);
        output << "1 2.0 3.0 4.0 extra\n";
    }
    EXPECT_EXIT(ModuleIO::read_td_vector_pot(output_prefix_), testing::ExitedWithCode(1), "");
}

} // namespace
