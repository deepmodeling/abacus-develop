#include "source_io/module_efield/td_efield_io.h"

#include "source_base/constants.h"
#include "source_estate/module_pot/td_field_manager.h"
#include "source_io/module_parameter/input_parameter.h"

#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace
{

class TDEFieldIOTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        clean_output_files();
    }

    void TearDown() override
    {
        clean_output_files();
    }

    std::string output_path(const int field_index) const
    {
        return output_prefix_ + "efield_" + std::to_string(field_index) + ".txt";
    }

    void clean_output_files() const
    {
        for (int field_index = 0; field_index <= 3; ++field_index)
        {
            std::remove(output_path(field_index).c_str());
        }
        std::remove(blocking_path_.c_str());
    }

    std::shared_ptr<elecstate::TDFieldManager> create_manager(const std::vector<double>& amplitudes) const
    {
        Input_para input;
        input.td_vext = true;
        input.td_stype = 0;
        input.td_tstart = 0;
        input.td_tend = 10;
        input.td_dt = ModuleBase::AU_to_FS;
        input.td_ttype.clear();
        input.td_vext_dire.clear();
        input.td_heavi_t0.clear();
        input.td_heavi_amp.clear();
        for (std::size_t field_index = 0; field_index < amplitudes.size(); ++field_index)
        {
            input.td_ttype.push_back(3);
            input.td_vext_dire.push_back(static_cast<int>(field_index % 3) + 1);
            input.td_heavi_t0.push_back(10.0);
            input.td_heavi_amp.push_back(amplitudes[field_index]);
        }
        return elecstate::create_td_field_manager(input);
    }

    std::vector<std::pair<double, double>> read_samples(const int field_index) const
    {
        std::ifstream input(output_path(field_index).c_str());
        std::vector<std::pair<double, double>> samples;
        double time = 0.0;
        double field = 0.0;
        while (input >> time >> field)
        {
            samples.push_back(std::make_pair(time, field));
        }
        return samples;
    }

    const std::string output_prefix_ = "td_efield_io_test_";
    const std::string blocking_path_ = "td_efield_io_test_blocking_path";
};

TEST_F(TDEFieldIOTest, FreshCalculationTruncatesAndUsesOneBasedFiles)
{
    {
        std::ofstream output_1(output_path(1).c_str());
        std::ofstream output_2(output_path(2).c_str());
        output_1 << "old data\n";
        output_2 << "old data\n";
    }

    std::shared_ptr<elecstate::TDFieldManager> manager = create_manager({2.0, 3.0});
    ModuleIO::prepare_td_field_output(output_prefix_, manager->fields().size(), false);
    EXPECT_TRUE(read_samples(1).empty());
    EXPECT_TRUE(read_samples(2).empty());

    manager->advance_length_gauge();
    ModuleIO::write_td_field_values(*manager, output_prefix_);
    manager->advance_length_gauge();
    ModuleIO::write_td_field_values(*manager, output_prefix_);

    const std::vector<std::pair<double, double>> samples_1 = read_samples(1);
    const std::vector<std::pair<double, double>> samples_2 = read_samples(2);
    ASSERT_EQ(samples_1.size(), 2U);
    ASSERT_EQ(samples_2.size(), 2U);
    EXPECT_DOUBLE_EQ(samples_1[0].first, 0.0);
    EXPECT_DOUBLE_EQ(samples_2[0].first, 0.0);
    EXPECT_NEAR(samples_1[0].second, 2.0, 1.0e-12);
    EXPECT_NEAR(samples_2[0].second, 3.0, 1.0e-12);
    EXPECT_NEAR(samples_1[1].first, ModuleBase::AU_to_FS, 1.0e-7);
    EXPECT_NEAR(samples_2[1].first, ModuleBase::AU_to_FS, 1.0e-7);

    std::ifstream zero_based_file(output_path(0).c_str());
    EXPECT_FALSE(zero_based_file.good());
}

TEST_F(TDEFieldIOTest, RestartPreservesExistingSamples)
{
    {
        std::ofstream output(output_path(1).c_str());
        output << "7 8\n";
    }

    std::shared_ptr<elecstate::TDFieldManager> manager = create_manager({4.0});
    ModuleIO::prepare_td_field_output(output_prefix_, manager->fields().size(), true);
    manager->advance_length_gauge();
    ModuleIO::write_td_field_values(*manager, output_prefix_);

    const std::vector<std::pair<double, double>> samples = read_samples(1);
    ASSERT_EQ(samples.size(), 2U);
    EXPECT_DOUBLE_EQ(samples[0].first, 7.0);
    EXPECT_DOUBLE_EQ(samples[0].second, 8.0);
    EXPECT_DOUBLE_EQ(samples[1].first, 0.0);
    EXPECT_NEAR(samples[1].second, 4.0, 1.0e-12);
}

TEST_F(TDEFieldIOTest, ReportsOutputOpenFailures)
{
    {
        std::ofstream blocking_file(blocking_path_.c_str());
        blocking_file << "not a directory\n";
    }
    const std::string invalid_directory = blocking_path_ + "/";
    std::shared_ptr<elecstate::TDFieldManager> manager = create_manager({4.0});

    EXPECT_EXIT(ModuleIO::prepare_td_field_output(invalid_directory, manager->fields().size(), false), testing::ExitedWithCode(1), "");

    manager->advance_length_gauge();
    EXPECT_EXIT(ModuleIO::write_td_field_values(*manager, invalid_directory), testing::ExitedWithCode(1), "");
}

} // namespace
