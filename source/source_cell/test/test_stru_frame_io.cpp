/**
 * @file test_stru_frame_io.cpp
 * @brief Unit test for stru_frame_io: last-frame extraction from a fixed-name
 * cumulative STRU file written in out_app_flag == 1 mode.
 */
#include "source_cell/stru_frame_io.h"

#include <fstream>
#include <sstream>
#include <string>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_base/fs_compat.h"

class StruFrameIoTest : public testing::Test
{
  protected:
    std::string case_dir;

    void SetUp() override
    {
        case_dir = testing::TempDir() + "/abacus_stru_frame_io/";
        ModuleBase::make_directory(case_dir);
    }

    std::string read_all(const std::string& path)
    {
        std::ifstream ifs(path.c_str());
        std::stringstream content;
        content << ifs.rdbuf();
        return content.str();
    }
};

TEST_F(StruFrameIoTest, ExtractLastFramePicksNewest)
{
    const std::string cumulative = case_dir + "STRU";
    std::ofstream ofs(cumulative.c_str());
    ofs << "# MD step 1\n\nATOMIC_SPECIES\nSi 28.0855\n";
    ofs << "# MD step 2\n\nATOMIC_SPECIES\nSi 28.0855\n\nLATTICE_CONSTANT\n10.2\n";
    ofs.close();

    const std::string frame = case_dir + "STRU_restart_frame";
    EXPECT_TRUE(stru_frame_io::extract_last_frame(cumulative, frame));

    const std::string content = read_all(frame);
    // the extracted file must hold only the last frame, not the earlier one
    EXPECT_EQ(content.find("# MD step"), std::string::npos);
    EXPECT_NE(content.find("LATTICE_CONSTANT"), std::string::npos);
}

TEST_F(StruFrameIoTest, ExtractLastFrameRewritesOnSecondCall)
{
    const std::string cumulative = case_dir + "STRU";
    std::ofstream ofs(cumulative.c_str());
    ofs << "# MD step 1\n\nATOMIC_SPECIES\nSi 28.0855\n";
    ofs.close();

    const std::string frame = case_dir + "STRU_restart_frame";
    EXPECT_TRUE(stru_frame_io::extract_last_frame(cumulative, frame));

    // simulate a later restart: append another frame and re-extract
    std::ofstream app(cumulative.c_str(), std::ios::app);
    app << "# MD step 2\n\nLATTICE_CONSTANT\n10.9\n";
    app.close();
    EXPECT_TRUE(stru_frame_io::extract_last_frame(cumulative, frame));

    const std::string content = read_all(frame);
    EXPECT_NE(content.find("10.9"), std::string::npos);
}

TEST_F(StruFrameIoTest, ExtractLastFrameFailsOnMissingFile)
{
    const std::string missing = case_dir + "NO_SUCH_STRU";
    const std::string frame = case_dir + "STRU_restart_frame";
    EXPECT_FALSE(stru_frame_io::extract_last_frame(missing, frame));
}

TEST_F(StruFrameIoTest, ExtractLastFrameFailsWithoutMarker)
{
    const std::string cumulative = case_dir + "STRU";
    std::ofstream ofs(cumulative.c_str());
    ofs << "ATOMIC_SPECIES\nSi 28.0855\n";
    ofs.close();

    const std::string frame = case_dir + "STRU_restart_frame";
    EXPECT_FALSE(stru_frame_io::extract_last_frame(cumulative, frame));
}
