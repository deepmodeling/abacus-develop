#include "source_estate/scf_diagnostics.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <sstream>

TEST(ScfDiagnostics, PreservePrecisionMissingInputAndAppend)
{
    ScfDiagnostics::ScfDiagnosticStep step;
    step.output_dir = testing::TempDir();
    step.mixing_type = "broyden";
    step.threshold = 1e-6;
    step.charge_residual = 1.2345678901234567e-7;
    step.magnetic_residual = 2.5e-6;
    step.restart_count = 1;
    step.restart_executed = true;
    ScfDiagnostics::ScfOccupationBlock block;
    block.atom = 1;
    block.l = 2;
    block.raw_output.assign(25, 0.0);
    block.raw_output[0] = 0.12345678901234567;
    std::vector<ScfDiagnostics::ScfOccupationBlock> blocks(1, block);
    std::ostringstream log;
    ScfDiagnostics::write_scf_diagnostics(step, blocks, log);
    step.iteration = 2;
    blocks[0].input = block.raw_output;
    blocks[0].raw_before_occupation_mixing = {0.8};
    blocks[0].after_occupation_mixing = {0.6};
    blocks[0].after_dmr_mixing = {0.4};
    step.uramping_eV = 0.1;
    step.dmr_mixing_executed = true;
    ScfDiagnostics::write_scf_diagnostics(step, blocks, log);
    const std::string path = step.output_dir + "scf_diagnostics.jsonl";
    std::ifstream in(path.c_str());
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(text.find("1.234567890123456"), std::string::npos);
    EXPECT_NE(text.find("\"n_in\":null"), std::string::npos);
    EXPECT_NE(text.find("\"iteration\":2"), std::string::npos);
    EXPECT_NE(text.find("\"fundamental_gap_eV\":null"), std::string::npos);
    EXPECT_NE(text.find("\"charge_residual_L1_per_electron\":1.234567890123456"), std::string::npos);
    EXPECT_NE(text.find("\"mixing_restart_count\":1"), std::string::npos);
    EXPECT_NE(text.find("\"mixing_restart_executed\":1"), std::string::npos);
    EXPECT_EQ(std::count(text.begin(), text.end(), '\n'), 2);
    EXPECT_EQ(blocks[0].raw_output, block.raw_output);
    EXPECT_NE(text.find("\"schema\":4"), std::string::npos);
    EXPECT_NE(text.find("\"occupation_mixing_source_iteration\":1"), std::string::npos);
    EXPECT_NE(text.find("\"n_after_dmr_mixing\":[4.000000000000000"), std::string::npos);
    EXPECT_NE(text.find("\"n_after_occupation_mixing\":[5.999999999999999"), std::string::npos);
    EXPECT_NE(log.str().find(path), std::string::npos);
    in.close();
    std::remove(path.c_str());
}

TEST(ScfDiagnostics, NonWriterRankDoesNotOpenAFile)
{
    ScfDiagnostics::ScfDiagnosticStep step;
    step.rank = 1;
    step.output_dir = "/nonexistent/";
    std::ostringstream log;
    const std::vector<ScfDiagnostics::ScfOccupationBlock> blocks;
    EXPECT_NO_THROW(ScfDiagnostics::write_scf_diagnostics(step, blocks, log));
    EXPECT_TRUE(log.str().empty());
}
