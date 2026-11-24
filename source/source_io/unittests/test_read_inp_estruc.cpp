#include "source_io/module_parameter/read_input.h"

#include <gtest/gtest.h>
#include <cstdio>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace
{
class PPCGInputTest : public testing::Test
{
protected:
    std::string filename;

    void SetUp() override
    {
        const std::string pattern = testing::TempDir() + "ppcg-input-XXXXXX";
        std::vector<char> path(pattern.begin(), pattern.end());
        path.push_back('\0');
        const int descriptor = mkstemp(path.data());
        ASSERT_GE(descriptor, 0);
        close(descriptor);
        filename = path.data();
    }

    void TearDown() override
    {
        std::remove(filename.c_str());
    }

    void write_input(const std::string& solver, const std::string& calculation)
    {
        std::ofstream out(filename);
        out << "INPUT_PARAMETERS\nbasis_type pw\nntype 1\necutwfc 20\n"
            << "ks_solver " << solver << "\ncalculation " << calculation << '\n';
        ASSERT_TRUE(out.good());
    }
};

TEST_F(PPCGInputTest, RejectsNscfBeforeDiagonalization)
{
    write_input("ppcg", "nscf");
    ModuleIO::ReadInput reader(0);
    reader.check_ntype_flag = false;
    Parameter input;
    testing::internal::CaptureStdout();
    EXPECT_EXIT(reader.read_parameters(input, filename), testing::ExitedWithCode(1), "");
    const std::string output = testing::internal::GetCapturedStdout();
    EXPECT_NE(output.find("calculation!=nscf"), std::string::npos);
}

TEST_F(PPCGInputTest, AcceptsScfAndPreservesOtherNscfSolvers)
{
    for (const std::string solver : {"ppcg", "cg", "dav", "dav_subspace"})
    {
        const std::string calculation = solver == "ppcg" ? "scf" : "nscf";
        SCOPED_TRACE(solver);
        write_input(solver, calculation);
        ModuleIO::ReadInput reader(0);
        reader.check_ntype_flag = false;
        Parameter input;
        reader.read_parameters(input, filename);
        EXPECT_EQ(input.inp.ks_solver, solver);
        EXPECT_EQ(input.inp.calculation, calculation);
    }
}
} // namespace
