#include "source_base/module_out/csr_reader.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <complex>

namespace
{
class CompactCsr : public testing::Test
{
protected:
    std::string path;
    void SetUp() override { path = testing::TempDir() + "/abacus_compact_csr_test.csr"; }
    void TearDown() override { std::remove(path.c_str()); }
    void write(const std::string& blocks, const std::string& dimension, const std::string& count)
    {
        std::ofstream out(path);
        out << "STEP: 7\nMatrix Dimension of Hexxs_0(R): " << dimension
            << "\nMatrix number of Hexxs_0(R): " << count << '\n' << blocks;
    }
};
TEST_F(CompactCsr, RealMultipleRAndWrappedValues)
{
    write("0 0 0 2\n 1.25\n -2\n 0 1\n 0 1 2\n0 0 1 1\n 3\n 1\n 0 1 1\n", "2", "2");
    ModuleIO::csrFileReader<double> reader(path);
    EXPECT_EQ(reader.getStep(), 7);
    EXPECT_EQ(reader.getMatrixDimension(), 2);
    EXPECT_EQ(reader.getNumberOfR(), 2);
    EXPECT_DOUBLE_EQ(reader.getMatrix(0)(0, 0), 1.25);
    EXPECT_DOUBLE_EQ(reader.getMatrix(0)(1, 1), -2);
    EXPECT_DOUBLE_EQ(reader.getMatrix(0, 0, 1)(0, 1), 3);
}
TEST_F(CompactCsr, ComplexValues)
{
    write("0 0 0 2\n (1,2) (-3,-4)\n0 1\n0 1 2\n", "2", "1");
    ModuleIO::csrFileReader<std::complex<double>> reader(path);
    EXPECT_EQ(reader.getMatrix(0)(0, 0), std::complex<double>(1, 2));
    EXPECT_EQ(reader.getMatrix(0)(1, 1), std::complex<double>(-3, -4));
}
TEST_F(CompactCsr, EmptyRankMatrix)
{
    write("", "2", "0");
    ModuleIO::csrFileReader<double> reader(path);
    EXPECT_EQ(reader.getMatrixDimension(), 2);
    EXPECT_EQ(reader.getNumberOfR(), 0);
}
TEST_F(CompactCsr, RejectsInvalidDimension)
{
    write("", "0", "0");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_F(CompactCsr, RejectsMalformedHeader)
{
    write("", "wrong", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_F(CompactCsr, RejectsTruncatedValues)
{
    write("0 0 0 2\n1\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_F(CompactCsr, RejectsInvalidColumn)
{
    write("0 0 0 1\n1\n2\n0 1 1\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_F(CompactCsr, RejectsInvalidRowPointers)
{
    write("0 0 0 1\n1\n0\n0 1 0\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_F(CompactCsr, RejectsMalformedNumber)
{
    write("0 0 0 1\nbad\n0\n0 1 1\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
}
