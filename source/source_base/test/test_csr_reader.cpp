#include "source_base/module_out/csr_reader.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <complex>
#include <iterator>

namespace
{
class CsrFormats : public testing::TestWithParam<bool>
{
protected:
    std::string path;
    void SetUp() override { path = testing::TempDir() + "/abacus_compact_csr_test.csr"; }
    void TearDown() override { std::remove(path.c_str()); }
    void write(const std::string& blocks, const std::string& dimension, const std::string& count)
    {
        std::ofstream out(path);
        if (GetParam())
        {
            out << "STEP: 7\nMatrix Dimension of Hexxs_0(R): " << dimension
                << "\nMatrix number of Hexxs_0(R): " << count << '\n';
        }
        else
        {
            out << " --- Ionic Step 7 ---\n # matrix\n1 # spins\n1 # spin\n"
                << dimension << " # basis\n" << count << " # R blocks\n\n"
                << "user_defined_lattice\n1\n1 0 0\n0 1 0\n0 0 1\nSi\n1\nDirect\n0 0 0\n";
        }
        out << blocks;
    }
};
TEST_P(CsrFormats, RealMultipleRAndWrappedValues)
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
TEST_P(CsrFormats, ComplexValues)
{
    write("0 0 0 2\n (1,2) (-3,-4)\n0 1\n0 1 2\n", "2", "1");
    ModuleIO::csrFileReader<std::complex<double>> reader(path);
    EXPECT_EQ(reader.getMatrix(0)(0, 0), std::complex<double>(1, 2));
    EXPECT_EQ(reader.getMatrix(0)(1, 1), std::complex<double>(-3, -4));
}
TEST_P(CsrFormats, EmptyRankMatrix)
{
    write("", "2", "0");
    ModuleIO::csrFileReader<double> reader(path);
    EXPECT_EQ(reader.getMatrixDimension(), 2);
    EXPECT_EQ(reader.getNumberOfR(), 0);
}
TEST_P(CsrFormats, RejectsInvalidDimension)
{
    write("", "0", "0");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_P(CsrFormats, RejectsMalformedHeader)
{
    write("", "wrong", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_P(CsrFormats, RejectsTruncatedValues)
{
    write("0 0 0 2\n1\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_P(CsrFormats, RejectsInvalidColumn)
{
    write("0 0 0 1\n1\n2\n0 1 1\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_P(CsrFormats, RejectsInvalidRowPointers)
{
    write("0 0 0 1\n1\n0\n0 1 0\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_P(CsrFormats, RejectsMalformedNumber)
{
    write("0 0 0 1\nbad\n0\n0 1 1\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}

TEST_P(CsrFormats, CommentsAndEmptyBlock)
{
    write("# blocks\n0 0 0 0\n# values\n\n# columns\n\n# pointers\n0 0 0\n"
          "\n0 0 1 1\n# values\n2.5\n# columns\n1\n# pointers\n0 1 1\n", "2", "2");
    ModuleIO::csrFileReader<double> reader(path);
    EXPECT_TRUE(reader.getMatrix(0).getElements().empty());
    EXPECT_DOUBLE_EQ(reader.getMatrix(1)(0, 1), 2.5);
}
TEST_P(CsrFormats, RejectsDuplicateR)
{
    write("0 0 0 1\n1\n0\n0 1 1\n0 0 0 1\n2\n0\n0 1 1\n", "2", "2");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_P(CsrFormats, RejectsExtraValueShiftingRowPointers)
{
    write("0 0 0 1\n1 0\n0\n0 1 1\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_P(CsrFormats, AcceptsCommentAfterFinalRowPointer)
{
    write("0 0 0 1\n1\n0\n0 1 1 # end of block\n", "2", "1");
    ModuleIO::csrFileReader<double> reader(path);
    EXPECT_DOUBLE_EQ(reader.getMatrix(0)(0, 0), 1.0);
}
TEST_P(CsrFormats, RejectsWrappedLeftoverRowPointer)
{
    write("0 0 0 1\n1 0\n0\n0 1\n1\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
TEST_P(CsrFormats, PreservesNextFrameAfterComments)
{
    write("0 0 0 1\n1\n0\n0\n1 1\n\n# frame end\n", "2", "1");
    std::ifstream input(path);
    const std::istreambuf_iterator<char> begin(input);
    const std::istreambuf_iterator<char> end;
    const std::string frame(begin, end);
    input.close();
    std::ofstream append(path, std::ios::app);
    append << frame;
    append.close();
    ModuleIO::csrFileReader<double> reader(path);
    EXPECT_DOUBLE_EQ(reader.getMatrix(0)(0, 0), 1.0);
    reader.parseFile();
    EXPECT_EQ(reader.getStep(), 7);
    EXPECT_DOUBLE_EQ(reader.getMatrix(0)(0, 0), 1.0);
}
TEST_P(CsrFormats, RejectsDataAfterTrailingComment)
{
    write("0 0 0 1\n1\n0\n0 1 1\n# end\n\n9\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
INSTANTIATE_TEST_SUITE_P(BothHeaders, CsrFormats, testing::Bool());
}
