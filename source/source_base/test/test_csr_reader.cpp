#include "source_base/module_out/csr_reader.h"
#include "gtest/gtest.h"
#include <cstdio>
#include <fstream>
#include <complex>
#include <iterator>

namespace
{
// Header spellings covered by the parameter: 0 = legacy compact "STEP:",
// 1 = current compact "Ionic_Step", 2 = legacy cell-aware
// "--- Ionic Step ---", 3 = current cell-aware "Ionic_Step".
class CsrFormats : public testing::TestWithParam<int>
{
protected:
    std::string path;
    void SetUp() override { path = testing::TempDir() + "/abacus_compact_csr_test.csr"; }
    void TearDown() override { std::remove(path.c_str()); }
    void write(const std::string& blocks, const std::string& dimension, const std::string& count)
    {
        std::ofstream out(path);
        const int format = GetParam();
        if (format == 0)
        {
            out << "STEP: 7\nMatrix Dimension of Hexxs_0(R): " << dimension
                << "\nMatrix number of Hexxs_0(R): " << count << '\n';
        }
        else if (format == 1)
        {
            out << "Ionic_Step 7\nMatrix Dimension of Hexxs_0(R): " << dimension
                << "\nMatrix number of Hexxs_0(R): " << count << '\n';
        }
        else if (format == 2)
        {
            out << " --- Ionic Step 7 ---\n # matrix\n1 # spins\n1 # spin\n"
                << dimension << " # basis\n" << count << " # R blocks\n\n"
                << "user_defined_lattice\n1\n1 0 0\n0 1 0\n0 0 1\nSi\n1\nDirect\n0 0 0\n";
        }
        else
        {
            out << "Ionic_Step 7\n # matrix\n1 # spins\n1 # spin\n"
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
    std::string frame(begin, end);
    input.close();
    const std::size_t step_position = frame.find('7');
    ASSERT_NE(step_position, std::string::npos);
    frame.replace(step_position, 1, "8");
    const std::string first_block = "0 0 0 1\n1\n";
    const std::string second_block = "0 0 0 1\n2.5\n";
    const std::size_t block_position = frame.find(first_block);
    const std::size_t block_size = first_block.size();
    ASSERT_NE(block_position, std::string::npos);
    frame.replace(block_position, block_size, second_block);
    std::ofstream append(path, std::ios::app);
    append << frame;
    append.close();
    ModuleIO::csrFileReader<double> reader(path);
    EXPECT_DOUBLE_EQ(reader.getMatrix(0)(0, 0), 1.0);
    EXPECT_EQ(reader.getStep(), 7);
    reader.parseFile();
    EXPECT_EQ(reader.getStep(), 8);
    EXPECT_DOUBLE_EQ(reader.getMatrix(0)(0, 0), 2.5);
}
TEST_P(CsrFormats, RejectsDataAfterTrailingComment)
{
    write("0 0 0 1\n1\n0\n0 1 1\n# end\n\n9\n", "2", "1");
    EXPECT_EXIT({ModuleIO::csrFileReader<double> reader(path);}, testing::ExitedWithCode(1), "");
}
INSTANTIATE_TEST_SUITE_P(AllHeaders, CsrFormats, testing::Values(0, 1, 2, 3));
}
