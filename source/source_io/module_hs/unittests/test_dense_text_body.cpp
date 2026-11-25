#include "../dense_text_body.h"

#include <complex>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

namespace
{
/// @brief Read a written matrix body back and tokenize every number in it.
/// @details Parentheses and commas of the complex layout are turned into
///          whitespace so that a single tokenizer serves both overloads; the
///          value order is preserved (real, imaginary, real, ...).
std::vector<double> read_numbers(const std::string& fname)
{
    std::ifstream ifs(fname.c_str());
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    std::string text = buffer.str();
    for (std::size_t i = 0; i < text.size(); ++i)
    {
        const char c = text[i];
        if (c == '(' || c == ')' || c == ',')
        {
            text[i] = ' ';
        }
    }
    std::istringstream iss(text);
    std::vector<double> values;
    double value = 0.0;
    while (iss >> value)
    {
        values.push_back(value);
    }
    return values;
}

std::string read_text(const std::string& fname)
{
    std::ifstream ifs(fname.c_str());
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    return buffer.str();
}

std::string body_path(const std::string& name)
{
    const std::string temp_dir = testing::TempDir();
    return temp_dir + name;
}
} // namespace

// The file body is row-major even when the gathered buffer is column-major,
// which is what the column-major eigensolvers leave behind. A non-symmetric
// probe makes an accidental transpose visible.
TEST(DenseTextBody, DoubleColumnMajorSourceWritesRowMajorBody)
{
    const int dim = 3;
    // Column-major buffer of the row-major matrix 1..9.
    const std::vector<double> mat = {1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0};
    const std::string fname = body_path("dense_text_body_double_cm.txt");

    std::ofstream ofs(fname.c_str());
    ASSERT_TRUE(ofs.is_open());
    ModuleIO::write_dense_text_body(ofs, mat, dim, 8, true);
    ofs.close();

    const std::vector<double> values = read_numbers(fname);
    ASSERT_EQ(values.size(), 9u);
    for (int i = 0; i < dim * dim; ++i)
    {
        EXPECT_NEAR(values[i], static_cast<double>(i + 1), 1e-12);
    }
    EXPECT_EQ(std::remove(fname.c_str()), 0);
}

TEST(DenseTextBody, DoubleRowMajorSourceIsUnchanged)
{
    const int dim = 2;
    const std::vector<double> mat = {1.0, 2.0, 3.0, 4.0};
    const std::string fname = body_path("dense_text_body_double_rm.txt");

    std::ofstream ofs(fname.c_str());
    ASSERT_TRUE(ofs.is_open());
    ModuleIO::write_dense_text_body(ofs, mat, dim, 8, false);
    ofs.close();

    const std::vector<double> values = read_numbers(fname);
    ASSERT_EQ(values.size(), 4u);
    EXPECT_NEAR(values[0], 1.0, 1e-12);
    EXPECT_NEAR(values[1], 2.0, 1e-12);
    EXPECT_NEAR(values[2], 3.0, 1e-12);
    EXPECT_NEAR(values[3], 4.0, 1e-12);
    EXPECT_EQ(std::remove(fname.c_str()), 0);
}

TEST(DenseTextBody, ComplexColumnMajorSourceWritesRowMajorBody)
{
    const int dim = 2;
    // Row-major probe: (0,0)=1, (0,1)=2+3i, (1,0)=4+5i, (1,1)=6.
    // The column-major buffer swaps the two off-diagonal elements, so an
    // accidental transpose would emit 4+5i where 2+3i belongs.
    const std::vector<std::complex<double>> mat
        = {std::complex<double>(1.0, 0.0), std::complex<double>(4.0, 5.0),
           std::complex<double>(2.0, 3.0), std::complex<double>(6.0, 0.0)};
    const std::string fname = body_path("dense_text_body_complex_cm.txt");

    std::ofstream ofs(fname.c_str());
    ASSERT_TRUE(ofs.is_open());
    ModuleIO::write_dense_text_body(ofs, mat, dim, 8, true);
    ofs.close();

    const std::vector<double> values = read_numbers(fname);
    ASSERT_EQ(values.size(), 8u);
    EXPECT_NEAR(values[0], 1.0, 1e-12);
    EXPECT_NEAR(values[1], 0.0, 1e-12);
    EXPECT_NEAR(values[2], 2.0, 1e-12);
    EXPECT_NEAR(values[3], 3.0, 1e-12);
    EXPECT_NEAR(values[4], 4.0, 1e-12);
    EXPECT_NEAR(values[5], 5.0, 1e-12);
    EXPECT_NEAR(values[6], 6.0, 1e-12);
    EXPECT_NEAR(values[7], 0.0, 1e-12);
    EXPECT_EQ(std::remove(fname.c_str()), 0);
}

// Components below the internal threshold are flushed to zero, and the record
// ends with a single newline so that appended records start on a fresh line.
TEST(DenseTextBody, ComplexNoiseIsFlushedAndRecordIsTerminated)
{
    const int dim = 2;
    const std::vector<std::complex<double>> mat
        = {std::complex<double>(1.0, 1.0e-18), std::complex<double>(1.0e-20, 0.0),
           std::complex<double>(0.0, 0.0), std::complex<double>(2.0, -3.0)};
    const std::string fname = body_path("dense_text_body_complex_thr.txt");

    std::ofstream ofs(fname.c_str());
    ASSERT_TRUE(ofs.is_open());
    ModuleIO::write_dense_text_body(ofs, mat, dim, 8, false);
    ofs.close();

    const std::vector<double> values = read_numbers(fname);
    ASSERT_EQ(values.size(), 8u);
    EXPECT_NEAR(values[1], 0.0, 1e-30);
    EXPECT_NEAR(values[2], 0.0, 1e-30);
    EXPECT_NEAR(values[7], -3.0, 1e-12);

    const std::string text = read_text(fname);
    ASSERT_FALSE(text.empty());
    EXPECT_EQ(text[text.size() - 1], '\n');
    EXPECT_EQ(std::remove(fname.c_str()), 0);
}
