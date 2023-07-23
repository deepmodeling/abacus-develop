#include "module_io/sparse_matrix.h"

#include <complex>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of sparse_matrix.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *  - addValue()
 *   - Add value to the matrix with row and column indices
 * - convertToCSR()
 *  - Convert to CSR format
 * - printCSR()
 * - Print data in CSR format
 * - readCSR()
 * - Read CSR data from arrays
 */

class SparseMatrixTest : public testing::Test
{
  protected:
    ModuleIO::SparseMatrix<std::complex<double>> smc0;
    ModuleIO::SparseMatrix<double> smd = ModuleIO::SparseMatrix<double>(4, 4);
    ModuleIO::SparseMatrix<std::complex<double>> smc = ModuleIO::SparseMatrix<std::complex<double>>(4, 4);
    std::string output;
};

TEST_F(SparseMatrixTest, AddValueDouble)
{
    // Add a value to the matrix with row and column indices
    smd.addValue(2, 2, 3.0);
    smd.addValue(3, 3, 4.0);

    EXPECT_EQ(smd.getData().size(), 2);
    EXPECT_EQ(smd.getCols(), 4);
    EXPECT_EQ(smd.getRows(), 4);

    // Check if the value was added correctly
    EXPECT_EQ(std::make_tuple(2, 2, 3.0), smd.getData()[0]);
    EXPECT_EQ(std::make_tuple(3, 3, 4.0), smd.getData()[1]);
}

TEST_F(SparseMatrixTest, AddValueDoubleWarning)
{
    // case 1
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smd.addValue(2, -1, 3.0), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 2
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smd.addValue(-1, 0, 3.0), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 3
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smd.addValue(2, 4, 3.0), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 4
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smd.addValue(4, 2, 3.0), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
}

TEST_F(SparseMatrixTest, AddValueComplex)
{
    // Add a value to the matrix with row and column indices
    smc.addValue(2, 2, std::complex<double>(1.0, 0.0));
    smc.addValue(3, 3, std::complex<double>(0.0, 0.2));

    EXPECT_EQ(smc.getData().size(), 2);
    EXPECT_EQ(smc.getCols(), 4);
    EXPECT_EQ(smc.getRows(), 4);

    // Check if the value was added correctly
    EXPECT_EQ(std::make_tuple(2, 2, std::complex<double>(1.0, 0.0)), smc.getData()[0]);
    EXPECT_EQ(std::make_tuple(3, 3, std::complex<double>(0.0, 0.2)), smc.getData()[1]);
}

TEST_F(SparseMatrixTest, AddValueComplexWarning)
{
    // case 1
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smc.addValue(2, -1, std::complex<double>(1.0, 0.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 2
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smc.addValue(-1, 0, std::complex<double>(1.0, 0.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 3
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smc.addValue(2, 4, std::complex<double>(1.0, 0.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 4
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smc.addValue(4, 2, std::complex<double>(1.0, 0.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
}

TEST_F(SparseMatrixTest, CSRDouble)
{
    // test the printCSR function for the following matrix
    // 1.0 2.0 0.0 0.0 0.0 0.0
    // 0.0 3.0 0.0 4.0 0.0 0.0
    // 0.0 0.0 5.0 6.0 7.0 0.0
    // 0.0 0.0 0.0 0.0 0.0 8.0
    // the expect output is
    // csr_values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    // csr_col_ind = [0, 1, 1, 3, 2, 3, 4, 5]
    // csr_row_ptr = [0, 2, 4, 7, 8]
    ModuleIO::SparseMatrix<double> smd0;
    ModuleIO::SparseMatrix<double> smd1;
    smd0.setRows(4);
    smd0.setCols(6);
    smd0.addValue(0, 0, 1.0);
    smd0.addValue(0, 1, 2.0);
    smd0.addValue(1, 1, 3.0);
    smd0.addValue(1, 3, 4.0);
    smd0.addValue(2, 2, 5.0);
    smd0.addValue(2, 3, 6.0);
    smd0.addValue(2, 4, 7.0);
    smd0.addValue(3, 5, 8.0);
    smd0.convertToCSR(0.0);
    testing::internal::CaptureStdout();
    smd0.printCSR(std::cout, 2);
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("1.00e+00 2.00e+00 3.00e+00 4.00e+00 5.00e+00 6.00e+00 7.00e+00 8.00e+00"));
    EXPECT_THAT(output, testing::HasSubstr("0 1 1 3 2 3 4 5"));
    EXPECT_THAT(output, testing::HasSubstr("0 2 4 7 8"));

    std::vector<double> expected_csr_values = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    std::vector<int> expected_csr_col_ind = {0, 1, 1, 3, 2, 3, 4, 5};
    std::vector<int> expected_csr_row_ptr = {0, 2, 4, 7, 8};

    smd1.setRows(4);
    smd1.setCols(6);
    smd1.readCSR(expected_csr_values, expected_csr_col_ind, expected_csr_row_ptr);

    for (int i = 0; i < 8; i++)
    {
        EXPECT_DOUBLE_EQ(smd0.getCSRValues()[i], smd1.getCSRValues()[i]);
        EXPECT_EQ(smd0.getCSRColInd()[i], smd1.getCSRColInd()[i]);
    }
    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(smd0.getCSRRowPtr()[i], smd1.getCSRRowPtr()[i]);
    }
}

TEST_F(SparseMatrixTest, CSRComplex)
{
    // test the printCSR function for the following matrix
    // 1.0 2.0 0.0 0.0 0.0 0.0
    // 0.0 3.0 0.0 4.0 0.0 0.0
    // 0.0 0.0 5.0 6.0 7.0 0.0
    // 0.0 0.0 0.0 0.0 0.0 8.0
    // the expect output is
    // csr_values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    // csr_col_ind = [0, 1, 1, 3, 2, 3, 4, 5]
    // csr_row_ptr = [0, 2, 4, 7, 8]
    ModuleIO::SparseMatrix<std::complex<double>> smc0;
    ModuleIO::SparseMatrix<std::complex<double>> smc1;
    smc0.setRows(4);
    smc0.setCols(6);
    smc0.addValue(0, 0, std::complex<double>(1.0, 0.0));
    smc0.addValue(0, 1, std::complex<double>(2.0, 0.0));
    smc0.addValue(1, 1, std::complex<double>(3.0, 0.0));
    smc0.addValue(1, 3, std::complex<double>(4.0, 0.0));
    smc0.addValue(2, 2, std::complex<double>(5.0, 0.0));
    smc0.addValue(2, 3, std::complex<double>(6.0, 0.0));
    smc0.addValue(2, 4, std::complex<double>(7.0, 0.0));
    smc0.addValue(3, 5, std::complex<double>(8.0, 0.0));
    smc0.convertToCSR(0.5);
    testing::internal::CaptureStdout();
    smc0.printCSR(std::cout, 1);
    output = testing::internal::GetCapturedStdout();
    EXPECT_THAT(output, testing::HasSubstr("(1.0e+00,0.0e+00) (2.0e+00,0.0e+00) (3.0e+00,0.0e+00) (4.0e+00,0.0e+00) (5.0e+00,0.0e+00) (6.0e+00,0.0e+00) (7.0e+00,0.0e+00) (8.0e+00,0.0e+00)"));
    EXPECT_THAT(output, testing::HasSubstr("0 1 1 3 2 3 4 5"));
    EXPECT_THAT(output, testing::HasSubstr("0 2 4 7 8"));

    std::vector<std::complex<double>> expected_csr_values = {std::complex<double>(1.0, 0.0),
                                                             std::complex<double>(2.0, 0.0),
                                                             std::complex<double>(3.0, 0.0),
                                                             std::complex<double>(4.0, 0.0),
                                                             std::complex<double>(5.0, 0.0),
                                                             std::complex<double>(6.0, 0.0),
                                                             std::complex<double>(7.0, 0.0),
                                                             std::complex<double>(8.0, 0.0)};
    std::vector<int> expected_csr_col_ind = {0, 1, 1, 3, 2, 3, 4, 5};
    std::vector<int> expected_csr_row_ptr = {0, 2, 4, 7, 8};

    smc1.setRows(4);
    smc1.setCols(6);
    smc1.readCSR(expected_csr_values, expected_csr_col_ind, expected_csr_row_ptr);

    for (int i = 0; i < 8; i++)
    {
        EXPECT_DOUBLE_EQ(smc0.getCSRValues()[i].real(), smc1.getCSRValues()[i].real());
        EXPECT_DOUBLE_EQ(smc0.getCSRValues()[i].imag(), smc1.getCSRValues()[i].imag());
        EXPECT_EQ(smc0.getCSRColInd()[i], smc1.getCSRColInd()[i]);
    }
    for (int i = 0; i < 5; i++)
    {
        EXPECT_EQ(smc0.getCSRRowPtr()[i], smc1.getCSRRowPtr()[i]);
    }
}