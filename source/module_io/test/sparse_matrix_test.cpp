#include "module_io/sparse_matrix.h"

#include <complex>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

/************************************************
 *  unit test of sparse_matrix.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *  - insert()
 *   - Add value to the matrix with row and column indices
 * - printToCSR()
 *   - Print data in CSR format
 * - readCSR()
 *   - Read CSR data from arrays
 */

class SparseMatrixTest : public testing::Test
{
  protected:
    ModuleIO::SparseMatrix<double> smd = ModuleIO::SparseMatrix<double>(4, 4);
    ModuleIO::SparseMatrix<std::complex<double>> smc = ModuleIO::SparseMatrix<std::complex<double>>(4, 4);
    std::string output;
};

TEST_F(SparseMatrixTest, InsertDouble)
{
    // Add a value to the matrix with row and column indices
    smd.insert(2, 2, 3.0);
    smd.insert(3, 3, 4.0);

    EXPECT_EQ(smd.getNNZ(), 2);
    EXPECT_EQ(smd.getCols(), 4);
    EXPECT_EQ(smd.getRows(), 4);

    // Check if the value was added correctly
    EXPECT_DOUBLE_EQ(smd(2, 2), 3.0);
    EXPECT_DOUBLE_EQ(smd(3, 3), 4.0);
}

TEST_F(SparseMatrixTest, InsertDoubleWarning)
{
    // case 1
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smd.insert(2, -1, 3.0), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 2
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smd.insert(-1, 0, 3.0), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 3
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smd.insert(2, 4, 3.0), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 4
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smd.insert(4, 2, 3.0), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
}

TEST_F(SparseMatrixTest, InsertComplex)
{
    // Add a value to the matrix with row and column indices
    smc.insert(2, 2, std::complex<double>(1.0, 0.0));
    smc.insert(3, 3, std::complex<double>(0.0, 0.2));

    EXPECT_EQ(smc.getNNZ(), 2);
    EXPECT_EQ(smc.getCols(), 4);
    EXPECT_EQ(smc.getRows(), 4);

    // Check if the value was added correctly
    EXPECT_DOUBLE_EQ(smc(2, 2).real(), 1.0);
    EXPECT_DOUBLE_EQ(smc(2, 2).imag(), 0.0);
    EXPECT_DOUBLE_EQ(smc(3, 3).real(), 0.0);
    EXPECT_DOUBLE_EQ(smc(3, 3).imag(), 0.2);
}

TEST_F(SparseMatrixTest, InsertComplexWarning)
{
    // case 1
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smc.insert(2, -1, std::complex<double>(1.0, 0.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 2
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smc.insert(-1, 0, std::complex<double>(1.0, 0.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 3
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smc.insert(2, 4, std::complex<double>(1.0, 0.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 4
    testing::internal::CaptureStdout();
    EXPECT_EXIT(smc.insert(4, 2, std::complex<double>(1.0, 0.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
}

TEST_F(SparseMatrixTest, CSRDouble)
{
    // test the printCSR function for the following matrix
    // 1.0 2.0 0.5 0.4 0.0 0.0
    // 0.0 3.0 0.0 4.0 0.0 0.0
    // 0.0 0.0 5.0 6.0 7.0 0.0
    // 0.0 0.0 0.0 0.0 0.0 8.0
    // the expect output is
    // csr_values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    // csr_col_ind = [0, 1, 1, 3, 2, 3, 4, 5]
    // csr_row_ptr = [0, 2, 4, 7, 8]
    // set the sparse threshold to 0.6 and the precision to 2
    ModuleIO::SparseMatrix<double> smd0;
    ModuleIO::SparseMatrix<double> smd1;
    smd0.setRows(4);
    smd0.setCols(6);
    smd0.setSparseThreshold(0.6);
    smd0.insert(0, 0, 1.0);
    smd0.insert(0, 1, 2.0);
    smd0.insert(0, 2, 0.5);
    smd0.insert(0, 3, 0.4);
    smd0.insert(1, 1, 3.0);
    smd0.insert(1, 3, 4.0);
    smd0.insert(2, 2, 5.0);
    smd0.insert(2, 3, 6.0);
    smd0.insert(2, 4, 7.0);
    smd0.insert(3, 5, 8.0);
    EXPECT_DOUBLE_EQ(smd0.getSparseThreshold(), 0.6);
    testing::internal::CaptureStdout();
    smd0.printToCSR(std::cout, 2);
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

    EXPECT_EQ(smd0.getNNZ(), 8);
    EXPECT_EQ(smd1.getNNZ(), 8);

    for (const auto& element : smd0.getElements())
    {
        auto it = smd1.getElements().find(element.first);
        EXPECT_DOUBLE_EQ(it->second, element.second);
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
    // set the sparse threshold to 0.6 and the precision to 1
    ModuleIO::SparseMatrix<std::complex<double>> smc0;
    ModuleIO::SparseMatrix<std::complex<double>> smc1;
    smc0.setRows(4);
    smc0.setCols(6);
    smc0.setSparseThreshold(0.6);
    EXPECT_DOUBLE_EQ(smc0.getSparseThreshold(), 0.6);
    smc0.insert(0, 0, std::complex<double>(1.0, 0.0));
    smc0.insert(0, 1, std::complex<double>(2.0, 0.0));
    smc0.insert(0, 2, std::complex<double>(0.5, 0.0));
    smc0.insert(0, 3, std::complex<double>(0.0, 0.4));
    smc0.insert(1, 1, std::complex<double>(3.0, 0.0));
    smc0.insert(1, 3, std::complex<double>(4.0, 0.0));
    smc0.insert(2, 2, std::complex<double>(5.0, 0.0));
    smc0.insert(2, 3, std::complex<double>(6.0, 0.0));
    smc0.insert(2, 4, std::complex<double>(7.0, 0.0));
    smc0.insert(3, 5, std::complex<double>(8.0, 0.0));
    testing::internal::CaptureStdout();
    smc0.printToCSR(std::cout, 1);
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

    for (const auto& element : smc0.getElements())
    {
        auto it = smc1.getElements().find(element.first);
        EXPECT_DOUBLE_EQ(it->second.real(), element.second.real());
        EXPECT_DOUBLE_EQ(it->second.imag(), element.second.imag());
    }
    // index matrix elements
    EXPECT_DOUBLE_EQ(smc1(0, 1).imag(), 0.0);
    EXPECT_DOUBLE_EQ(smc1(0, 1).real(), 2.0);
    EXPECT_DOUBLE_EQ(smc1(1, 1).imag(), 0.0);
    EXPECT_DOUBLE_EQ(smc1(1, 1).real(), 3.0);
    EXPECT_DOUBLE_EQ(smc1(1, 2).imag(), 0.0);
    EXPECT_DOUBLE_EQ(smc1(1, 2).real(), 0.0);
}

#include "mpi.h"
int main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        testing::InitGoogleTest(&argc, argv);

        int nproc;
        int my_rank;

        MPI_Comm_size(MPI_COMM_WORLD,&nproc);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

        int result = RUN_ALL_TESTS();

        MPI_Finalize();

        return result;
}