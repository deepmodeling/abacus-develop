#include "module_io/sparse_matrix.h"

#include <complex>
#include <fstream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

template <typename T>
T get_value(const T& val)
{
    return val;
}

template <typename T>
T get_value(const std::complex<T>& val)
{
    return val.real();
}

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

template <typename T>
class SparseMatrixTest : public ::testing::Test
{
  protected:
    ModuleIO::SparseMatrix<T> sm = ModuleIO::SparseMatrix<T>(4, 4);
};

using MyTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(SparseMatrixTest, MyTypes);

TYPED_TEST(SparseMatrixTest, Insert)
{
    // Add a value to the matrix with row and column indices
    this->sm.insert(2, 2, static_cast<TypeParam>(3.0));
    this->sm.insert(3, 3, static_cast<TypeParam>(4.0));

    EXPECT_EQ(this->sm.getNNZ(), 2);
    EXPECT_EQ(this->sm.getCols(), 4);
    EXPECT_EQ(this->sm.getRows(), 4);

    // Check if the value was added correctly
    EXPECT_DOUBLE_EQ(get_value(this->sm(2, 2)), 3.0);
    EXPECT_DOUBLE_EQ(get_value(this->sm(3, 3)), 4.0);
}

TYPED_TEST(SparseMatrixTest, InsertWarning)
{
    // case 1
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sm.insert(2, -1, static_cast<TypeParam>(3.0)), ::testing::ExitedWithCode(0), "");
    std::string output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 2
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sm.insert(-1, 0, static_cast<TypeParam>(3.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 3
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sm.insert(2, 4, static_cast<TypeParam>(3.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
    // case 4
    testing::internal::CaptureStdout();
    EXPECT_EXIT(this->sm.insert(4, 2, static_cast<TypeParam>(3.0)), ::testing::ExitedWithCode(0), "");
    output = testing::internal::GetCapturedStdout();
    // test output on screen
    EXPECT_THAT(output, testing::HasSubstr("row or col index out of range"));
}

TYPED_TEST(SparseMatrixTest, CSR)
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
    ModuleIO::SparseMatrix<TypeParam> sm0;
    ModuleIO::SparseMatrix<TypeParam> sm1;
    sm0.setRows(4);
    sm0.setCols(6);
    sm0.setSparseThreshold(0.6);
    sm0.insert(0, 0, static_cast<TypeParam>(1.0));
    sm0.insert(0, 1, static_cast<TypeParam>(2.0));
    sm0.insert(0, 2, static_cast<TypeParam>(0.5));
    sm0.insert(0, 3, static_cast<TypeParam>(0.4));
    sm0.insert(1, 1, static_cast<TypeParam>(3.0));
    sm0.insert(1, 3, static_cast<TypeParam>(4.0));
    sm0.insert(2, 2, static_cast<TypeParam>(5.0));
    sm0.insert(2, 3, static_cast<TypeParam>(6.0));
    sm0.insert(2, 4, static_cast<TypeParam>(7.0));
    sm0.insert(3, 5, static_cast<TypeParam>(8.0));
    EXPECT_DOUBLE_EQ(sm0.getSparseThreshold(), 0.6);
    testing::internal::CaptureStdout();
    // 2 is the precision
    sm0.printToCSR(std::cout, 2);
    std::string output = testing::internal::GetCapturedStdout();
    if (std::is_same<TypeParam, double>::value)
    {
        EXPECT_THAT(output,
                    testing::HasSubstr("1.00e+00 2.00e+00 3.00e+00 4.00e+00 5.00e+00 6.00e+00 7.00e+00 8.00e+00"));
        EXPECT_THAT(output, testing::HasSubstr("0 1 1 3 2 3 4 5"));
        EXPECT_THAT(output, testing::HasSubstr("0 2 4 7 8"));

        std::vector<TypeParam> expected_csr_values = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        std::vector<int> expected_csr_col_ind = {0, 1, 1, 3, 2, 3, 4, 5};
        std::vector<int> expected_csr_row_ptr = {0, 2, 4, 7, 8};

        sm1.setRows(4);
        sm1.setCols(6);
        sm1.readCSR(expected_csr_values, expected_csr_col_ind, expected_csr_row_ptr);

        EXPECT_EQ(sm0.getNNZ(), 8);
        EXPECT_EQ(sm1.getNNZ(), 8);

        for (const auto& element: sm0.getElements())
        {
            auto it = sm1.getElements().find(element.first);
            EXPECT_DOUBLE_EQ(get_value(it->second), get_value(element.second));
        }
    }
    else if (std::is_same<TypeParam, std::complex<double>>::value)
    {
        EXPECT_THAT(
            output,
            testing::HasSubstr("(1.00e+00,0.00e+00) (2.00e+00,0.00e+00) (3.00e+00,0.00e+00) (4.00e+00,0.00e+00) "
                               "(5.00e+00,0.00e+00) (6.00e+00,0.00e+00) (7.00e+00,0.00e+00) (8.00e+00,0.00e+00)"));
        EXPECT_THAT(output, testing::HasSubstr("0 1 1 3 2 3 4 5"));
        EXPECT_THAT(output, testing::HasSubstr("0 2 4 7 8"));

        std::vector<TypeParam> expected_csr_values = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
        std::vector<int> expected_csr_col_ind = {0, 1, 1, 3, 2, 3, 4, 5};
        std::vector<int> expected_csr_row_ptr = {0, 2, 4, 7, 8};

        sm1.setRows(4);
        sm1.setCols(6);
        sm1.readCSR(expected_csr_values, expected_csr_col_ind, expected_csr_row_ptr);

        for (const auto& element: sm0.getElements())
        {
            auto it = sm1.getElements().find(element.first);
            EXPECT_DOUBLE_EQ(get_value(it->second), get_value(element.second));
        }
        // index matrix elements
        EXPECT_DOUBLE_EQ(get_value(sm1(0, 1)), 2.0);
        EXPECT_DOUBLE_EQ(get_value(sm1(1, 1)), 3.0);
        EXPECT_DOUBLE_EQ(get_value(sm1(1, 2)), 0.0);
    }
}
