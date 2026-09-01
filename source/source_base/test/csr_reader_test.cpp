#include "source_base/module_out/csr_reader.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <complex>

/************************************************
 *  unit test of csr_reader.cpp
 ***********************************************/

/**
 * - Tested Functions:
 * - csrFileReader()
 *   - Constructor
 * - parseFile()
 *   - Read the file and parse it
 * - getNumberOfR()
 *   - Get number of R
 * - getMatrix(rx, ry, rz)
 *   - Get sparse matrix of a specific R coordinate
 * - getMatrix(index)
 *   - Get matrix by using index
 * - getRCoordinate()
 *   - Get R coordinate using index
 * - getStep()
 *   - Get step
 * - getMatrixDimension()
 *   - Get matrix dimension
 */

class csrFileReaderTest : public testing::Test
{
  protected:
    std::string filename = "./support/SR.csr";
};

TEST_F(csrFileReaderTest, CsrReader)
{
    ModuleIO::csrFileReader<double> csr(filename);
    // Check if file is open
    EXPECT_TRUE(csr.isOpen());
    // get step
    EXPECT_EQ(csr.getStep(), 1);
    // get matrix dimension
    EXPECT_EQ(csr.getMatrixDimension(), 4);
    // get number of R
    EXPECT_EQ(csr.getNumberOfR(), 2);
    // get R coordinate using index
    std::vector<int> RCoord;
    // get R coordinate using index
    RCoord = csr.getRCoordinate(0);
    EXPECT_EQ(RCoord[0], 0);
    EXPECT_EQ(RCoord[1], 1);
    EXPECT_EQ(RCoord[2], 1);
    RCoord = csr.getRCoordinate(1);
    EXPECT_EQ(RCoord[0], 0);
    EXPECT_EQ(RCoord[1], 0);
    EXPECT_EQ(RCoord[2], 0);

    // get matrix by using index
    ModuleIO::SparseMatrix<double> sparse_matrix;
    ModuleIO::SparseMatrix<double> sparse_matrix1;
    // the first matrix should be
    // 0 0 0 4
    // 0 0 7 0
    // 0 0 0 0
    // 0 0 0 0
    // the second matrix should be
    // 0 0 0 0
    // 0 0 0 0
    // 0 0 5 6
    // 0 0 0 10
    sparse_matrix = csr.getMatrix(0);
    sparse_matrix1 = csr.getMatrix(0, 1, 1);
    for (const auto& element: sparse_matrix.getElements())
    {
        auto it = sparse_matrix1.getElements().find(element.first);
        EXPECT_EQ(it->first.first, element.first.first);
        EXPECT_EQ(it->first.second, element.first.second);
        EXPECT_DOUBLE_EQ(it->second, element.second);
        // std::cout << "element( " << element.first.first << ", " << element.first.second << " ) = " << element.second
        // << std::endl;
    }
    EXPECT_DOUBLE_EQ(sparse_matrix(0, 3), 4.0);
    EXPECT_DOUBLE_EQ(sparse_matrix(1, 2), 7.0);
    EXPECT_DOUBLE_EQ(sparse_matrix(0, 0), 0.0);
    // the second R
    sparse_matrix = csr.getMatrix(1);
    sparse_matrix1 = csr.getMatrix(0, 0, 0);
    for (const auto& element: sparse_matrix.getElements())
    {
        auto it = sparse_matrix1.getElements().find(element.first);
        EXPECT_EQ(it->first.first, element.first.first);
        EXPECT_EQ(it->first.second, element.first.second);
        EXPECT_DOUBLE_EQ(it->second, element.second);
        // std::cout << "element( " << element.first.first << ", " << element.first.second << " ) = " << element.second
        // << std::endl;
    }
    EXPECT_DOUBLE_EQ(sparse_matrix(2, 2), 5.0);
    EXPECT_DOUBLE_EQ(sparse_matrix(2, 3), 6.0);
    EXPECT_DOUBLE_EQ(sparse_matrix(3, 3), 10.0);
    EXPECT_DOUBLE_EQ(sparse_matrix(0, 0), 0.0);
}

TEST_F(csrFileReaderTest, ComplexCsrReader)
{
    ModuleIO::csrFileReader<std::complex<double>> csr(filename);

    EXPECT_TRUE(csr.isOpen());
    EXPECT_EQ(csr.getStep(), 1);
    EXPECT_EQ(csr.getMatrixDimension(), 4);
    EXPECT_EQ(csr.getNumberOfR(), 2);
    EXPECT_EQ(csr.getRCoordinate(0), std::vector<int>({0, 1, 1}));
    EXPECT_EQ(csr.getRCoordinate(1), std::vector<int>({0, 0, 0}));

    const ModuleIO::SparseMatrix<std::complex<double>> first_by_index = csr.getMatrix(0);
    const ModuleIO::SparseMatrix<std::complex<double>> first_by_coordinate = csr.getMatrix(0, 1, 1);
    EXPECT_EQ(first_by_index.getElements(), first_by_coordinate.getElements());
    EXPECT_EQ(first_by_index(0, 3), std::complex<double>(4.0, 0.0));
    EXPECT_EQ(first_by_index(1, 2), std::complex<double>(7.0, 0.0));
    EXPECT_EQ(first_by_index(0, 0), std::complex<double>(0.0, 0.0));

    const ModuleIO::SparseMatrix<std::complex<double>> second_by_index = csr.getMatrix(1);
    const ModuleIO::SparseMatrix<std::complex<double>> second_by_coordinate = csr.getMatrix(0, 0, 0);
    EXPECT_EQ(second_by_index.getElements(), second_by_coordinate.getElements());
    EXPECT_EQ(second_by_index(2, 2), std::complex<double>(5.0, 0.0));
    EXPECT_EQ(second_by_index(2, 3), std::complex<double>(6.0, 0.0));
    EXPECT_EQ(second_by_index(3, 3), std::complex<double>(10.0, 0.0));
    EXPECT_EQ(second_by_index(0, 0), std::complex<double>(0.0, 0.0));
}
