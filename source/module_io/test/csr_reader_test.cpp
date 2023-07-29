#include "module_io/csr_reader.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

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
    EXPECT_EQ(csr.getStep(), 0);
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
    for (const auto& element : sparse_matrix.getElements())
    {
        auto it = sparse_matrix1.getElements().find(element.first);
        EXPECT_EQ(it->first.first, element.first.first);
        EXPECT_EQ(it->first.second, element.first.second);
        EXPECT_DOUBLE_EQ(it->second, element.second);
        //std::cout << "element( " << element.first.first << ", " << element.first.second << " ) = " << element.second << std::endl;
    }
    EXPECT_DOUBLE_EQ(sparse_matrix(0, 3), 4.0);
    EXPECT_DOUBLE_EQ(sparse_matrix(1, 2), 7.0);
    EXPECT_DOUBLE_EQ(sparse_matrix(0, 0), 0.0);
    // the second R
    sparse_matrix = csr.getMatrix(1);
    sparse_matrix1 = csr.getMatrix(0, 0, 0);
    for (const auto& element : sparse_matrix.getElements())
    {
        auto it = sparse_matrix1.getElements().find(element.first);
        EXPECT_EQ(it->first.first, element.first.first);
        EXPECT_EQ(it->first.second, element.first.second);
        EXPECT_DOUBLE_EQ(it->second, element.second);
        //std::cout << "element( " << element.first.first << ", " << element.first.second << " ) = " << element.second << std::endl;
    }
    EXPECT_DOUBLE_EQ(sparse_matrix(2, 2), 5.0);
    EXPECT_DOUBLE_EQ(sparse_matrix(2, 3), 6.0);
    EXPECT_DOUBLE_EQ(sparse_matrix(3, 3), 10.0);
    EXPECT_DOUBLE_EQ(sparse_matrix(0, 0), 0.0);
}