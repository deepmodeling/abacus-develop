#ifndef CSR_READER_H
#define CSR_READER_H

#include <fstream>

#include "file_reader.h"
#include "sparse_matrix.h"

namespace ModuleIO
{

/**
 * @brief Class to read CSR file
 * @details This class is used to read CSR file
 * the default format is like:
 * ```
 * STEP: 0
 * Matrix Dimension of S(R): 4
 * Matrix number of S(R): 2
 * 0 1 1 2            # (0,1,1) is the R coordinate, 2 is the number of non-zero elements
 *  4.00e+00 7.00e+00 # non-zero elements
 *  3 2               # column indices
 *  0 1 2 2 2         # row pointer
 * 0 0 0 3            # the second R coordinate and number of non-zero elements
 * 5.00e+00 6.00e+00 1.00e+01
 * 2 3 3
 * 0 0 0 2 3
 * ```
 * It will store the R coordinates and sparse matrices as two vectors.
 * One can use getter functions to get the R coordinates and sparse matrices,
 * and related info including step, matrix dimension, number of R.
 */
template <typename T>
class csrFileReader : public FileReader
{
  public:
    // Constructor
    csrFileReader(const std::string& filename);

    // read all matrices of all R coordinates
    void parseFile();

    // get number of R
    int getNumberOfR() const;

    // get sparse matrix of a specific R coordinate
    SparseMatrix<T> getMatrix(int Rx, int Ry, int Rz);

    // get matrix by using index
    SparseMatrix<T> getMatrix(int index) const;

    // get R coordinate using index
    std::vector<int> getRCoordinate(int index) const;

    // get step
    int getStep() const;

    // get matrix dimension
    int getMatrixDimension() const;

  private:
    std::vector<std::vector<int>> RCoordinates;
    std::vector<SparseMatrix<T>> sparse_matrices;
    int step;
    int matrixDimension;
    int numberOfR;
};

} // namespace ModuleIO

#endif // READ_CSR_H