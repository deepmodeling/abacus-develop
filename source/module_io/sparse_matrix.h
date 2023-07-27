#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <iostream>
#include <map>
#include <utility>
#include <vector>

namespace ModuleIO
{

/**
 * @brief Sparse matrix class designed mainly for csr format input and output.
 * @details
 *   The sparse matrix is stored in a map.
 *   The map key is a pair of row and column indices.
 *   The map value is the matrix element.
 *   The matrix element is stored only if its absolute value is greater than the threshold.
 *   The threshold is set to 1.0e-10 by default.
 * @tparam T data type, it can be double or std::complex<double>
 */
template <typename T>
class SparseMatrix
{
  public:
    // Default constructor
    SparseMatrix() : _rows(0), _cols(0)
    {
    }

    SparseMatrix(int rows, int cols) : _rows(rows), _cols(cols)
    {
    }

    // add value to the matrix with row and column indices
    void insert(int row, int col, T value);

    // print data in CSR (Compressed Sparse Row) format
    void printToCSR(std::ostream& ofs, int precision = 8);

    // read CSR data from arrays
    void readCSR(const std::vector<T>& values, const std::vector<int>& col_ind, const std::vector<int>& row_ptr);

    // set number of rows
    void setRows(int rows)
    {
        _rows = rows;
    }

    // set number of columns
    void setCols(int cols)
    {
        _cols = cols;
    }

    // get number of rows
    int getRows() const
    {
        return _rows;
    }

    // get number of columns
    int getCols() const
    {
        return _cols;
    }

    // define the operator to index a matrix element
    T operator()(int row, int col);

    // set the threshold
    void setSparseThreshold(double sparse_threshold)
    {
        _sparse_threshold = sparse_threshold;
    }

    // get the threshold
    double getSparseThreshold() const
    {
        return _sparse_threshold;
    }

    // get the number of non-zero elements
    int getNNZ() const
    {
        return elements.size();
    }

    // get elements
    const std::map<std::pair<int, int>, T>& getElements() const
    {
        return elements;
    }

  private:
    int _rows;
    int _cols;
    std::map<std::pair<int, int>, T> elements;
    double _sparse_threshold = 1.0e-10;
}; // class SparseMatrix

} // namespace ModuleIO

#endif // SPARSE_MATRIX_H