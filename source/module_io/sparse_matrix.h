#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

namespace ModuleIO
{

/**
 * @brief Sparse matrix class
 *
 * @tparam T
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
    void addValue(int row, int col, T value);

    // get data vector
    const std::vector<std::tuple<int, int, T>>& getData() const
    {
        return data;
    }

    // print data in CSR (Compressed Sparse Row) format
    void printToCSR(std::ostream& ofs, double threshold, int precision = 8);

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

  private:
    int _rows;
    int _cols;
    std::vector<std::tuple<int, int, T>> data;
    std::vector<int> csr_row_ptr;
}; // class SparseMatrix

} // namespace ModuleIO

#endif // SPARSE_MATRIX_H