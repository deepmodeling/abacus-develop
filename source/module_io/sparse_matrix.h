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

    // get CSR_values
    const std::vector<T>& getCSRValues() const
    {
        return csr_values;
    }

    // get CSR_col_ind
    const std::vector<int>& getCSRColInd() const
    {
        return csr_col_ind;
    }

    // get CSR_row_ptr
    const std::vector<int>& getCSRRowPtr() const
    {
        return csr_row_ptr;
    }

    // convert to CSR format
    void convertToCSR(double threshold);

    // print data in CSR format  CSR: Compressed Sparse Row
    void printCSR(std::ostream& ofs, int precision = 8);

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
    std::vector<int> csr_col_ind;
    std::vector<T> csr_values;
}; // class SparseMatrix

} // namespace ModuleIO

#endif // SPARSE_MATRIX_H