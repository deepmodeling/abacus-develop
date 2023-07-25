#include "sparse_matrix.h"

#include <complex>

#include "module_base/tool_quit.h"

namespace ModuleIO
{

/**
 * @brief Add value to the matrix with row and column indices
 */
template <typename T>
void SparseMatrix<T>::addValue(int row, int col, T value)
{
    if (row < 0 || row >= _rows || col < 0 || col >= _cols)
    {
        ModuleBase::WARNING_QUIT("SparseMatrix::addValue", "row or col index out of range");
    }
    data.push_back(std::make_tuple(row, col, value));
}

/**
 * @brief Print to CSR format
 */
template <typename T>
void SparseMatrix<T>::printToCSR(std::ostream& ofs, double threshold, int precision)
{
    // Filter elements greater than the threshold
    auto it = std::remove_if(data.begin(), data.end(), [threshold](const std::tuple<int, int, T>& elem) {
        return std::abs(std::get<2>(elem)) <= threshold;
    });
    data.erase(it, data.end());

    // Sort data by row, then by column
    std::sort(data.begin(), data.end(), [](const std::tuple<int, int, T>& a, const std::tuple<int, int, T>& b) {
        if (std::get<0>(a) == std::get<0>(b))
        {
            return std::get<1>(a) < std::get<1>(b);
        }
        return std::get<0>(a) < std::get<0>(b);
    });

    // Initialize the CSR arrays
    std::vector<int> csr_row_ptr;
    csr_row_ptr.assign(_rows + 1, 0);

    // print the CSR values
    for (int i = 0; i < data.size(); i++)
    {
        ofs << " " << std::fixed << std::scientific << std::setprecision(precision) << std::get<2>(data[i]);
    }
    ofs << std::endl;
    // print the CSR column indices
    for (int i = 0; i < data.size(); i++)
    {
        ofs << " " << std::get<1>(data[i]);
        int row = std::get<0>(data[i]);
        csr_row_ptr[row + 1]++;
    }
    ofs << std::endl;

    // Compute the row pointers
    for (int i = 1; i <= _rows; i++)
    {
        csr_row_ptr[i] += csr_row_ptr[i - 1];
    }

    // print the CSR row pointers
    for (int i = 0; i < csr_row_ptr.size(); i++)
    {
        ofs << " " << csr_row_ptr[i];
    }
    ofs << std::endl;
}

/**
 * @brief Read CSR data from arrays
 */
template <typename T>
void SparseMatrix<T>::readCSR(const std::vector<T>& values,
                              const std::vector<int>& col_ind,
                              const std::vector<int>& row_ptr)
{
    if (row_ptr.size() != static_cast<size_t>(_rows) + 1)
    {
        ModuleBase::WARNING_QUIT("SparseMatrix::readCSR", "Invalid row_ptr size");
    }
    if (col_ind.size() != values.size())
    {
        ModuleBase::WARNING_QUIT("SparseMatrix::readCSR", "Column indices and values size mismatch");
    }

    data.clear();
    for (int row = 0; row < _rows; row++)
    {
        for (int idx = row_ptr[row]; idx < row_ptr[row + 1]; idx++)
        {
            data.push_back(std::make_tuple(row, col_ind[idx], values[idx]));
        }
    }
}

// Explicit instantiation of template classes
template class SparseMatrix<double>;
template class SparseMatrix<std::complex<double>>;

} // namespace ModuleIO