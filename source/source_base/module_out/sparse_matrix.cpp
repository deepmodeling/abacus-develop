#include "sparse_matrix.h"

#include <algorithm>
#include <complex>

#include "source_base/tool_quit.h"

namespace ModuleIO
{

/**
 * @brief Add value to the matrix with row and column indices
 */
template <typename T>
void SparseMatrix<T>::insert(int row, int col, T value)
{
    if (row < 0 || row >= _rows || col < 0 || col >= _cols)
    {
        ModuleBase::WARNING_QUIT("SparseMatrix::addValue", "row or col index out of range");
    }
    if (std::abs(value) > _sparse_threshold)
    {
        elements[std::make_pair(row, col)] = value;
    }
}

/**
 * @brief Print to CSR format
 */
template <typename T>
void SparseMatrix<T>::printToCSR(std::ostream& ofs, int precision)
{
    // Initialize the CSR arrays
    std::vector<int> csr_row_ptr;
    csr_row_ptr.assign(_rows + 1, 0);

    ofs << std::scientific << std::setprecision(precision);

    // print the CSR values
    ofs << " # CSR values";
    size_t count1 = 0;
    for (const auto &element : elements)
    {
	if(count1%6==0) ofs << std::endl;
	count1++;
        ofs << " " << element.second;
    }
    ofs << std::endl;
    // print the CSR column indices
    ofs << " # CSR column indices";
    size_t count2 = 0;
    for (const auto &element : elements)
    {
	if(count2%16==0) ofs << std::endl;
	count2++;
        ofs << " " << element.first.second;
        int row = element.first.first;
        csr_row_ptr[row + 1]++;
    }
    ofs << std::endl;

    // Compute the row pointers
    for (int i = 1; i <= _rows; i++)
    {
        csr_row_ptr[i] += csr_row_ptr[i - 1];
    }

    // print the CSR row pointers
    ofs << " # CSR row pointers";
    for (int i = 0; i < csr_row_ptr.size(); i++)
    {
        if(i%16==0) ofs << std::endl;
        ofs << " " << csr_row_ptr[i];
    }
    ofs << std::endl << std::endl;
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

    elements.clear();
    for (int row = 0; row < _rows; row++)
    {
        for (int idx = row_ptr[row]; idx < row_ptr[row + 1]; idx++)
        {
            elements[std::make_pair(row, col_ind[idx])] = values[idx];
        }
    }
}

// define the operator to index a matrix element
template <typename T>
T SparseMatrix<T>::operator()(int row, int col) const
{
    if (row < 0 || row >= _rows || col < 0 || col >= _cols)
    {
        ModuleBase::WARNING_QUIT("SparseMatrix::operator()", "row or col index out of range");
    }
    auto it = elements.find(std::make_pair(row, col));
    if (it != elements.end())
    {
        return it->second;
    }
    else
    {
        return static_cast<T>(0);
    }
}

// Explicit instantiation of template classes
template class SparseMatrix<double>;
template class SparseMatrix<std::complex<double>>;

} // namespace ModuleIO
