#ifndef MODULE_BASE_MATRIX_BLOCK_H
#define MODULE_BASE_MATRIX_BLOCK_H

#include <cstddef>

namespace ModuleBase
{

/**
 * @brief A non-owning description of a matrix stored in memory.
 *
 * It records only where the matrix lives and how it is laid out: the data
 * pointer, the local number of rows and columns, and the BLACS array
 * descriptor when the matrix is block-cyclically distributed. It carries no
 * physical meaning, so it can be shared between the code that fills a matrix
 * and the code that diagonalizes it without either side depending on the
 * other.
 *
 * @note This is an aggregate on purpose; several call sites brace-initialize
 * it as MatrixBlock<T>{p, row, col, desc}.
 */
template <typename T> struct MatrixBlock
{
    /* would change to Eigen in the future */
    T* p;
    size_t row;
    size_t col;
    const int* desc;
};

} // namespace ModuleBase

#endif // MODULE_BASE_MATRIX_BLOCK_H
