#ifndef MATRIXBLOCK_H
#define MATRIXBLOCK_H

#include <cstddef>
namespace hamilt
{

template <typename T> struct MatrixBlock
{
    /* this is a simple template block of a matrix
       would change to Eigen in the future */
    T* p;
    size_t row;
    size_t col;
    int* desc;
};

} // namespace hamilt
#endif