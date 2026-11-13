#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_COLLECTIVE_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_COLLECTIVE_H

#include "source_base/matrix2d.h"

#include <memory>

#ifdef __MPI
namespace matrix2d_impl
{
template <typename T>
Matrix2d<T> allgather(const Matrix2d<T>& matrix, const MPI_Comm comm);

template <typename T>
std::unique_ptr<Matrix2d<T>> gather_on_root(const Matrix2d<T>& matrix,
                                            const int iproot,
                                            const MPI_Comm comm);

template <typename T>
void scatter_from_host(const Matrix2d<T>& matglb, Matrix2d<T>& matloc);
} // namespace matrix2d_impl
#endif

#endif
