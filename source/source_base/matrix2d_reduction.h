#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_REDUCTION_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_REDUCTION_H

#include "source_base/matrix2d.h"

namespace matrix2d_impl
{
template <typename T>
T sum(const Matrix2d<T>& matrix);

template <typename T>
double sum_squares(const Matrix2d<T>& matrix);

template <typename T>
T dot(const Matrix2d<T>& left, const Matrix2d<T>& right);

template <typename T>
double norm2(const Matrix2d<T>& matrix);

template <typename T>
double max_abs(const Matrix2d<T>& matrix);

template <typename T>
std::vector<T> diagonal(const Matrix2d<T>& matrix);
} // namespace matrix2d_impl

#endif
