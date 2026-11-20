#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_LAPACK_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_LAPACK_H

#include <complex>

template <typename T>
class Matrix2dView;

namespace matrix2d_impl
{
namespace lapack
{
void inv(Matrix2dView<double> matrix);
void inv(Matrix2dView<std::complex<double>> matrix);

void solve(Matrix2dView<double> matrix,
           Matrix2dView<double> rhs);
void solve(Matrix2dView<std::complex<double>> matrix,
           Matrix2dView<std::complex<double>> rhs);
} // namespace lapack
} // namespace matrix2d_impl

#endif
