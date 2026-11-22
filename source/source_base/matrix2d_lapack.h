#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_LAPACK_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_LAPACK_H

/**
 * @file matrix2d_lapack.h
 * @internal
 * @brief Matrix2d LAPACK backend declarations.
 *
 * Callers should include `source_base/matrix2d_linalg.h` instead of this
 * header.
 */

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
