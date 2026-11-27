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
/**
 * @brief Invert a square host matrix in place with GETRF and GETRI.
 *
 * @warning Explicit inversion has poorer numerical stability than solving a
 * linear system directly. Use this function only when the inverse matrix
 * itself is required. When solving `A * X = B`, call `solve()` instead.
 */
void inv(Matrix2dView<double> matrix);
/** @copydoc inv(Matrix2dView<double>) */
void inv(Matrix2dView<std::complex<double>> matrix);

void solve(Matrix2dView<double> matrix,
           Matrix2dView<double> rhs);
void solve(Matrix2dView<std::complex<double>> matrix,
           Matrix2dView<std::complex<double>> rhs);
} // namespace lapack
} // namespace matrix2d_impl

#endif
