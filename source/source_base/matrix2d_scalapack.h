#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_SCALAPACK_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_SCALAPACK_H

/**
 * @file matrix2d_scalapack.h
 * @internal
 * @brief Matrix2d ScaLAPACK backend declarations.
 *
 * Callers should include `source_base/matrix2d_linalg.h` instead of this
 * header.
 */

#include <complex>

template <typename T>
class Matrix2dView;

#ifdef __MPI
namespace matrix2d_impl
{
namespace scalapack
{
void assert_compatible_layout(const int* left_desc, const int* right_desc);

void transpose(Matrix2dView<double> result,
               Matrix2dView<const double> source);
void transpose(Matrix2dView<std::complex<double>> result,
               Matrix2dView<const std::complex<double>> source);
void conjugate_transpose(Matrix2dView<double> result,
                         Matrix2dView<const double> source);
void conjugate_transpose(
    Matrix2dView<std::complex<double>> result,
    Matrix2dView<const std::complex<double>> source);
/**
 * @brief Invert a square distributed matrix in place with P?GETRF and
 * P?GETRI.
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
void matmul(Matrix2dView<double> result,
            Matrix2dView<const double> left,
            Matrix2dView<const double> right);
void matmul(Matrix2dView<std::complex<double>> result,
            Matrix2dView<const std::complex<double>> left,
            Matrix2dView<const std::complex<double>> right);

} // namespace scalapack
} // namespace matrix2d_impl
#endif

#endif
