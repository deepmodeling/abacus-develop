#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_LINALG_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_LINALG_H

#include "source_base/matrix2d.h"
#include "source_base/matrix2d_view.h"

#include <complex>

namespace matrix2d_linalg
{
/**
 * @brief Invert a square matrix in place using an explicit matrix inverse.
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

/**
 * @brief Multiply two matrix views into result.
 *
 * The result view must not overlap either input view. The host implementation
 * writes result while still reading left and right.
 */
void matmul(Matrix2dView<double> result,
            Matrix2dView<const double> left,
            Matrix2dView<const double> right);
/**
 * @copydoc matmul(Matrix2dView<double>,Matrix2dView<const double>,Matrix2dView<const double>)
 */
void matmul(Matrix2dView<std::complex<double>> result,
            Matrix2dView<const std::complex<double>> left,
            Matrix2dView<const std::complex<double>> right);

/**
 * @brief Compute `result = alpha * left + beta * right`.
 *
 * All views must have the same global shape and compatible host or distributed
 * layouts. The result may overlap either input because each element is read
 * before it is written.
 */
void axpby(Matrix2dView<double> result,
           const double alpha,
           Matrix2dView<const double> left,
           const double beta,
           Matrix2dView<const double> right);
/**
 * @copydoc axpby(Matrix2dView<double>,double,Matrix2dView<const double>,double,Matrix2dView<const double>)
 */
void axpby(
    Matrix2dView<std::complex<double>> result,
    const std::complex<double> alpha,
    Matrix2dView<const std::complex<double>> left,
    const std::complex<double> beta,
    Matrix2dView<const std::complex<double>> right);

void scale(Matrix2dView<double> matrix, const double factor);
void scale(Matrix2dView<std::complex<double>> matrix,
           const std::complex<double> factor);

/**
 * @brief Write the identity matrix into an existing square view.
 *
 * This operation mutates `matrix` in place and does not allocate a new matrix.
 */
void set_identity(Matrix2dView<double> matrix);
/** @copydoc set_identity(Matrix2dView<double>) */
void set_identity(Matrix2dView<std::complex<double>> matrix);

/**
 * @brief Create a new `n` by `n` identity matrix.
 *
 * The returned matrix owns host storage. Use `set_identity()` when an existing
 * view should be overwritten in place.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
Matrix2d<T> eye(const int n);
} // namespace matrix2d_linalg

#endif
