#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_LINALG_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_LINALG_H

#include "source_base/matrix2d.h"
#include "source_base/matrix2d_view.h"

#include <complex>

namespace matrix2d_linalg
{
/**
 * Value-returning overloads return `Matrix2d<T>` by value. `Matrix2d<T>` has
 * defaulted copy and move operations, so C++11 return-value optimization or
 * the move constructor transfers ownership without a deep result copy.
 */

/**
 * @brief Return the inverse of a square matrix.
 *
 * The input matrix is unchanged. The result has the same host or frozen
 * distributed layout as the input and is collective when the input is an
 * active MPI distribution.
 *
 * @warning Explicit inversion has poorer numerical stability than solving a
 * linear system directly. Use this function only when the inverse matrix
 * itself is required. When solving `A * X = B`, call `solve()` instead.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
Matrix2d<T> inv(const Matrix2d<T>& matrix);

/**
 * @brief Write the inverse of a square matrix into an existing matrix.
 *
 * The input matrix is unchanged. `inverse` receives a complete result with
 * the same host or frozen distributed layout as the input.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
void inv(const Matrix2d<T>& matrix, Matrix2d<T>& inverse);

/**
 * @brief Invert a square matrix in place without copying.
 *
 * @param matrix `inout` square matrix. It is overwritten by its inverse.
 *
 * This is a low-level escape hatch over caller-owned storage. Prefer `inv()`
 * unless avoiding the input copy is essential.
 */
void inv_inplace(Matrix2dView<double> matrix);
/** @copydoc inv_inplace(Matrix2dView<double>) */
void inv_inplace(Matrix2dView<std::complex<double>> matrix);

/**
 * @brief Return the solution X of `A * X = B`.
 *
 * `matrix` and `rhs` are unchanged. The result has the same host or frozen
 * distributed layout as `rhs`. When `rhs` is actively distributed, this
 * operation is collective and every rank must call it.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
Matrix2d<T> solve(const Matrix2d<T>& matrix,
                  const Matrix2d<T>& rhs);

/**
 * @brief Write the solution X of `A * X = B` into an existing matrix.
 *
 * `matrix` and `rhs` are unchanged. `solution` receives a complete result
 * with the same host or frozen distributed layout as `rhs`. The output may
 * alias either input.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
void solve(const Matrix2d<T>& matrix,
           const Matrix2d<T>& rhs,
           Matrix2d<T>& solution);

/**
 * @brief Solve `A * X = B` in place without copying.
 *
 * @param matrix `inout` coefficient matrix. It is overwritten by its LU
 * factorization.
 * @param rhs `inout` right-hand side. It is overwritten by the solution X.
 *
 * This is a low-level escape hatch for callers that explicitly accept the
 * mutation. Prefer `solve()` unless avoiding both input copies is essential.
 * The operation is collective for active MPI distributions.
 */
void solve_inplace(Matrix2dView<double> matrix,
                   Matrix2dView<double> rhs);
/** @copydoc solve_inplace(Matrix2dView<double>,Matrix2dView<double>) */
void solve_inplace(Matrix2dView<std::complex<double>> matrix,
                   Matrix2dView<std::complex<double>> rhs);

/**
 * @brief Return the matrix product `A * B`.
 *
 * Both inputs are unchanged. The result has a host layout when the inputs are
 * host matrices. For distributed inputs it carries the left input's process
 * grid, block size, and BLACS context and the operation is collective.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
Matrix2d<T> matmul(const Matrix2d<T>& left,
                   const Matrix2d<T>& right);

/**
 * @brief Write the matrix product `A * B` into an existing matrix.
 *
 * Both inputs are unchanged. `result` receives a result with the layout rules
 * described by `matmul(left, right)`. The output may alias either input.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
void matmul(const Matrix2d<T>& left,
            const Matrix2d<T>& right,
            Matrix2d<T>& result);

/**
 * @brief Multiply two matrix views into result.
 *
 * The result storage must not overlap either input. Exact alias as well as
 * shifted or partial overlap is unsupported because the implementation writes
 * result while still reading left and right.
 *
 * @param result `out` matrix receiving the product.
 * @param left Read-only left input.
 * @param right Read-only right input.
 */
void matmul_inplace(Matrix2dView<double> result,
                    Matrix2dView<const double> left,
                    Matrix2dView<const double> right);
/**
 * @copydoc matmul_inplace(Matrix2dView<double>,Matrix2dView<const double>,Matrix2dView<const double>)
 */
void matmul_inplace(Matrix2dView<std::complex<double>> result,
                    Matrix2dView<const std::complex<double>> left,
                    Matrix2dView<const std::complex<double>> right);

/**
 * @brief Return `alpha * left + beta * right`.
 *
 * The result is newly allocated with the same layout as `left`. Inputs are
 * unchanged.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
Matrix2d<T> axpby(const T alpha,
                  const Matrix2d<T>& left,
                  const T beta,
                  const Matrix2d<T>& right);

/**
 * @brief Write `alpha * left + beta * right` into an existing matrix.
 *
 * The output may alias either input. It receives the layout of `left`.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
void axpby(const T alpha,
           const Matrix2d<T>& left,
           const T beta,
           const Matrix2d<T>& right,
           Matrix2d<T>& result);

/**
 * @brief Compute `result = alpha * left + beta * right` through views.
 *
 * All views must have the same global shape and compatible host or distributed
 * layouts. The result may alias either input only when the aliased views use
 * the same element mapping for every logical coordinate, including a fully
 * in-place call such as `result == left`. Shifted and partial overlap are
 * unsupported because a write may invalidate an input read by another element
 * or OpenMP thread.
 *
 * @param result `inout` matrix receiving the combination.
 */
void axpby_inplace(Matrix2dView<double> result,
                   const double alpha,
                   Matrix2dView<const double> left,
                   const double beta,
                   Matrix2dView<const double> right);
/**
 * @copydoc axpby_inplace(Matrix2dView<double>,double,Matrix2dView<const double>,double,Matrix2dView<const double>)
 */
void axpby_inplace(
    Matrix2dView<std::complex<double>> result,
    const std::complex<double> alpha,
    Matrix2dView<const std::complex<double>> left,
    const std::complex<double> beta,
    Matrix2dView<const std::complex<double>> right);

/**
 * @brief Return a newly allocated copy of `matrix` scaled by `factor`.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
Matrix2d<T> scale(const Matrix2d<T>& matrix, const T factor);

/**
 * @brief Write a scaled copy of `matrix` into `result`.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
void scale(const Matrix2d<T>& matrix,
           const T factor,
           Matrix2d<T>& result);

/**
 * @brief Scale a matrix view in place without copying.
 *
 * @param matrix `inout` matrix to scale.
 * @param factor Scalar multiplier.
 */
void scale_inplace(Matrix2dView<double> matrix, const double factor);
/** @copydoc scale_inplace(Matrix2dView<double>,double) */
void scale_inplace(Matrix2dView<std::complex<double>> matrix,
                   const std::complex<double> factor);

/**
 * @brief Write the identity matrix into an existing square view.
 *
 * @param matrix `inout` square matrix receiving the identity.
 *
 * This is a low-level zero-copy operation. Prefer `eye()` when a new identity
 * matrix is required.
 */
void set_identity_inplace(Matrix2dView<double> matrix);
/** @copydoc set_identity_inplace(Matrix2dView<double>) */
void set_identity_inplace(Matrix2dView<std::complex<double>> matrix);

/**
 * @brief Create a new `n` by `n` identity matrix.
 *
 * The returned matrix owns host storage. Use `set_identity_inplace()` only
 * when an existing view should be overwritten without allocating.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
Matrix2d<T> eye(const int n);

/**
 * @brief Return the non-conjugate transpose of `matrix`.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
Matrix2d<T> transpose(const Matrix2d<T>& matrix);

/**
 * @brief Return the conjugate transpose of `matrix`.
 *
 * @tparam T Either `double` or `std::complex<double>`.
 */
template <typename T>
Matrix2d<T> dagger(const Matrix2d<T>& matrix);
} // namespace matrix2d_linalg

#endif
