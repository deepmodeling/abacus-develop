#include "source_base/matrix2d_linalg.h"

#include "source_base/matrix2d_lapack.h"
#include "source_base/matrix2d_view.h"
#include "source_base/parallel_blkcyc.h"
#ifdef __MPI
#include "source_base/matrix2d_scalapack.h"
#endif

#include <complex>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace
{
/**
 * @brief Return true when a distribution handle exists and is not the serial
 * one-by-one layout.
 */
template <typename T>
bool distribution_is_active(const Matrix2dView<T>& matrix)
{
    return matrix.is_active();
}

/**
 * @brief Return true when a distribution handle exists and describes the
 * serial one-by-one layout.
 */
template <typename T>
bool distribution_is_serial(const Matrix2dView<T>& matrix)
{
    return matrix.has_distribution() && matrix.distribution()->is_serial();
}

template <typename T, typename U>
void assert_same_shape(const Matrix2dView<T>& left,
                       const Matrix2dView<U>& right)
{
    if (left.nrglb() != right.nrglb() || left.ncglb() != right.ncglb())
    {
        throw std::invalid_argument("matrix shapes are incompatible");
    }
}

template <typename T, typename U>
void assert_compatible_layout(const Matrix2dView<T>& left,
                              const Matrix2dView<U>& right)
{
    if (left.has_distribution() != right.has_distribution())
    {
        throw std::invalid_argument(
            "matrix operations require matching host or distributed layouts");
    }
    if (!left.has_distribution())
    {
        return;
    }
    const bool serial = distribution_is_serial(left);
    const bool serial_ = distribution_is_serial(right);
    if (serial != serial_)
    {
        throw std::invalid_argument(
            "matrix operations require matching serial or parallel distributions");
    }
    if (serial)
    {
        return;
    }
#ifdef __MPI
    matrix2d_impl::scalapack::assert_compatible_layout(left.desc(),
                                                       right.desc());
#else
    throw std::invalid_argument(
        "matrix operations require MPI-compatible distributed layouts");
#endif
}

template <typename T>
void assert_matmul_views(Matrix2dView<T> result,
                         Matrix2dView<const T> left,
                         Matrix2dView<const T> right)
{
    if (left.ncglb() != right.nrglb()
        || result.nrglb() != left.nrglb()
        || result.ncglb() != right.ncglb())
    {
        throw std::invalid_argument("matrix multiplication dimensions are incompatible");
    }
    assert_compatible_layout(result, left);
    assert_compatible_layout(result, right);
}

template <typename T>
void matmul_host(Matrix2dView<T> result,
                 Matrix2dView<const T> left,
                 Matrix2dView<const T> right)
{
    for (int i = 0; i < result.nrglb(); ++i)
    {
        for (int j = 0; j < result.ncglb(); ++j)
        {
            T value = T();
            for (int k = 0; k < left.ncglb(); ++k)
            {
                value += left.atloc(i, k) * right.atloc(k, j);
            }
            result.atloc(i, j) = value;
        }
    }
}

template <typename T>
void axpby_impl(Matrix2dView<T> result,
                const T alpha,
                Matrix2dView<const T> left,
                const T beta,
                Matrix2dView<const T> right)
{
    assert_same_shape(result, left);
    assert_same_shape(result, right);
    assert_compatible_layout(result, left);
    assert_compatible_layout(result, right);
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i = 0; i < result.nrloc(); ++i)
    {
        for (int j = 0; j < result.ncloc(); ++j)
        {
            result.atloc(i, j)
                = alpha * left.atloc(i, j) + beta * right.atloc(i, j);
        }
    }
}

template <typename T>
void scale_impl(Matrix2dView<T> matrix, const T factor)
{
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static)
#endif
    for (int i = 0; i < matrix.nrloc(); ++i)
    {
        for (int j = 0; j < matrix.ncloc(); ++j)
        {
            matrix.atloc(i, j) *= factor;
        }
    }
}

template <typename T>
void set_identity_impl(Matrix2dView<T> matrix)
{
    if (matrix.nrglb() != matrix.ncglb())
    {
        throw std::invalid_argument("identity requires a square matrix");
    }
    const std::shared_ptr<const Parallel_Blkcyclic> distrib
        = matrix.distribution();
    for (int i = 0; i < matrix.nrloc(); ++i)
    {
        for (int j = 0; j < matrix.ncloc(); ++j)
        {
            const int irglb = distrib == nullptr ? i : distrib->irglb(i);
            const int icglb = distrib == nullptr ? j : distrib->icglb(j);
            matrix.atloc(i, j) = irglb == icglb ? T(1.0) : T();
        }
    }
}

void assert_inverse_view(const Matrix2dView<double>& matrix)
{
    if (matrix.nrglb() <= 0 || matrix.nrglb() != matrix.ncglb())
    {
        throw std::invalid_argument("matrix inverse requires a non-empty square matrix");
    }
}

void assert_inverse_view(const Matrix2dView<std::complex<double>>& matrix)
{
    if (matrix.nrglb() <= 0 || matrix.nrglb() != matrix.ncglb())
    {
        throw std::invalid_argument("matrix inverse requires a non-empty square matrix");
    }
}

void assert_solve_views(const Matrix2dView<double>& matrix,
                        const Matrix2dView<double>& rhs)
{
    if (matrix.nrglb() <= 0
        || matrix.nrglb() != matrix.ncglb()
        || rhs.nrglb() != matrix.nrglb())
    {
        throw std::invalid_argument(
            "linear solve requires a square matrix and matching right-hand side");
    }
    if (matrix.has_distribution() != rhs.has_distribution())
    {
        throw std::invalid_argument(
            "linear solve requires matching host or distributed layouts");
    }
    if (distribution_is_serial(matrix) != distribution_is_serial(rhs))
    {
        throw std::invalid_argument(
            "linear solve requires matching serial or parallel distributions");
    }
    if (distribution_is_serial(matrix))
    {
        return;
    }
    if (!matrix.has_distribution())
    {
        return;
    }
#ifdef __MPI
    matrix2d_impl::scalapack::assert_compatible_layout(matrix.desc(),
                                                       rhs.desc());
#else
    throw std::invalid_argument(
        "linear solve requires MPI-compatible distributed layouts");
#endif
}

void assert_solve_views(const Matrix2dView<std::complex<double>>& matrix,
                        const Matrix2dView<std::complex<double>>& rhs)
{
    if (matrix.nrglb() <= 0
        || matrix.nrglb() != matrix.ncglb()
        || rhs.nrglb() != matrix.nrglb())
    {
        throw std::invalid_argument(
            "linear solve requires a square matrix and matching right-hand side");
    }
    if (matrix.has_distribution() != rhs.has_distribution())
    {
        throw std::invalid_argument(
            "linear solve requires matching host or distributed layouts");
    }
    if (distribution_is_serial(matrix) != distribution_is_serial(rhs))
    {
        throw std::invalid_argument(
            "linear solve requires matching serial or parallel distributions");
    }
    if (distribution_is_serial(matrix))
    {
        return;
    }
    if (!matrix.has_distribution())
    {
        return;
    }
#ifdef __MPI
    matrix2d_impl::scalapack::assert_compatible_layout(matrix.desc(),
                                                       rhs.desc());
#else
    throw std::invalid_argument(
        "linear solve requires MPI-compatible distributed layouts");
#endif
}

template <typename T>
void assert_supported_linalg_type()
{
    static_assert(std::is_same<T, double>::value
                      || std::is_same<T, std::complex<double>>::value,
                  "linear algebra supports double and std::complex<double>");
}

template <typename T>
Matrix2d<T> make_matrix_like(const Matrix2d<T>& source,
                             const int nrglb,
                             const int ncglb)
{
    if (!source.has_distribution())
    {
        return Matrix2d<T>(nrglb, ncglb);
    }
    if (!source.is_active())
    {
        return Matrix2d<T>(make_frozen_distribution(nrglb, ncglb));
    }
#ifdef __MPI
    if (source.nrglb() == nrglb && source.ncglb() == ncglb)
    {
        return Matrix2d<T>(source.distribution());
    }
    const int nb = source.distribution()->nb();
    const std::shared_ptr<const Parallel_Blkcyclic> result_distribution
        = make_borrowed_distribution(nrglb,
                                     ncglb,
                                     nb,
                                     source.distribution());
    return Matrix2d<T>(result_distribution);
#else
    static_cast<void>(nrglb);
    static_cast<void>(ncglb);
    throw std::invalid_argument("distributed matrices require MPI support");
#endif
}
} // namespace

namespace matrix2d_linalg
{
void inv_inplace(Matrix2dView<double> matrix)
{
    assert_inverse_view(matrix);
    if (distribution_is_active(matrix))
    {
#ifdef __MPI
        matrix2d_impl::scalapack::inv(matrix);
#else
        throw std::invalid_argument("distributed inverse requires MPI support");
#endif
        return;
    }
    matrix2d_impl::lapack::inv(matrix);
}

void inv_inplace(Matrix2dView<std::complex<double>> matrix)
{
    assert_inverse_view(matrix);
    if (distribution_is_active(matrix))
    {
#ifdef __MPI
        matrix2d_impl::scalapack::inv(matrix);
#else
        throw std::invalid_argument("distributed inverse requires MPI support");
#endif
        return;
    }
    matrix2d_impl::lapack::inv(matrix);
}

void solve_inplace(Matrix2dView<double> matrix,
                   Matrix2dView<double> rhs)
{
    assert_solve_views(matrix, rhs);
    if (distribution_is_active(matrix))
    {
#ifdef __MPI
        matrix2d_impl::scalapack::solve(matrix, rhs);
#else
        throw std::invalid_argument("distributed solve requires MPI support");
#endif
        return;
    }
    matrix2d_impl::lapack::solve(matrix, rhs);
}

void solve_inplace(Matrix2dView<std::complex<double>> matrix,
                   Matrix2dView<std::complex<double>> rhs)
{
    assert_solve_views(matrix, rhs);
    if (distribution_is_active(matrix))
    {
#ifdef __MPI
        matrix2d_impl::scalapack::solve(matrix, rhs);
#else
        throw std::invalid_argument("distributed solve requires MPI support");
#endif
        return;
    }
    matrix2d_impl::lapack::solve(matrix, rhs);
}

void matmul_inplace(Matrix2dView<double> result,
                    Matrix2dView<const double> left,
                    Matrix2dView<const double> right)
{
    assert_matmul_views(result, left, right);
    if (distribution_is_active(result))
    {
#ifdef __MPI
        matrix2d_impl::scalapack::matmul(result, left, right);
#else
        throw std::invalid_argument("distributed matmul requires MPI support");
#endif
        return;
    }
    matmul_host(result, left, right);
}

void matmul_inplace(Matrix2dView<std::complex<double>> result,
                    Matrix2dView<const std::complex<double>> left,
                    Matrix2dView<const std::complex<double>> right)
{
    assert_matmul_views(result, left, right);
    if (distribution_is_active(result))
    {
#ifdef __MPI
        matrix2d_impl::scalapack::matmul(result, left, right);
#else
        throw std::invalid_argument("distributed matmul requires MPI support");
#endif
        return;
    }
    matmul_host(result, left, right);
}

void axpby_inplace(Matrix2dView<double> result,
                   const double alpha,
                   Matrix2dView<const double> left,
                   const double beta,
                   Matrix2dView<const double> right)
{
    axpby_impl(result, alpha, left, beta, right);
}

void axpby_inplace(
    Matrix2dView<std::complex<double>> result,
    const std::complex<double> alpha,
    Matrix2dView<const std::complex<double>> left,
    const std::complex<double> beta,
    Matrix2dView<const std::complex<double>> right)
{
    axpby_impl(result, alpha, left, beta, right);
}

void scale_inplace(Matrix2dView<double> matrix, const double factor)
{
    scale_impl(matrix, factor);
}

void scale_inplace(Matrix2dView<std::complex<double>> matrix,
                   const std::complex<double> factor)
{
    scale_impl(matrix, factor);
}

void set_identity_inplace(Matrix2dView<double> matrix)
{
    set_identity_impl(matrix);
}

void set_identity_inplace(Matrix2dView<std::complex<double>> matrix)
{
    set_identity_impl(matrix);
}

template <typename T>
Matrix2d<T> inv(const Matrix2d<T>& matrix)
{
    assert_supported_linalg_type<T>();
    Matrix2d<T> inverse = matrix;
    inv_inplace(inverse.view());
    return inverse;
}

template <typename T>
void inv(const Matrix2d<T>& matrix, Matrix2d<T>& inverse)
{
    Matrix2d<T> inverse_result = inv(matrix);
    inverse = std::move(inverse_result);
}

template <typename T>
Matrix2d<T> solve(const Matrix2d<T>& matrix,
                  const Matrix2d<T>& rhs)
{
    assert_supported_linalg_type<T>();
    Matrix2d<T> matrix_copy = matrix;
    Matrix2d<T> solution = rhs;
    solve_inplace(matrix_copy.view(), solution.view());
    return solution;
}

template <typename T>
void solve(const Matrix2d<T>& matrix,
           const Matrix2d<T>& rhs,
           Matrix2d<T>& solution)
{
    Matrix2d<T> solution_result = solve(matrix, rhs);
    solution = std::move(solution_result);
}

template <typename T>
Matrix2d<T> matmul(const Matrix2d<T>& left,
                   const Matrix2d<T>& right)
{
    assert_supported_linalg_type<T>();
    const int nrglb = left.nrglb();
    const int ncglb = right.ncglb();
    Matrix2d<T> result = make_matrix_like(left, nrglb, ncglb);
    matmul_inplace(result.view(), left.view(), right.view());
    return result;
}

template <typename T>
void matmul(const Matrix2d<T>& left,
            const Matrix2d<T>& right,
            Matrix2d<T>& result)
{
    Matrix2d<T> product = matmul(left, right);
    result = std::move(product);
}

template <typename T>
Matrix2d<T> axpby(const T alpha,
                  const Matrix2d<T>& left,
                  const T beta,
                  const Matrix2d<T>& right)
{
    assert_supported_linalg_type<T>();
    const int nrglb = left.nrglb();
    const int ncglb = left.ncglb();
    Matrix2d<T> result
        = make_matrix_like(left, nrglb, ncglb);
    axpby_inplace(result.view(),
                  alpha,
                  left.view(),
                  beta,
                  right.view());
    return result;
}

template <typename T>
void axpby(const T alpha,
           const Matrix2d<T>& left,
           const T beta,
           const Matrix2d<T>& right,
           Matrix2d<T>& result)
{
    Matrix2d<T> combination = axpby(alpha, left, beta, right);
    result = std::move(combination);
}

template <typename T>
Matrix2d<T> scale(const Matrix2d<T>& matrix, const T factor)
{
    assert_supported_linalg_type<T>();
    Matrix2d<T> result = matrix;
    scale_inplace(result.view(), factor);
    return result;
}

template <typename T>
void scale(const Matrix2d<T>& matrix,
           const T factor,
           Matrix2d<T>& result)
{
    Matrix2d<T> scaled = scale(matrix, factor);
    result = std::move(scaled);
}

template <typename T>
Matrix2d<T> eye(const int n)
{
    static_assert(std::is_same<T, double>::value
                      || std::is_same<T, std::complex<double>>::value,
                  "eye supports double and std::complex<double>");
    Matrix2d<T> matrix(n, n);
    set_identity_inplace(matrix.view());
    return matrix;
}

template <typename T>
Matrix2d<T> transpose(const Matrix2d<T>& matrix)
{
    return matrix.transpose();
}

template <typename T>
Matrix2d<T> dagger(const Matrix2d<T>& matrix)
{
    return matrix.dagger();
}

template Matrix2d<double> inv<double>(const Matrix2d<double>& matrix);
template void inv<double>(const Matrix2d<double>& matrix,
                          Matrix2d<double>& inverse);
template Matrix2d<std::complex<double>> inv<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix);
template void inv<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix,
    Matrix2d<std::complex<double>>& inverse);

template Matrix2d<double> solve<double>(const Matrix2d<double>& matrix,
                                        const Matrix2d<double>& rhs);
template void solve<double>(const Matrix2d<double>& matrix,
                            const Matrix2d<double>& rhs,
                            Matrix2d<double>& solution);
template Matrix2d<std::complex<double>> solve<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix,
    const Matrix2d<std::complex<double>>& rhs);
template void solve<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix,
    const Matrix2d<std::complex<double>>& rhs,
    Matrix2d<std::complex<double>>& solution);

template Matrix2d<double> matmul<double>(const Matrix2d<double>& left,
                                         const Matrix2d<double>& right);
template void matmul<double>(const Matrix2d<double>& left,
                             const Matrix2d<double>& right,
                             Matrix2d<double>& result);
template Matrix2d<std::complex<double>> matmul<std::complex<double>>(
    const Matrix2d<std::complex<double>>& left,
    const Matrix2d<std::complex<double>>& right);
template void matmul<std::complex<double>>(
    const Matrix2d<std::complex<double>>& left,
    const Matrix2d<std::complex<double>>& right,
    Matrix2d<std::complex<double>>& result);

template Matrix2d<double> axpby<double>(const double alpha,
                                        const Matrix2d<double>& left,
                                        const double beta,
                                        const Matrix2d<double>& right);
template void axpby<double>(const double alpha,
                            const Matrix2d<double>& left,
                            const double beta,
                            const Matrix2d<double>& right,
                            Matrix2d<double>& result);
template Matrix2d<std::complex<double>> axpby<std::complex<double>>(
    const std::complex<double> alpha,
    const Matrix2d<std::complex<double>>& left,
    const std::complex<double> beta,
    const Matrix2d<std::complex<double>>& right);
template void axpby<std::complex<double>>(
    const std::complex<double> alpha,
    const Matrix2d<std::complex<double>>& left,
    const std::complex<double> beta,
    const Matrix2d<std::complex<double>>& right,
    Matrix2d<std::complex<double>>& result);

template Matrix2d<double> scale<double>(const Matrix2d<double>& matrix,
                                        const double factor);
template void scale<double>(const Matrix2d<double>& matrix,
                            const double factor,
                            Matrix2d<double>& result);
template Matrix2d<std::complex<double>> scale<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix,
    const std::complex<double> factor);
template void scale<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix,
    const std::complex<double> factor,
    Matrix2d<std::complex<double>>& result);

template Matrix2d<double> eye<double>(const int n);
template Matrix2d<std::complex<double>> eye<std::complex<double>>(const int n);
template Matrix2d<double> transpose<double>(const Matrix2d<double>& matrix);
template Matrix2d<std::complex<double>> transpose<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix);
template Matrix2d<double> dagger<double>(const Matrix2d<double>& matrix);
template Matrix2d<std::complex<double>> dagger<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix);
} // namespace matrix2d_linalg
