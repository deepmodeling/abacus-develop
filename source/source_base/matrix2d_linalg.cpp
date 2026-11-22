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
} // namespace

namespace matrix2d_linalg
{
void inv(Matrix2dView<double> matrix)
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

void inv(Matrix2dView<std::complex<double>> matrix)
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

void solve(Matrix2dView<double> matrix,
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

void solve(Matrix2dView<std::complex<double>> matrix,
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

void matmul(Matrix2dView<double> result,
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

void matmul(Matrix2dView<std::complex<double>> result,
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

void axpby(Matrix2dView<double> result,
           const double alpha,
           Matrix2dView<const double> left,
           const double beta,
           Matrix2dView<const double> right)
{
    axpby_impl(result, alpha, left, beta, right);
}

void axpby(
    Matrix2dView<std::complex<double>> result,
    const std::complex<double> alpha,
    Matrix2dView<const std::complex<double>> left,
    const std::complex<double> beta,
    Matrix2dView<const std::complex<double>> right)
{
    axpby_impl(result, alpha, left, beta, right);
}

void scale(Matrix2dView<double> matrix, const double factor)
{
    scale_impl(matrix, factor);
}

void scale(Matrix2dView<std::complex<double>> matrix,
           const std::complex<double> factor)
{
    scale_impl(matrix, factor);
}

void set_identity(Matrix2dView<double> matrix)
{
    set_identity_impl(matrix);
}

void set_identity(Matrix2dView<std::complex<double>> matrix)
{
    set_identity_impl(matrix);
}

template <typename T>
Matrix2d<T> eye(const int n)
{
    static_assert(std::is_same<T, double>::value
                      || std::is_same<T, std::complex<double>>::value,
                  "eye supports double and std::complex<double>");
    Matrix2d<T> matrix(n, n);
    set_identity_impl(matrix.view());
    return matrix;
}

template Matrix2d<double> eye<double>(const int n);
template Matrix2d<std::complex<double>> eye<std::complex<double>>(const int n);
} // namespace matrix2d_linalg
