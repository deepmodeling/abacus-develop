#include "source_base/matrix2d_linalg.h"

#include "source_base/matrix2d_lapack.h"
#include "source_base/matrix2d_view.h"
#include "source_base/parallel_blkcyc.h"
#ifdef __MPI
#include "source_base/matrix2d_scalapack.h"
#endif

#include <complex>
#include <stdexcept>

namespace
{
template <typename T>
bool has_active_distribution(const Matrix2dView<T>& matrix)
{
    return matrix.is_distributed() && !matrix.distribution()->is_serial();
}

template <typename T>
bool has_serial_distribution(const Matrix2dView<T>& matrix)
{
    return matrix.is_distributed() && matrix.distribution()->is_serial();
}

template <typename T, typename U>
void check_same_shape(const Matrix2dView<T>& left,
                      const Matrix2dView<U>& right)
{
    if (left.nrglb() != right.nrglb() || left.ncglb() != right.ncglb())
    {
        throw std::invalid_argument("matrix shapes are incompatible");
    }
}

template <typename T, typename U>
void check_compatible_layout(const Matrix2dView<T>& left,
                             const Matrix2dView<U>& right)
{
    if (left.is_distributed() != right.is_distributed())
    {
        throw std::invalid_argument(
            "matrix operations require matching host or distributed layouts");
    }
    if (!left.is_distributed())
    {
        return;
    }
    const bool left_is_serial = has_serial_distribution(left);
    const bool right_is_serial = has_serial_distribution(right);
    if (left_is_serial != right_is_serial)
    {
        throw std::invalid_argument(
            "matrix operations require matching serial or parallel distributions");
    }
    if (left_is_serial)
    {
        return;
    }
    if (left.desc() == nullptr
        || right.desc() == nullptr
        || left.desc()[1] != right.desc()[1]
        || left.desc()[4] != right.desc()[4]
        || left.desc()[5] != right.desc()[5]
        || left.desc()[6] != right.desc()[6]
        || left.desc()[7] != right.desc()[7])
    {
        throw std::invalid_argument(
            "matrix operations require compatible distributed layouts");
    }
}

template <typename T>
void check_matmul_views(Matrix2dView<T> result,
                        Matrix2dView<const T> left,
                        Matrix2dView<const T> right)
{
    if (left.ncglb() != right.nrglb()
        || result.nrglb() != left.nrglb()
        || result.ncglb() != right.ncglb())
    {
        throw std::invalid_argument("matrix multiplication dimensions are incompatible");
    }
    check_compatible_layout(result, left);
    check_compatible_layout(result, right);
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
                value += left.at_local(i, k) * right.at_local(k, j);
            }
            result.at_local(i, j) = value;
        }
    }
}

template <typename T>
void linear_combination_impl(Matrix2dView<T> result,
                             const T alpha,
                             Matrix2dView<const T> left,
                             const T beta,
                             Matrix2dView<const T> right)
{
    check_same_shape(result, left);
    check_same_shape(result, right);
    check_compatible_layout(result, left);
    check_compatible_layout(result, right);
    for (int i = 0; i < result.nrloc(); ++i)
    {
        for (int j = 0; j < result.ncloc(); ++j)
        {
            result.at_local(i, j)
                = alpha * left.at_local(i, j) + beta * right.at_local(i, j);
        }
    }
}

template <typename T>
void scale_impl(Matrix2dView<T> matrix, const T factor)
{
    for (int i = 0; i < matrix.nrloc(); ++i)
    {
        for (int j = 0; j < matrix.ncloc(); ++j)
        {
            matrix.at_local(i, j) *= factor;
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
    const std::shared_ptr<const Parallel_Blkcyclic> dist
        = matrix.distribution();
    for (int i = 0; i < matrix.nrloc(); ++i)
    {
        for (int j = 0; j < matrix.ncloc(); ++j)
        {
            const int irglb = dist == nullptr ? i : dist->irglb(i);
            const int icglb = dist == nullptr ? j : dist->icglb(j);
            matrix.at_local(i, j) = irglb == icglb ? T(1.0) : T();
        }
    }
}

void check_inverse_view(const Matrix2dView<double>& matrix)
{
    if (matrix.nrglb() <= 0 || matrix.nrglb() != matrix.ncglb())
    {
        throw std::invalid_argument("matrix inverse requires a non-empty square matrix");
    }
}

void check_inverse_view(const Matrix2dView<std::complex<double>>& matrix)
{
    if (matrix.nrglb() <= 0 || matrix.nrglb() != matrix.ncglb())
    {
        throw std::invalid_argument("matrix inverse requires a non-empty square matrix");
    }
}

void check_solve_views(const Matrix2dView<double>& matrix,
                       const Matrix2dView<double>& rhs)
{
    if (matrix.nrglb() <= 0
        || matrix.nrglb() != matrix.ncglb()
        || rhs.nrglb() != matrix.nrglb())
    {
        throw std::invalid_argument(
            "linear solve requires a square matrix and matching right-hand side");
    }
    if (matrix.is_distributed() != rhs.is_distributed())
    {
        throw std::invalid_argument(
            "linear solve requires matching host or distributed layouts");
    }
    if (has_serial_distribution(matrix) != has_serial_distribution(rhs))
    {
        throw std::invalid_argument(
            "linear solve requires matching serial or parallel distributions");
    }
    if (has_serial_distribution(matrix))
    {
        return;
    }
    if (matrix.is_distributed()
        && (matrix.desc() == nullptr
            || rhs.desc() == nullptr
            || matrix.desc()[1] != rhs.desc()[1]
            || matrix.desc()[4] != rhs.desc()[4]
            || matrix.desc()[5] != rhs.desc()[5]
            || matrix.desc()[6] != rhs.desc()[6]
            || matrix.desc()[7] != rhs.desc()[7]))
    {
        throw std::invalid_argument(
            "linear solve requires compatible distributed layouts");
    }
}

void check_solve_views(const Matrix2dView<std::complex<double>>& matrix,
                       const Matrix2dView<std::complex<double>>& rhs)
{
    if (matrix.nrglb() <= 0
        || matrix.nrglb() != matrix.ncglb()
        || rhs.nrglb() != matrix.nrglb())
    {
        throw std::invalid_argument(
            "linear solve requires a square matrix and matching right-hand side");
    }
    if (matrix.is_distributed() != rhs.is_distributed())
    {
        throw std::invalid_argument(
            "linear solve requires matching host or distributed layouts");
    }
    if (has_serial_distribution(matrix) != has_serial_distribution(rhs))
    {
        throw std::invalid_argument(
            "linear solve requires matching serial or parallel distributions");
    }
    if (has_serial_distribution(matrix))
    {
        return;
    }
    if (matrix.is_distributed()
        && (matrix.desc() == nullptr
            || rhs.desc() == nullptr
            || matrix.desc()[1] != rhs.desc()[1]
            || matrix.desc()[4] != rhs.desc()[4]
            || matrix.desc()[5] != rhs.desc()[5]
            || matrix.desc()[6] != rhs.desc()[6]
            || matrix.desc()[7] != rhs.desc()[7]))
    {
        throw std::invalid_argument(
            "linear solve requires compatible distributed layouts");
    }
}
} // namespace

namespace matrix2d_linalg
{
void inv(Matrix2dView<double> matrix)
{
    check_inverse_view(matrix);
    if (has_active_distribution(matrix))
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
    check_inverse_view(matrix);
    if (has_active_distribution(matrix))
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
    check_solve_views(matrix, rhs);
    if (has_active_distribution(matrix))
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
    check_solve_views(matrix, rhs);
    if (has_active_distribution(matrix))
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
    check_matmul_views(result, left, right);
    if (has_active_distribution(result))
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
    check_matmul_views(result, left, right);
    if (has_active_distribution(result))
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

void linear_combination(Matrix2dView<double> result,
                        const double alpha,
                        Matrix2dView<const double> left,
                        const double beta,
                        Matrix2dView<const double> right)
{
    linear_combination_impl(result, alpha, left, beta, right);
}

void linear_combination(
    Matrix2dView<std::complex<double>> result,
    const std::complex<double> alpha,
    Matrix2dView<const std::complex<double>> left,
    const std::complex<double> beta,
    Matrix2dView<const std::complex<double>> right)
{
    linear_combination_impl(result, alpha, left, beta, right);
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
} // namespace matrix2d_linalg
