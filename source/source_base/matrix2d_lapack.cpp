#include "source_base/matrix2d_lapack.h"

#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/matrix2d_view.h"

#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void check_host_layout(const Matrix2dView<double>& matrix)
{
    if (matrix.is_distributed())
    {
        throw std::invalid_argument("LAPACK requires a host matrix view");
    }
}

void check_host_layout(const Matrix2dView<std::complex<double>>& matrix)
{
    if (matrix.is_distributed())
    {
        throw std::invalid_argument("LAPACK requires a host matrix view");
    }
}

void check_inverse_dimensions(const int nr, const int nc)
{
    if (nr <= 0 || nr != nc)
    {
        throw std::invalid_argument("matrix inverse requires a non-empty square matrix");
    }
}

void check_solve_dimensions(const int nr,
                            const int nc,
                            const int rhs_nr)
{
    if (nr <= 0 || nr != nc || rhs_nr != nr)
    {
        throw std::invalid_argument(
            "linear solve requires a square matrix and matching right-hand side");
    }
}

void throw_lapack_error(const char* operation, const int info)
{
    if (info < 0)
    {
        throw std::invalid_argument(std::string(operation)
                                    + " received an invalid argument");
    }
    if (info > 0)
    {
        throw std::runtime_error(std::string(operation)
                                 + " failed because the matrix is singular");
    }
}
} // namespace

namespace matrix2d_impl
{
namespace lapack
{
void inv(Matrix2dView<double> matrix)
{
    check_host_layout(matrix);
    const int n = matrix.nrglb();
    check_inverse_dimensions(n, matrix.ncglb());

    const int lda = matrix.lld();
    std::vector<int> ipiv(static_cast<std::size_t>(n));
    int info = 0;
    container::lapackConnector::getrf(
        /* m= */ n,
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* info= */ info);
    throw_lapack_error("LAPACK getrf", info);

    int lwork = -1;
    std::vector<double> work(1);
    container::lapackConnector::getri(
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* work= */ work.data(),
        /* lwork= */ lwork,
        /* info= */ info);
    throw_lapack_error("LAPACK getri workspace query", info);
    lwork = static_cast<int>(work[0]);
    work.resize(static_cast<std::size_t>(lwork));
    container::lapackConnector::getri(
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* work= */ work.data(),
        /* lwork= */ lwork,
        /* info= */ info);
    throw_lapack_error("LAPACK getri", info);
}

void inv(Matrix2dView<std::complex<double>> matrix)
{
    check_host_layout(matrix);
    const int n = matrix.nrglb();
    check_inverse_dimensions(n, matrix.ncglb());

    const int lda = matrix.lld();
    std::vector<int> ipiv(static_cast<std::size_t>(n));
    int info = 0;
    container::lapackConnector::getrf(
        /* m= */ n,
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* info= */ info);
    throw_lapack_error("LAPACK getrf", info);

    int lwork = -1;
    std::vector<std::complex<double>> work(1);
    container::lapackConnector::getri(
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* work= */ work.data(),
        /* lwork= */ lwork,
        /* info= */ info);
    throw_lapack_error("LAPACK getri workspace query", info);
    lwork = static_cast<int>(work[0].real());
    work.resize(static_cast<std::size_t>(lwork));
    container::lapackConnector::getri(
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* work= */ work.data(),
        /* lwork= */ lwork,
        /* info= */ info);
    throw_lapack_error("LAPACK getri", info);
}

void solve(Matrix2dView<double> matrix,
           Matrix2dView<double> rhs)
{
    check_host_layout(matrix);
    check_host_layout(rhs);
    const int n = matrix.nrglb();
    const int nrhs = rhs.ncglb();
    check_solve_dimensions(n, matrix.ncglb(), rhs.nrglb());

    const int lda = matrix.lld();
    std::vector<int> ipiv(static_cast<std::size_t>(n));
    int info = 0;
    container::lapackConnector::getrf(
        /* m= */ n,
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* info= */ info);
    throw_lapack_error("LAPACK getrf", info);

    std::vector<double> rhs_column_major(
        static_cast<std::size_t>(n) * static_cast<std::size_t>(nrhs));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < nrhs; ++j)
        {
            rhs_column_major[static_cast<std::size_t>(i)
                             + static_cast<std::size_t>(j)
                                 * static_cast<std::size_t>(n)]
                = rhs.at_local(i, j);
        }
    }

    const char trans = 'T';
    container::lapackConnector::getrs(
        /* trans= */ trans,
        /* n= */ n,
        /* nrhs= */ nrhs,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* b= */ rhs_column_major.data(),
        /* ldb= */ n,
        /* info= */ info);
    throw_lapack_error("LAPACK getrs", info);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < nrhs; ++j)
        {
            rhs.at_local(i, j)
                = rhs_column_major[static_cast<std::size_t>(i)
                                   + static_cast<std::size_t>(j)
                                       * static_cast<std::size_t>(n)];
        }
    }
}

void solve(Matrix2dView<std::complex<double>> matrix,
           Matrix2dView<std::complex<double>> rhs)
{
    check_host_layout(matrix);
    check_host_layout(rhs);
    const int n = matrix.nrglb();
    const int nrhs = rhs.ncglb();
    check_solve_dimensions(n, matrix.ncglb(), rhs.nrglb());

    const int lda = matrix.lld();
    std::vector<int> ipiv(static_cast<std::size_t>(n));
    int info = 0;
    container::lapackConnector::getrf(
        /* m= */ n,
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* info= */ info);
    throw_lapack_error("LAPACK getrf", info);

    std::vector<std::complex<double>> rhs_column_major(
        static_cast<std::size_t>(n) * static_cast<std::size_t>(nrhs));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < nrhs; ++j)
        {
            rhs_column_major[static_cast<std::size_t>(i)
                             + static_cast<std::size_t>(j)
                                 * static_cast<std::size_t>(n)]
                = rhs.at_local(i, j);
        }
    }

    const char trans = 'T';
    container::lapackConnector::getrs(
        /* trans= */ trans,
        /* n= */ n,
        /* nrhs= */ nrhs,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* b= */ rhs_column_major.data(),
        /* ldb= */ n,
        /* info= */ info);
    throw_lapack_error("LAPACK getrs", info);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < nrhs; ++j)
        {
            rhs.at_local(i, j)
                = rhs_column_major[static_cast<std::size_t>(i)
                                   + static_cast<std::size_t>(j)
                                       * static_cast<std::size_t>(n)];
        }
    }
}
} // namespace lapack
} // namespace matrix2d_impl
