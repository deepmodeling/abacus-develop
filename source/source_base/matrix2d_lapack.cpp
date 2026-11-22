#include "source_base/matrix2d_lapack.h"

#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/matrix2d_view.h"

#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
void assert_host_layout(const Matrix2dView<double>& matrix)
{
    if (matrix.is_active())
    {
        throw std::invalid_argument("LAPACK requires a host matrix view");
    }
}

void assert_host_layout(const Matrix2dView<std::complex<double>>& matrix)
{
    if (matrix.is_active())
    {
        throw std::invalid_argument("LAPACK requires a host matrix view");
    }
}

void assert_inverse_dimensions(const int nr, const int nc)
{
    if (nr <= 0 || nr != nc)
    {
        throw std::invalid_argument("matrix inverse requires a non-empty square matrix");
    }
}

void assert_solve_dimensions(const int nr,
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
    assert_host_layout(matrix);
    const int n = matrix.nrglb();
    assert_inverse_dimensions(n, matrix.ncglb());
    const std::size_t n_ = static_cast<std::size_t>(n);

    const int lda = matrix.lld();
    std::vector<int> ipiv(n_);
    int info = 0;
    container::lapackConnector::getrf(
        /* m= */ n,
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* info= */ info);
    throw_lapack_error("LAPACK getrf", info);

    int lws = -1;
    std::vector<double> ws(1);
    container::lapackConnector::getri(
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* work= */ ws.data(),
        /* lwork= */ lws,
        /* info= */ info);
    throw_lapack_error("LAPACK getri workspace query", info);
    lws = static_cast<int>(ws[0]);
    const std::size_t lws_ = static_cast<std::size_t>(lws);
    ws.resize(lws_);
    container::lapackConnector::getri(
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* work= */ ws.data(),
        /* lwork= */ lws,
        /* info= */ info);
    throw_lapack_error("LAPACK getri", info);
}

void inv(Matrix2dView<std::complex<double>> matrix)
{
    assert_host_layout(matrix);
    const int n = matrix.nrglb();
    assert_inverse_dimensions(n, matrix.ncglb());
    const std::size_t n_ = static_cast<std::size_t>(n);

    const int lda = matrix.lld();
    std::vector<int> ipiv(n_);
    int info = 0;
    container::lapackConnector::getrf(
        /* m= */ n,
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* info= */ info);
    throw_lapack_error("LAPACK getrf", info);

    int lws = -1;
    std::vector<std::complex<double>> ws(1);
    container::lapackConnector::getri(
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* work= */ ws.data(),
        /* lwork= */ lws,
        /* info= */ info);
    throw_lapack_error("LAPACK getri workspace query", info);
    lws = static_cast<int>(ws[0].real());
    const std::size_t lws_ = static_cast<std::size_t>(lws);
    ws.resize(lws_);
    container::lapackConnector::getri(
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* work= */ ws.data(),
        /* lwork= */ lws,
        /* info= */ info);
    throw_lapack_error("LAPACK getri", info);
}

void solve(Matrix2dView<double> matrix,
           Matrix2dView<double> rhs)
{
    assert_host_layout(matrix);
    assert_host_layout(rhs);
    const int n = matrix.nrglb();
    const int nrhs = rhs.ncglb();
    assert_solve_dimensions(n, matrix.ncglb(), rhs.nrglb());
    const std::size_t n_ = static_cast<std::size_t>(n);
    const std::size_t nrhs_ = static_cast<std::size_t>(nrhs);

    const int lda = matrix.lld();
    std::vector<int> ipiv(n_);
    int info = 0;
    container::lapackConnector::getrf(
        /* m= */ n,
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* info= */ info);
    throw_lapack_error("LAPACK getrf", info);

    // Column-major LAPACK right-hand side for A * X = B.
    std::vector<double> b(n_ * nrhs_);
    for (int i = 0; i < n; ++i)
    {
        const std::size_t i_ = static_cast<std::size_t>(i);
        for (int j = 0; j < nrhs; ++j)
        {
            const std::size_t j_ = static_cast<std::size_t>(j);
            b[i_ + j_ * n_] = rhs.atloc(i, j);
        }
    }

    const char trans = matrix.has_distribution() ? 'N' : 'T';
    container::lapackConnector::getrs(
        /* trans= */ trans,
        /* n= */ n,
        /* nrhs= */ nrhs,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* b= */ b.data(),
        /* ldb= */ n,
        /* info= */ info);
    throw_lapack_error("LAPACK getrs", info);

    for (int i = 0; i < n; ++i)
    {
        const std::size_t i_ = static_cast<std::size_t>(i);
        for (int j = 0; j < nrhs; ++j)
        {
            const std::size_t j_ = static_cast<std::size_t>(j);
            rhs.atloc(i, j) = b[i_ + j_ * n_];
        }
    }
}

void solve(Matrix2dView<std::complex<double>> matrix,
           Matrix2dView<std::complex<double>> rhs)
{
    assert_host_layout(matrix);
    assert_host_layout(rhs);
    const int n = matrix.nrglb();
    const int nrhs = rhs.ncglb();
    assert_solve_dimensions(n, matrix.ncglb(), rhs.nrglb());
    const std::size_t n_ = static_cast<std::size_t>(n);
    const std::size_t nrhs_ = static_cast<std::size_t>(nrhs);

    const int lda = matrix.lld();
    std::vector<int> ipiv(n_);
    int info = 0;
    container::lapackConnector::getrf(
        /* m= */ n,
        /* n= */ n,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* info= */ info);
    throw_lapack_error("LAPACK getrf", info);

    // Column-major LAPACK right-hand side for A * X = B.
    std::vector<std::complex<double>> b(n_ * nrhs_);
    for (int i = 0; i < n; ++i)
    {
        const std::size_t i_ = static_cast<std::size_t>(i);
        for (int j = 0; j < nrhs; ++j)
        {
            const std::size_t j_ = static_cast<std::size_t>(j);
            b[i_ + j_ * n_] = rhs.atloc(i, j);
        }
    }

    const char trans = matrix.has_distribution() ? 'N' : 'T';
    container::lapackConnector::getrs(
        /* trans= */ trans,
        /* n= */ n,
        /* nrhs= */ nrhs,
        /* a= */ matrix.data(),
        /* lda= */ lda,
        /* ipiv= */ ipiv.data(),
        /* b= */ b.data(),
        /* ldb= */ n,
        /* info= */ info);
    throw_lapack_error("LAPACK getrs", info);

    for (int i = 0; i < n; ++i)
    {
        const std::size_t i_ = static_cast<std::size_t>(i);
        for (int j = 0; j < nrhs; ++j)
        {
            const std::size_t j_ = static_cast<std::size_t>(j);
            rhs.atloc(i, j) = b[i_ + j_ * n_];
        }
    }
}
} // namespace lapack
} // namespace matrix2d_impl
