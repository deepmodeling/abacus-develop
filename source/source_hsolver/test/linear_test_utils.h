#ifndef HSOLVER_LINEAR_TEST_UTILS_H
#define HSOLVER_LINEAR_TEST_UTILS_H

#include "source_base/module_container/base/third_party/lapack.h"
#include "source_base/parallel_comm.h"
#include "source_hsolver/diag_comm_info.h"

#include <complex>
#include <stdexcept>
#include <vector>

namespace linear_test
{
using Complex = std::complex<double>;

inline hsolver::diag_comm_info world_comm()
{
#ifdef __MPI
    MPICommGroup group(MPI_COMM_WORLD);
    return hsolver::diag_comm_info(MPI_COMM_WORLD, group.grank, group.gsize);
#else
    return hsolver::diag_comm_info(0, 1);
#endif
}

// Column-major matrices with a decoupled first eigenvector and a dense complex block.
template <typename T>
std::vector<T> hermitian(const int n)
{
    std::vector<T> h(n * n, T(0));
    for (int j = 0; j < n; ++j)
    {
        h[j * n + j] = T(1.0 + 0.4 * j);
        for (int i = 1; i < j; ++i)
        {
            const T value(0.12 / (j - i + 1), 0.07 / (i + j + 1));
            h[j * n + i] = value;
            h[i * n + j] = std::conj(value);
        }
    }
    return h;
}

template <typename T>
std::vector<T> multiply(const std::vector<T>& a, const std::vector<T>& x, const int n, const int nvec)
{
    std::vector<T> b(n * nvec, T(0));
    for (int band = 0; band < nvec; ++band)
    {
        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                b[band * n + i] += a[j * n + i] * x[band * n + j];
            }
        }
    }
    return b;
}

// Extract the local rows of column-major reference data, leaving padding zero.
template <typename T>
std::vector<T> local_columns(const std::vector<T>& global, int n, int columns, int start, int dim, int ld)
{
    std::vector<T> local(ld * columns, T(0));
    for (int j = 0; j < columns; ++j)
    {
        for (int i = 0; i < dim; ++i)
        {
            local[j * ld + i] = global[j * n + start + i];
        }
    }
    return local;
}

// Solve in double precision even when the iterative input was rounded to float.
template <typename T>
std::vector<Complex> lapack_solve(const std::vector<T>& a, const std::vector<T>& b, const int n, const int nvec)
{
    std::vector<Complex> lu(a.begin(), a.end());
    std::vector<Complex> x(b.begin(), b.end());
    std::vector<int> pivots(n);
    int info = 0;
    container::lapackConnector::getrf(n, n, lu.data(), n, pivots.data(), info);
    if (info != 0)
    {
        throw std::runtime_error("LAPACK reference factorization failed");
    }
    container::lapackConnector::getrs('N', n, nvec, lu.data(), n, pivots.data(), x.data(), n, info);
    if (info != 0)
    {
        throw std::runtime_error("LAPACK reference solve failed");
    }
    return x;
}
} // namespace linear_test
#endif
