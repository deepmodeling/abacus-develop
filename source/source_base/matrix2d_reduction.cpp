#include "source_base/matrix2d_reduction.h"

#include "source_base/parallel_common.h"
#ifdef __MPI
#include "source_base/matrix2d_mpi.h"
#endif

#include <algorithm>
#include <cmath>
#include <cstring>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
// Width in bytes of a native int in the homogeneous diagonal wire format.
constexpr int W_INT = sizeof(int);

#ifdef __MPI
namespace mpi = matrix2d_impl::mpi;
#endif

namespace elementwise
{
template <typename T>
T conj(const T& value)
{
    return value;
}

template <typename T>
std::complex<T> conj(const std::complex<T>& value)
{
    return std::conj(value);
}
} // namespace elementwise

template <typename T>
double absolute_square(const T& value)
{
    const double real_value = static_cast<double>(value);
    return real_value * real_value;
}

template <typename T>
double absolute_square(const std::complex<T>& value)
{
    const double real = static_cast<double>(value.real());
    const double imag = static_cast<double>(value.imag());
    return real * real + imag * imag;
}

template <typename T>
double absolute_value(const T& value)
{
    const double widened_value = static_cast<double>(value);
    return std::abs(widened_value);
}

template <typename T>
double absolute_value(const std::complex<T>& value)
{
    return static_cast<double>(std::abs(value));
}

template <typename Value, typename Element>
void reduce_sum_over_processes(Value& value, const Matrix2d<Element>& matrix)
{
#ifdef __MPI
    const std::shared_ptr<const Parallel_Blkcyclic> distrib
        = matrix.distribution();
    if (distrib != nullptr && !distrib->is_serial())
    {
        const MPI_Comm comm = distrib->comm();
        if (comm != MPI_COMM_NULL)
        {
            Parallel_Common::allreduce_sum(value, comm);
        }
    }
#else
    static_cast<void>(value);
    static_cast<void>(matrix);
#endif
}

template <typename T>
void reduce_max_over_processes(double& value, const Matrix2d<T>& matrix)
{
#ifdef __MPI
    const std::shared_ptr<const Parallel_Blkcyclic> distrib
        = matrix.distribution();
    if (distrib != nullptr && !distrib->is_serial())
    {
        const MPI_Comm comm = distrib->comm();
        if (comm != MPI_COMM_NULL)
        {
            Parallel_Common::allreduce_max(value, comm);
        }
    }
#else
    static_cast<void>(value);
    static_cast<void>(matrix);
#endif
}

} // namespace

namespace matrix2d_impl
{
template <typename T>
T sum(const Matrix2d<T>& matrix)
{
    T value = T();
    const int nrloc = matrix.nrloc();
    const int ncloc = matrix.ncloc();
#ifdef _OPENMP
    std::vector<T> partial_sums(omp_get_max_threads(), T());
#pragma omp parallel
    {
        const int ip = omp_get_thread_num();
#pragma omp for collapse(2) schedule(static)
        for (int i = 0; i < nrloc; ++i)
        {
            for (int j = 0; j < ncloc; ++j)
            {
                partial_sums[ip] += matrix.atloc(i, j);
            }
        }
    }
    for (const T& partial_sum : partial_sums)
    {
        value += partial_sum;
    }
#else
    for (int i = 0; i < nrloc; ++i)
    {
        for (int j = 0; j < ncloc; ++j)
        {
            value += matrix.atloc(i, j);
        }
    }
#endif
    reduce_sum_over_processes(value, matrix);
    return value;
}

template <typename T>
double norm(const Matrix2d<T>& matrix)
{
    double value = 0.0;
    const int nrloc = matrix.nrloc();
    const int ncloc = matrix.ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(+ : value) schedule(static)
#endif
    for (int i = 0; i < nrloc; ++i)
    {
        for (int j = 0; j < ncloc; ++j)
        {
            value += absolute_square(matrix.atloc(i, j));
        }
    }
    reduce_sum_over_processes(value, matrix);
    return std::sqrt(value);
}

template <typename T>
T dot(const Matrix2d<T>& left, const Matrix2d<T>& right)
{
    T value = T();
    const int nrloc = left.nrloc();
    const int ncloc = left.ncloc();
#ifdef _OPENMP
    std::vector<T> partial_sums(omp_get_max_threads(), T());
#pragma omp parallel
    {
        const int ip = omp_get_thread_num();
#pragma omp for collapse(2) schedule(static)
        for (int i = 0; i < nrloc; ++i)
        {
            for (int j = 0; j < ncloc; ++j)
            {
                partial_sums[ip] += elementwise::conj(left.atloc(i, j))
                    * right.atloc(i, j);
            }
        }
    }
    for (const T& partial_sum : partial_sums)
    {
        value += partial_sum;
    }
#else
    for (int i = 0; i < nrloc; ++i)
    {
        for (int j = 0; j < ncloc; ++j)
        {
            value += elementwise::conj(left.atloc(i, j))
                * right.atloc(i, j);
        }
    }
#endif
    reduce_sum_over_processes(value, left);
    return value;
}

template <typename T>
double max_abs(const Matrix2d<T>& matrix)
{
    double max_value = 0.0;
    const int nrloc = matrix.nrloc();
    const int ncloc = matrix.ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(max : max_value) schedule(static)
#endif
    for (int i = 0; i < nrloc; ++i)
    {
        for (int j = 0; j < ncloc; ++j)
        {
            const T& element = matrix.atloc(i, j);
            const double element_abs = absolute_value(element);
            max_value = std::max(max_value, element_abs);
        }
    }
    reduce_max_over_processes(max_value, matrix);
    return max_value;
}

template <typename T>
std::vector<T> diagonal(const Matrix2d<T>& matrix)
{
    const int diagonal_size = std::min(matrix.nrglb(), matrix.ncglb());
    std::vector<T> values(static_cast<std::size_t>(diagonal_size), T());
    const int nrloc = matrix.nrloc();
    const int ncloc = matrix.ncloc();
    const std::shared_ptr<const Parallel_Blkcyclic> distrib
        = matrix.distribution();
    std::vector<std::pair<int, T>> local_entries;
    for (int irloc = 0; irloc < nrloc; ++irloc)
    {
        for (int icloc = 0; icloc < ncloc; ++icloc)
        {
            const int irglb = distrib == nullptr ? irloc : distrib->irglb(irloc);
            const int icglb = distrib == nullptr ? icloc : distrib->icglb(icloc);
            if (irglb == icglb && irglb < diagonal_size)
            {
                const T value = matrix.atloc(irloc, icloc);
                local_entries.push_back(std::pair<int, T>(irglb, value));
            }
        }
    }
#ifdef __MPI
    if (distrib != nullptr && !distrib->is_serial())
    {
        const MPI_Comm comm = distrib->comm();
        if (comm != MPI_COMM_NULL)
        {
            const std::size_t nbytes_entry = W_INT + sizeof(T);
            const std::size_t nbytes_local
                = W_INT + local_entries.size() * nbytes_entry;
            std::vector<unsigned char> local_diagonal_buffer(nbytes_local);
            const int nentry_local
                = mpi::checked_mpi_count(local_entries.size());
            std::memcpy(local_diagonal_buffer.data(),
                        &nentry_local,
                        W_INT);
            std::size_t offset = W_INT;
            for (const std::pair<int, T>& entry : local_entries)
            {
                std::memcpy(local_diagonal_buffer.data() + offset,
                            &entry.first,
                            sizeof(entry.first));
                offset += sizeof(entry.first);
                std::memcpy(local_diagonal_buffer.data() + offset,
                            &entry.second,
                            sizeof(entry.second));
                offset += sizeof(entry.second);
            }

            const int nbytes_local_mpi
                = mpi::checked_mpi_count(nbytes_local);
            const std::vector<int> nbytes_recv
                = mpi::allgather_nbytes(nbytes_local_mpi, comm);
            const std::vector<int> rdispls
                = mpi::make_byte_offsets(nbytes_recv);
            const int total_bytes
                = rdispls.back() + nbytes_recv.back();
            std::vector<unsigned char> received(
                static_cast<std::size_t>(total_bytes));
            Parallel_Common::allgatherv_bytes(local_diagonal_buffer.data(),
                                              nbytes_local_mpi,
                                              received.data(),
                                              nbytes_recv.data(),
                                              rdispls.data(),
                                              comm);

            for (std::size_t ip = 0; ip < nbytes_recv.size(); ++ip)
            {
                const std::size_t begin
                    = static_cast<std::size_t>(rdispls[ip]);
                const std::size_t end
                    = begin + static_cast<std::size_t>(nbytes_recv[ip]);
                if (end - begin < W_INT)
                {
                    throw std::runtime_error(
                        "Matrix2d diagonal data buffer is truncated");
                }
                int nentry = 0;
                std::memcpy(&nentry,
                            received.data() + begin,
                            sizeof(nentry));
                const std::size_t nbytes_expected
                    = static_cast<std::size_t>(nentry) * nbytes_entry;
                if (nentry < 0
                    || end - begin != W_INT + nbytes_expected)
                {
                    throw std::runtime_error(
                        "Matrix2d diagonal data buffer has an invalid size");
                }
                std::size_t entry_offset = begin + sizeof(nentry);
                for (int ientry = 0; ientry < nentry; ++ientry)
                {
                    int irglb = 0;
                    T value = T();
                    std::memcpy(&irglb,
                                received.data() + entry_offset,
                                W_INT);
                    entry_offset += W_INT;
                    std::memcpy(&value,
                                received.data() + entry_offset,
                                sizeof(value));
                    entry_offset += sizeof(value);
                    if (irglb < 0 || irglb >= diagonal_size)
                    {
                        throw std::runtime_error(
                            "Matrix2d diagonal entry is outside global shape");
                    }
                    values[static_cast<std::size_t>(irglb)] += value;
                }
            }
        }
    }
#endif
    if (distrib == nullptr || distrib->is_serial())
    {
        for (const std::pair<int, T>& entry : local_entries)
        {
            values[static_cast<std::size_t>(entry.first)] += entry.second;
        }
    }
    return values;
}

template int sum<int>(const Matrix2d<int>& matrix);
template float sum<float>(const Matrix2d<float>& matrix);
template double sum<double>(const Matrix2d<double>& matrix);
template std::complex<float> sum<std::complex<float>>(
    const Matrix2d<std::complex<float>>& matrix);
template std::complex<double> sum<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix);

template int dot<int>(const Matrix2d<int>& left, const Matrix2d<int>& right);
template float dot<float>(const Matrix2d<float>& left,
                          const Matrix2d<float>& right);
template double dot<double>(const Matrix2d<double>& left,
                            const Matrix2d<double>& right);
template std::complex<float> dot<std::complex<float>>(
    const Matrix2d<std::complex<float>>& left,
    const Matrix2d<std::complex<float>>& right);
template std::complex<double> dot<std::complex<double>>(
    const Matrix2d<std::complex<double>>& left,
    const Matrix2d<std::complex<double>>& right);

template double norm<int>(const Matrix2d<int>& matrix);
template double norm<float>(const Matrix2d<float>& matrix);
template double norm<double>(const Matrix2d<double>& matrix);
template double norm<std::complex<float>>(
    const Matrix2d<std::complex<float>>& matrix);
template double norm<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix);

template double max_abs<int>(const Matrix2d<int>& matrix);
template double max_abs<float>(const Matrix2d<float>& matrix);
template double max_abs<double>(const Matrix2d<double>& matrix);
template double max_abs<std::complex<float>>(
    const Matrix2d<std::complex<float>>& matrix);
template double max_abs<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix);

template std::vector<int> diagonal<int>(const Matrix2d<int>& matrix);
template std::vector<float> diagonal<float>(const Matrix2d<float>& matrix);
template std::vector<double> diagonal<double>(const Matrix2d<double>& matrix);
template std::vector<std::complex<float>> diagonal<std::complex<float>>(
    const Matrix2d<std::complex<float>>& matrix);
template std::vector<std::complex<double>> diagonal<std::complex<double>>(
    const Matrix2d<std::complex<double>>& matrix);
} // namespace matrix2d_impl
