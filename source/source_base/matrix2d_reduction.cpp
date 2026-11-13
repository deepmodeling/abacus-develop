#include "source_base/matrix2d_reduction.h"

#include "source_base/parallel_common.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
template <typename T>
T conjugate_value(const T& value)
{
    return value;
}

template <typename T>
std::complex<T> conjugate_value(const std::complex<T>& value)
{
    return std::conj(value);
}

template <typename T>
double absolute_square(const T& value)
{
    const double real_value = static_cast<double>(value);
    return real_value * real_value;
}

template <typename T>
double absolute_square(const std::complex<T>& value)
{
    const double real_value = static_cast<double>(value.real());
    const double imaginary_value = static_cast<double>(value.imag());
    return real_value * real_value + imaginary_value * imaginary_value;
}

template <typename Value, typename Element>
void reduce_sum_over_processes(Value& value, const Matrix2d<Element>& matrix)
{
#ifdef __MPI
    const std::shared_ptr<const Parallel_Blkcyclic> dist
        = matrix.distribution();
    if (dist != nullptr && !dist->is_serial())
    {
        const MPI_Comm comm = dist->comm();
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
    const std::shared_ptr<const Parallel_Blkcyclic> dist
        = matrix.distribution();
    if (dist != nullptr && !dist->is_serial())
    {
        const MPI_Comm comm = dist->comm();
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

int checked_mpi_count(const std::size_t count)
{
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Matrix2d diagonal payload exceeds MPI count range");
    }
    return static_cast<int>(count);
}

std::vector<int> make_byte_displacements(const std::vector<int>& counts)
{
    std::vector<int> displacements(counts.size(), 0);
    int displacement = 0;
    for (std::size_t ip = 0; ip < counts.size(); ++ip)
    {
        displacements[ip] = displacement;
        if (counts[ip] > std::numeric_limits<int>::max() - displacement)
        {
            throw std::invalid_argument("Matrix2d diagonal payload size overflow");
        }
        displacement += counts[ip];
    }
    return displacements;
}

#ifdef __MPI
std::vector<int> gather_byte_counts(const int local_count, const MPI_Comm comm)
{
    const int nproc = Parallel_Common::comm_size(comm);
    const std::vector<int> receive_counts(
        static_cast<std::size_t>(nproc),
        checked_mpi_count(sizeof(int)));
    const std::vector<int> receive_displacements
        = make_byte_displacements(receive_counts);
    std::vector<int> counts(static_cast<std::size_t>(nproc), 0);
    Parallel_Common::allgatherv_bytes(&local_count,
                                      checked_mpi_count(sizeof(local_count)),
                                      counts.data(),
                                      receive_counts.data(),
                                      receive_displacements.data(),
                                      comm);
    return counts;
}
#endif
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
                partial_sums[ip] += matrix.at_local(i, j);
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
            value += matrix.at_local(i, j);
        }
    }
#endif
    reduce_sum_over_processes(value, matrix);
    return value;
}

template <typename T>
double sum_squares(const Matrix2d<T>& matrix)
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
            value += absolute_square(matrix.at_local(i, j));
        }
    }
    reduce_sum_over_processes(value, matrix);
    return value;
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
                partial_sums[ip]
                    += conjugate_value(left.at_local(i, j))
                        * right.at_local(i, j);
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
            value += conjugate_value(left.at_local(i, j))
                * right.at_local(i, j);
        }
    }
#endif
    reduce_sum_over_processes(value, left);
    return value;
}

template <typename T>
double norm2(const Matrix2d<T>& matrix)
{
    return std::sqrt(sum_squares(matrix));
}

template <typename T>
double max_abs(const Matrix2d<T>& matrix)
{
    double value = 0.0;
    const int nrloc = matrix.nrloc();
    const int ncloc = matrix.ncloc();
#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(max : value) schedule(static)
#endif
    for (int i = 0; i < nrloc; ++i)
    {
        for (int j = 0; j < ncloc; ++j)
        {
            value = std::max(
                value,
                static_cast<double>(std::abs(matrix.at_local(i, j))));
        }
    }
    reduce_max_over_processes(value, matrix);
    return value;
}

template <typename T>
std::vector<T> diagonal(const Matrix2d<T>& matrix)
{
    const int diagonal_size = std::min(matrix.nrglb(), matrix.ncglb());
    std::vector<T> values(static_cast<std::size_t>(diagonal_size), T());
    const int nrloc = matrix.nrloc();
    const int ncloc = matrix.ncloc();
    const std::shared_ptr<const Parallel_Blkcyclic> dist
        = matrix.distribution();
    std::vector<std::pair<int, T>> local_entries;
    for (int irloc = 0; irloc < nrloc; ++irloc)
    {
        for (int icloc = 0; icloc < ncloc; ++icloc)
        {
            const int irglb = dist == nullptr ? irloc : dist->irglb(irloc);
            const int icglb = dist == nullptr ? icloc : dist->icglb(icloc);
            if (irglb == icglb && irglb < diagonal_size)
            {
                const T value = matrix.at_local(irloc, icloc);
                local_entries.push_back(std::pair<int, T>(irglb, value));
            }
        }
    }
#ifdef __MPI
    if (dist != nullptr && !dist->is_serial())
    {
        const MPI_Comm comm = dist->comm();
        if (comm != MPI_COMM_NULL)
        {
            const std::size_t entry_bytes = sizeof(int) + sizeof(T);
            const std::size_t local_payload_bytes
                = sizeof(int) + local_entries.size() * entry_bytes;
            std::vector<unsigned char> local_payload(local_payload_bytes);
            const int local_entry_count
                = checked_mpi_count(local_entries.size());
            std::memcpy(local_payload.data(),
                        &local_entry_count,
                        sizeof(local_entry_count));
            std::size_t offset = sizeof(local_entry_count);
            for (const std::pair<int, T>& entry : local_entries)
            {
                std::memcpy(local_payload.data() + offset,
                            &entry.first,
                            sizeof(entry.first));
                offset += sizeof(entry.first);
                std::memcpy(local_payload.data() + offset,
                            &entry.second,
                            sizeof(entry.second));
                offset += sizeof(entry.second);
            }

            const int local_count = checked_mpi_count(local_payload.size());
            const std::vector<int> payload_counts
                = gather_byte_counts(local_count, comm);
            const std::vector<int> payload_displacements
                = make_byte_displacements(payload_counts);
            const int total_bytes
                = payload_displacements.back() + payload_counts.back();
            std::vector<unsigned char> received(
                static_cast<std::size_t>(total_bytes));
            Parallel_Common::allgatherv_bytes(local_payload.data(),
                                              local_count,
                                              received.data(),
                                              payload_counts.data(),
                                              payload_displacements.data(),
                                              comm);

            for (std::size_t ip = 0; ip < payload_counts.size(); ++ip)
            {
                const std::size_t begin
                    = static_cast<std::size_t>(payload_displacements[ip]);
                const std::size_t end
                    = begin + static_cast<std::size_t>(payload_counts[ip]);
                if (end - begin < sizeof(int))
                {
                    throw std::runtime_error(
                        "Matrix2d diagonal payload is truncated");
                }
                int entry_count = 0;
                std::memcpy(&entry_count,
                            received.data() + begin,
                            sizeof(entry_count));
                const std::size_t expected_bytes
                    = static_cast<std::size_t>(entry_count) * entry_bytes;
                if (entry_count < 0
                    || end - begin != sizeof(int) + expected_bytes)
                {
                    throw std::runtime_error(
                        "Matrix2d diagonal payload has an invalid size");
                }
                std::size_t entry_offset = begin + sizeof(entry_count);
                for (int ientry = 0; ientry < entry_count; ++ientry)
                {
                    int index = 0;
                    T value = T();
                    std::memcpy(&index,
                                received.data() + entry_offset,
                                sizeof(index));
                    entry_offset += sizeof(index);
                    std::memcpy(&value,
                                received.data() + entry_offset,
                                sizeof(value));
                    entry_offset += sizeof(value);
                    if (index < 0 || index >= diagonal_size)
                    {
                        throw std::runtime_error(
                            "Matrix2d diagonal index is outside global shape");
                    }
                    values[static_cast<std::size_t>(index)] += value;
                }
            }
        }
    }
#endif
    if (dist == nullptr || dist->is_serial())
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

template double sum_squares<int>(const Matrix2d<int>& matrix);
template double sum_squares<float>(const Matrix2d<float>& matrix);
template double sum_squares<double>(const Matrix2d<double>& matrix);
template double sum_squares<std::complex<float>>(
    const Matrix2d<std::complex<float>>& matrix);
template double sum_squares<std::complex<double>>(
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

template double norm2<int>(const Matrix2d<int>& matrix);
template double norm2<float>(const Matrix2d<float>& matrix);
template double norm2<double>(const Matrix2d<double>& matrix);
template double norm2<std::complex<float>>(
    const Matrix2d<std::complex<float>>& matrix);
template double norm2<std::complex<double>>(
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
