#include "sternheimer_chi0_mpi.h"

#ifdef __MPI

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace module_ri
{
namespace sternheimer_chi0
{

void reduce_branch_to_root(std::vector<std::complex<double>>& branch, const int root, MPI_Comm communicator)
{
    if (communicator == MPI_COMM_NULL)
    {
        throw std::invalid_argument("Sternheimer chi0 reduction requires a valid MPI communicator.");
    }

    int rank = -1;
    int rank_count = 0;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &rank_count);
    if (root < 0 || root >= rank_count)
    {
        throw std::invalid_argument("Sternheimer chi0 reduction root is outside the MPI communicator.");
    }

    const unsigned long long local_size = static_cast<unsigned long long>(branch.size());
    unsigned long long minimum_size = local_size;
    unsigned long long maximum_size = local_size;
    MPI_Allreduce(MPI_IN_PLACE, &minimum_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, communicator);
    MPI_Allreduce(MPI_IN_PLACE, &maximum_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, communicator);
    if (minimum_size != maximum_size)
    {
        throw std::invalid_argument("Sternheimer chi0 branch sizes differ between MPI ranks.");
    }

    std::vector<std::complex<double>> reduced;
    if (rank == root)
    {
        reduced.assign(branch.size(), std::complex<double>(0.0, 0.0));
    }
    std::size_t offset = 0;
    while (offset != branch.size())
    {
        const std::size_t remaining = branch.size() - offset;
        const int count = static_cast<int>(
            std::min(remaining, static_cast<std::size_t>(std::numeric_limits<int>::max())));
        std::complex<double>* receive_buffer = rank == root ? reduced.data() + offset : nullptr;
        if (MPI_Reduce(branch.data() + offset,
                       receive_buffer,
                       count,
                       MPI_DOUBLE_COMPLEX,
                       MPI_SUM,
                       root,
                       communicator)
            != MPI_SUCCESS)
        {
            throw std::runtime_error("Sternheimer chi0 MPI reduction failed.");
        }
        offset += static_cast<std::size_t>(count);
    }

    if (rank == root)
    {
        branch.swap(reduced);
    }
    else
    {
        branch.clear();
        branch.shrink_to_fit();
    }
}

} // namespace sternheimer_chi0
} // namespace module_ri

#endif
