#include "source_base/matrix2d_mpi.h"

#ifdef __MPI
#include "source_base/module_parallel/para_mpi_func.h"
#include "source_base/module_parallel/para_world.h"

#include <limits>
#include <stdexcept>

namespace matrix2d_impl
{
namespace mpi
{
int checked_mpi_count(const std::size_t nbytes)
{
    // FIXME: BLACS and ScaLAPACK interfaces still expose native int
    // dimensions. Revisit this size_t-to-int boundary when those interfaces
    // are migrated together.
    if (nbytes > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument(
            "Matrix2d MPI buffer exceeds the native int count range");
    }
    return static_cast<int>(nbytes);
}

std::vector<int> make_byte_offsets(const std::vector<int>& nbytes)
{
    std::vector<int> offsets(nbytes.size(), 0);
    int offset = 0;
    for (std::size_t ip = 0; ip < nbytes.size(); ++ip)
    {
        offsets[ip] = offset;
        if (nbytes[ip] > std::numeric_limits<int>::max() - offset)
        {
            throw std::invalid_argument(
                "Matrix2d MPI buffer size overflow");
        }
        offset += nbytes[ip];
    }
    return offsets;
}

std::vector<int> allgather_nbytes(const int nbytes_local,
                                  const MPI_Comm comm)
{
    const Parallel::ParaWorld world
        = Parallel::ParaWorld::make_mpi("matrix2d_mpi_counts", comm);
    if (!world.valid())
    {
        throw std::invalid_argument(
            "Matrix2d MPI count collection requires a valid communicator");
    }

    int rank_nbytes = nbytes_local;
    std::vector<int> nbytes(static_cast<std::size_t>(world.size()), 0);
    Parallel::allgather_int(rank_nbytes, nbytes.data(), world);
    return nbytes;
}
} // namespace mpi
} // namespace matrix2d_impl
#endif
