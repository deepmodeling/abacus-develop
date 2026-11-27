#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_MPI_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_MPI_H

/**
 * @file matrix2d_mpi.h
 * @internal
 * @brief MPI bookkeeping shared by Matrix2d backend implementations.
 *
 * This header is not part of the public Matrix2d interface. Callers should
 * include `source_base/matrix2d.h` or `source_base/matrix2d_linalg.h` instead.
 */

#ifdef __MPI
#include "mpi.h"

#include <cstddef>
#include <vector>

namespace matrix2d_impl
{
namespace mpi
{
int checked_mpi_count(const std::size_t nbytes);

/**
 * @brief Check a rank-local byte count after collective agreement.
 *
 * Every rank must call this function with the same communicator. A local
 * overflow is reduced across ranks so either all ranks throw or all ranks
 * continue to the subsequent communication call.
 */
int checked_mpi_count_collective(const std::size_t nbytes,
                                 const MPI_Comm comm);

std::vector<int> make_byte_offsets(const std::vector<int>& nbytes);

std::vector<int> allgather_nbytes(const int nbytes_local,
                                  const MPI_Comm comm);
} // namespace mpi
} // namespace matrix2d_impl
#endif

#endif
