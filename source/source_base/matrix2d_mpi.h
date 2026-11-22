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

std::vector<int> make_byte_offsets(const std::vector<int>& nbytes);

std::vector<int> allgather_nbytes(const int nbytes_local,
                                  const MPI_Comm comm);
} // namespace mpi
} // namespace matrix2d_impl
#endif

#endif
