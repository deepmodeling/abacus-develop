#ifndef STERNHEIMER_SIAB_MPI_H
#define STERNHEIMER_SIAB_MPI_H

#include "sternheimer_siab_data.h"

#include <complex>
#include <vector>

#ifdef __MPI
#include <mpi.h>
#endif

namespace module_ri
{
namespace sternheimer_siab
{

struct PrimitiveSlab
{
    int startz = 0;
    int nplane = 0;
    /// Primitive-major values in PW local order [primitive][ixy][local iz].
    std::vector<std::vector<std::complex<double>>> values;
};

std::vector<std::vector<std::complex<double>>> assemble_full_primitive_grids(const std::vector<PrimitiveSlab>& slabs,
                                                                             int nxy,
                                                                             int nz);

#ifdef __MPI
std::vector<std::vector<std::complex<double>>> allgather_full_primitive_grids(
    const std::vector<std::vector<std::complex<double>>>& local_primitives,
    int nxy,
    int nz,
    int startz,
    int nplane,
    MPI_Comm communicator);

std::vector<ReferenceRow> gather_reference_rows_to_root(const std::vector<ReferenceRow>& local_rows,
                                                        std::size_t nprimitive,
                                                        int root,
                                                        MPI_Comm communicator);
#else
std::vector<std::vector<std::complex<double>>> allgather_full_primitive_grids(
    const std::vector<std::vector<std::complex<double>>>& local_primitives,
    int nxy,
    int nz,
    int startz,
    int nplane);

std::vector<ReferenceRow> gather_reference_rows_to_root(const std::vector<ReferenceRow>& local_rows,
                                                        std::size_t nprimitive,
                                                        int root);
#endif

} // namespace sternheimer_siab
} // namespace module_ri

#endif
