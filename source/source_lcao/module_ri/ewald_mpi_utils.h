#ifndef EWALD_MPI_UTILS_H
#define EWALD_MPI_UTILS_H

#include "ewald_tail_utils.h"

#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"

#include <RI/distribute/Distribute_Equally.h>
#include <array>
#include <cstddef>
#include <mpi.h>
#include <stdexcept>
#include <vector>

namespace EwaldVqDetail
{
template <typename TA, typename Tcell, std::size_t Ndim>
auto distribute_common_realspace_tasks(const MPI_Comm& mpi_comm,
                                       const std::vector<TA>& atoms,
                                       const std::array<Tcell, Ndim>& bare_period)
{
    return RI::Distribute_Equally::distribute_atoms_periods(mpi_comm, atoms, bare_period, 2, false);
}

inline std::vector<TailKey> distribute_tail_shell_keys(const std::vector<int>& atoms,
                                                       const TailCell& production_period,
                                                       const TailCell& probe_period,
                                                       const int mpi_rank,
                                                       const int mpi_size)
{
    if (mpi_rank < 0 || mpi_rank >= mpi_size)
    {
        throw std::invalid_argument("MPI rank is outside the communicator");
    }
    const auto shell_cells = make_shell_cells(production_period, probe_period);
    std::vector<TailKey> local_shell;
    const std::size_t global_key_count = atoms.size() * atoms.size() * shell_cells.size();
    local_shell.reserve(global_key_count / static_cast<std::size_t>(mpi_size) + 1);
    for (const int atom_i: atoms)
    {
        for (const int atom_j: atoms)
        {
            for (const TailCell& cell: shell_cells)
            {
                const TailKey key{atom_i, atom_j, cell};
                if (shell_owner(key, mpi_size) == mpi_rank) local_shell.push_back(key);
            }
        }
    }
    return local_shell;
}

inline TailStats reduce_tail_stats(const MPI_Comm mpi_comm, const TailStats& local)
{
    const int mpi_rank = Parallel_Common::my_rank(mpi_comm);

    TailStats global = local;
    unsigned long long counts[2] = {static_cast<unsigned long long>(local.production_blocks),
                                    static_cast<unsigned long long>(local.shell_blocks)};
    Parallel_Reduce::reduce_all(mpi_comm, counts, 2);
    global.production_blocks = static_cast<std::size_t>(counts[0]);
    global.shell_blocks = static_cast<std::size_t>(counts[1]);

    Parallel_Reduce::reduce_max(mpi_comm, global.reference_norm);
    Parallel_Reduce::reduce_all(mpi_comm, global.sum_norm);
    Parallel_Reduce::reduce_max(mpi_comm, global.hermitian_residual);

    // The worst offender is located together with the rank that owns it, so
    // every rank can then broadcast that rank's diagnostic values.
    double max_norm = local.max_norm;
    int max_rank = mpi_rank;
    Parallel_Reduce::reduce_max_loc(mpi_comm, max_norm, max_rank);
    global.max_norm = max_norm;

    int coverage_complete = local.coverage_complete ? 1 : 0;
    Parallel_Reduce::reduce_min(mpi_comm, coverage_complete);
    global.coverage_complete = coverage_complete != 0;

    int worst_key[5]
        = {local.worst_key.atom_i,
           local.worst_key.atom_j,
           local.worst_key.cell[0],
           local.worst_key.cell[1],
           local.worst_key.cell[2]};
    double worst_values[3]
        = {local.worst_bare_norm, local.worst_gaussian_norm, local.worst_distance};
    Parallel_Common::bcast_int(worst_key, 5, mpi_comm, max_rank);
    Parallel_Common::bcast_double(worst_values, 3, mpi_comm, max_rank);
    global.worst_key = TailKey{worst_key[0], worst_key[1], {worst_key[2], worst_key[3], worst_key[4]}};
    global.worst_bare_norm = worst_values[0];
    global.worst_gaussian_norm = worst_values[1];
    global.worst_distance = worst_values[2];
    return global;
}
} // namespace EwaldVqDetail

#endif
