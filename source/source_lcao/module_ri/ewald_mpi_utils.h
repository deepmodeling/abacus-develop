#ifndef EWALD_MPI_UTILS_H
#define EWALD_MPI_UTILS_H

#include <RI/distribute/Distribute_Equally.h>
#include <array>
#include <cstddef>
#include <mpi.h>
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
} // namespace EwaldVqDetail

#endif
