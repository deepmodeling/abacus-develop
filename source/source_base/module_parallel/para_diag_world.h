#ifndef PARA_DIAG_WORLD_H
#define PARA_DIAG_WORLD_H

#include "para_world.h"

namespace Parallel
{

/**
 * @brief diag parallel domain: diagonalization group topology.
 *
 * Self-contained replacement for DIAG_WORLD + GlobalV::DRANK/DSIZE/DCOLOR.
 * The diag domain is created by splitting MPI_COMM_WORLD into groups
 * for parallel diagonalization (PEXSI, ScaLAPACK).
 *
 * Tests only need this header; no parallel_comm.h or parallel_global.h.
 */
class ParaDiagWorld : public ParaWorld
{
public:
    /**
     * @brief Construct a serial diag domain (single-process group).
     */
    ParaDiagWorld();

#ifdef __MPI
    /**
     * @brief Construct a diag domain from an existing communicator.
     *
     * @param[in] comm   diag communicator (e.g. DIAG_WORLD)
     * @param[in] dcolor color used in MPI_Comm_split to create this group
     */
    ParaDiagWorld(const MPI_Comm& comm, int dcolor);
#endif

    /// Color used in MPI_Comm_split to create this diag group.
    int dcolor() const { return dcolor_; }

    /// Rank within the diag group (alias for rank()).
    int drank() const { return rank(); }

    /// Number of processes in the diag group (alias for size()).
    int dsize() const { return size(); }

private:
    int dcolor_ = 0;
};

} // namespace Parallel

#endif // PARA_DIAG_WORLD_H
