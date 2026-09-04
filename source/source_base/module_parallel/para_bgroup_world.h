#ifndef PARA_BGROUP_WORLD_H
#define PARA_BGROUP_WORLD_H

#include "para_world.h"

namespace Parallel
{

/**
 * @brief bgroup parallel domain: band group communication topology.
 *
 * Self-contained replacement for INT_BGROUP + BP_WORLD +
 * GlobalV::MY_BNDGROUP/NPROC_IN_BNDGROUP/RANK_IN_BPGROUP.
 *
 * The band group domain has two communicators:
 *   - intra: INT_BGROUP (same band group, different k/pw)
 *   - inter: BP_WORLD (different band groups, same k)
 *
 * Tests only need this header.
 */
class ParaBgroupWorld : public ParaWorld
{
public:
    /**
     * @brief Construct a serial bgroup domain (single band group).
     */
    ParaBgroupWorld();

#ifdef __MPI
    /**
     * @brief Construct a bgroup domain from intra and inter communicators.
     *
     * @param[in] intra_comm  intra-group communicator (e.g. INT_BGROUP)
     * @param[in] inter_comm  inter-group communicator (e.g. BP_WORLD)
     * @param[in] nbndgroup   number of band groups
     */
    ParaBgroupWorld(const MPI_Comm& intra_comm, const MPI_Comm& inter_comm, int nbndgroup);
#endif

    /// Band group index of this process.
    int my_bndgroup() const { return my_bndgroup_; }

    /// Number of band groups.
    int nbndgroup() const { return nbndgroup_; }

    /// Rank within the band group (alias for rank()).
    int rank_in_bpgroup() const { return rank(); }

    /// Number of processes in the band group (alias for size()).
    int nproc_in_bndgroup() const { return size(); }

#ifdef __MPI
    /// Inter-group communicator (BP_WORLD equivalent).
    MPI_Comm inter_comm() const { return inter_comm_; }
#endif

private:
    int my_bndgroup_ = 0;
    int nbndgroup_ = 1;
#ifdef __MPI
    MPI_Comm inter_comm_ = MPI_COMM_NULL;
#endif
};

} // namespace Parallel

#endif // PARA_BGROUP_WORLD_H
