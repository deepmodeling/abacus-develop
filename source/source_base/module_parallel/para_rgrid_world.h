#ifndef PARA_RGRID_WORLD_H
#define PARA_RGRID_WORLD_H

#include <vector>

#include "para_world.h"

namespace Parallel
{

/**
 * @brief rgrid parallel domain: real-space FFT grid distribution.
 *
 * Self-contained replacement for GRID_WORLD + GlobalV::GRANK/GSIZE +
 * Parallel_Grid's z-distribution tables. Owns grid dimensions and the
 * per-process z-plane allocation (numz/startz/whichpro).
 *
 * Cross-pool operations (reduce_across_pools, bcast, reduce) will accept
 * communicators as parameters rather than reading global POOL_WORLD/
 * KP_WORLD, breaking the cross-domain dependency.
 *
 * Tests only need this header.
 */
class ParaRgridWorld : public ParaWorld
{
public:
    /**
     * @brief Construct a serial rgrid domain (all z-planes on one process).
     *
     * @param[in] ncx, ncy, ncz  global grid dimensions
     */
    ParaRgridWorld(int ncx, int ncy, int ncz);

#ifdef __MPI
    /**
     * @brief Construct an rgrid domain on an existing communicator.
     *
     * @param[in] comm   grid communicator (e.g. GRID_WORLD)
     * @param[in] ncx, ncy, ncz  global grid dimensions
     */
    ParaRgridWorld(const MPI_Comm& comm, int ncx, int ncy, int ncz);
#endif

    /// Global grid dimension in x.
    int ncx() const { return ncx_; }

    /// Global grid dimension in y.
    int ncy() const { return ncy_; }

    /// Global grid dimension in z.
    int ncz() const { return ncz_; }

    /// Local z-plane count for this process.
    int nczp() const { return nczp_; }

    /// Total real-space grid points on this process (ncx * ncy * nczp).
    int nrxx() const { return ncx_ * ncy_ * nczp_; }

    /// Number of z-planes assigned to process p in this pool.
    int numz(int p) const;

    /// Starting global z-index for process p.
    int startz(int p) const;

    /// Which process owns global z-plane iz.
    int whichpro(int iz) const;

    // ===== Cross-domain operations =====

    /**
     * @brief Sum local grid data across all pools (KP_WORLD or INT_BGROUP).
     *
     * Replaces Parallel_Grid::reduce_across_pools.
     * In serial mode or single-pool, this is a no-op.
     *
     * @param[in,out] data  local grid buffer (nrxx elements), overwritten with sum
     * @param[in] kmesh_world  k-mesh domain providing the cross-pool communicator
     */
    void reduce_across_pools(double* data, const ParaWorld& kmesh_world) const;

    /**
     * @brief Broadcast global grid to local z-slabs (replaces Parallel_Grid::bcast).
     *
     * @param[in] data_global  global grid (ncxyz elements, only valid on root)
     * @param[out] data_local  local grid buffer (nrxx elements)
     * @param[in] comm_world   communicator for broadcast
     * @param[in] root         root rank in comm_world
     */
    void bcast_data(const double* data_global, double* data_local,
                    const ParaWorld& comm_world, int root = 0) const;

    /**
     * @brief Gather local z-slabs into a global grid (replaces Parallel_Grid::reduce).
     *
     * @param[out] rhotot    global grid (ncxyz elements, only valid on root)
     * @param[in] rhoin      local grid buffer (nrxx elements)
     * @param[in] comm_world   communicator for gather
     */
    void reduce_data(double* rhotot, const double* rhoin,
                     const ParaWorld& comm_world) const;

private:
    void distribute_z();

    int ncx_ = 0;
    int ncy_ = 0;
    int ncz_ = 0;
    int nczp_ = 0;  ///< local z-plane count

    std::vector<int> numz_;      ///< z-planes per process
    std::vector<int> startz_;     ///< start z-index per process
    std::vector<int> whichpro_;   ///< owner of each global z-plane
};

} // namespace Parallel

#endif // PARA_RGRID_WORLD_H
