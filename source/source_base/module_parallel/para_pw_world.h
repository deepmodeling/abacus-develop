#ifndef PARA_PW_WORLD_H
#define PARA_PW_WORLD_H

#include <vector>

#include "para_world.h"

namespace Parallel
{

/**
 * @brief pw parallel domain: plane-wave distribution within a pool.
 *
 * Self-contained replacement for poolnproc/poolrank/npw_per/npwtot
 * members scattered across PW_Basis and GlobalV::NPROC_IN_POOL.
 * Owns the pool-level parallel topology and plane-wave count distribution.
 *
 * The actual FFT-based distribution algorithm (method1/method2) stays
 * in PW_Basis; this class only holds the result: how many plane waves
 * each process in the pool gets.
 *
 * Tests only need this header; no PW_Basis, no parallel_comm.h.
 */
class ParaPwWorld : public ParaWorld
{
public:
    /**
     * @brief Construct a serial (single-process) pw domain.
     *
     * @param[in] npw  number of plane waves on this process
     */
    explicit ParaPwWorld(int npw);

#ifdef __MPI
    /**
     * @brief Construct a pw domain on an existing pool communicator.
     *
     * @param[in] comm   pool communicator (e.g. POOL_WORLD)
     * @param[in] npw_per  array of plane-wave counts per process (size = pool size)
     */
    ParaPwWorld(const MPI_Comm& comm, const std::vector<int>& npw_per);
#endif

    /// Number of plane waves on this process.
    int npw() const { return npw_; }

    /// Total number of plane waves in the pool.
    int npwtot() const { return npwtot_; }

    /// Number of plane waves on process p in the pool.
    int npw_per(int p) const;

    /// Number of processes in the pool (same as size()).
    int poolnproc() const { return size(); }

    /// Rank within the pool (same as rank()).
    int poolrank() const { return rank(); }

private:
    int npw_ = 0;               ///< local plane-wave count
    int npwtot_ = 0;            ///< total plane waves in pool
    std::vector<int> npw_per_;  ///< per-process plane-wave counts
};

} // namespace Parallel

#endif // PARA_PW_WORLD_H
