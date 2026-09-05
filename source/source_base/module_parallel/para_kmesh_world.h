#ifndef PARA_KMESH_WORLD_H
#define PARA_KMESH_WORLD_H

#include <complex>
#include <vector>

#include "para_world.h"

namespace Parallel
{

/**
 * @brief k-mesh parallel domain: pure k-point pool (k-pool) topology.
 *
 * Self-contained replacement for Parallel_Kpoints + KP_WORLD +
 * GlobalV::KPAR / MY_POOL. Owns the k-point pool layout and provides
 * query / collection operations.
 *
 * The k-pool split (kpar) is independent of bndpar: the domain knows
 * nothing about band groups. Reductions that must span both dimensions
 * (e.g. the total electron count under BPCG) therefore combine
 * ParaBdiffKsameWorld::reduce_across_bdiff_ksame (band dimension, run FIRST)
 * with reduce_across_pools (k dimension, run SECOND).
 *
 * In serial builds all operations degenerate to single-pool behavior.
 * Tests only need this header; no GlobalV, no parallel_comm.h.
 */
class ParaKmeshWorld : public ParaWorld
{
public:
    /**
     * @brief Construct a serial (single-pool) k-mesh domain.
     *
     * @param[in] nkstot  total number of k-points (without spin)
     * @param[in] nspin   number of spin components
     */
    ParaKmeshWorld(int nkstot, int nspin);

    /**
     * @brief Construct a reduce-only k-mesh domain with no k-point
     *        distribution data.
     *
     * kpar_/comm are set from the bridge globals so that
     * reduce_across_pools / reduce_max/min_across_pools work correctly.
     * The distribution arrays (nks_pool_, whichpool_, ...) are left
     * empty; calling pool_collection / gather_kvec is invalid.
     */
    ParaKmeshWorld();

#ifdef __MPI
    /**
     * @brief Construct a k-mesh domain on an existing communicator.
     *
     * @param[in] comm     communicator spanning every process of every
     *                     k-pool (MPI_COMM_WORLD in the current bridge)
     * @param[in] kpar     number of k-pools
     * @param[in] my_pool  k-pool index of this process
     * @param[in] nkstot   total number of k-points (without spin)
     * @param[in] nspin    number of spin components
     */
    ParaKmeshWorld(const MPI_Comm& comm, int kpar, int my_pool, int nkstot, int nspin);
#endif

    /// Number of k-pools.
    int kpar() const { return kpar_; }

    /// k-pool index of this process.
    int my_pool() const { return my_pool_; }

    /// Rank within the communicator (world rank in the bridge layout).
    int rank_in_pool() const { return rank_in_pool_; }

    /// Total number of processes.
    int nproc() const { return nproc_; }

    /// Number of spin components.
    int nspin() const { return nspin_; }

    /// Total number of k-points (without spin).
    int nkstot() const { return nkstot_; }

    /// Number of k-points in this k-pool.
    int nks_local() const { return nks_local_; }

    /// Global start index of this k-pool's k-points.
    int startk_global() const { return startk_global_; }

    /// Number of k-points in the given k-pool.
    int nks_pool(int pool) const;

    /// Global start index of the given k-pool's k-points.
    int startk_pool(int pool) const;

    /// Which k-pool owns the given global k-point index.
    int which_pool(int ik_global) const;

    /// First communicator rank of the given k-pool.
    int startpro_pool(int pool) const;

    /// Maximum number of k-points across all k-pools.
    int max_nks_pool() const;

    /// Whether this process is the first process of its k-pool.
    bool kpool_root() const { return kpool_root_; }

    // ===== Cross-pool reductions =====

    /**
     * @brief Sum a scalar across the k-pools: one contribution per pool.
     *
     * Replaces Parallel_Reduce::reduce_double_allpool. The first process
     * of each k-pool injects the partial sum, all other processes inject
     * zero, so a single world-wide MPI_Allreduce yields the sum of the
     * per-k-pool partial sums. Correct for uneven k-pool sizes and free
     * of the legacy normalization division (which divided by an average
     * pool size and double-counted pools of uneven layouts).
     *
     * Precondition: with band parallelism (bndpar > 1) the caller must
     * first combine the band-group partial sums (e.g.
     * ParaBdiffKsameWorld::reduce_across_bdiff_ksame) so that every process of a
     * k-pool holds one complete per-pool partial sum.
     *
     * No-op when kpar() <= 1.
     *
     * @param[in,out] value  local partial sum, overwritten with global total
     */
    void reduce_across_pools(double& value) const;

    /**
     * @brief Global max across all k-pools and band groups.
     *
     * Max/min are idempotent, so a plain world-wide Allreduce is correct
     * for every pool layout and covers both the k-pool and the band-group
     * dimension (band-parallel shards see different eigenvalue windows,
     * so the Fermi-level bounds must be extremized across both).
     * Replaces Parallel_Reduce::reduce_max (all of MPI_COMM_WORLD).
     * No-op when this domain spans a single process.
     *
     * @param[in,out] value  local value, overwritten with global max
     */
    void reduce_max_across_pools(double& value) const;

    /**
     * @brief Global min across all k-pools and band groups.
     *
     * @param[in,out] value  local value, overwritten with global min
     */
    void reduce_min_across_pools(double& value) const;

    // ===== Cross-domain operations =====

    /**
     * @brief Collect a scalar value from the pool that owns k-point ik.
     *
     * Pool 0 receives the value; other pools send. Only rank_in_pool==0
     * participates in the actual communication.
     *
     * @param[out] value  collected value (valid on pool 0 root)
     * @param[in]  wk     local k-point weights array
     * @param[in]  ik     global k-point index
     */
    void pool_collection(double& value, const double* wk, int ik) const;

    /**
     * @brief Collect an array slice from the pool that owns k-point ik.
     *
     * @param[out] value  output array (dim elements)
     * @param[in]  w      input array (nkstot * dim elements, row-major by k)
     * @param[in]  dim    number of elements per k-point
     * @param[in]  ik     global k-point index
     */
    template <typename T>
    void pool_collection(T* value, const T* w, int dim, int ik) const;

    /**
     * @brief Gather local k-point vectors to global array.
     *
     * Only pool-root processes contribute their local k-points;
     * the result is valid on all ranks after MPI_Allreduce.
     *
     * @param[in]  vec_local   local k-point vectors (nks_local elements)
     * @param[out] vec_global  global k-point vectors (nkstot elements)
     */
    void gather_kvec(const std::vector<double>& vec_local, std::vector<double>& vec_global) const;

private:
    void distribute_kpoints();

    int kpar_ = 1;
    int my_pool_ = 0;
    int rank_in_pool_ = 0;
    int nproc_ = 1;
    int nspin_ = 1;
    int nkstot_ = 0;
    int nks_local_ = 0;
    int startk_global_ = 0;
    bool kpool_root_ = true; ///< first process of my k-pool (reduction contributor)

    std::vector<int> nks_pool_;      ///< k-points per k-pool
    std::vector<int> startk_pool_;   ///< global start index per k-pool
    std::vector<int> whichpool_;     ///< k-pool index per k-point
    std::vector<int> startpro_pool_; ///< first communicator rank per k-pool
};

} // namespace Parallel

#endif // PARA_KMESH_WORLD_H
