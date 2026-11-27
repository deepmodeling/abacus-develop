#ifndef PARA_KMESH_WORLD_H
#define PARA_KMESH_WORLD_H

#include <complex>
#include <vector>

#include "para_world.h"

namespace Parallel
{

/**
 * @brief k-mesh parallel domain: k-point distribution across pools.
 *
 * Self-contained replacement for Parallel_Kpoints + KP_WORLD +
 * GlobalV::KPAR / MY_POOL / RANK_IN_POOL. Owns all k-point pool
 * topology data and provides query / collection operations.
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

#ifdef __MPI
    /**
     * @brief Construct a k-mesh domain on an existing communicator.
     *
     * @param[in] comm     k-point pool communicator (e.g. KP_WORLD)
     * @param[in] kpar     number of pools
     * @param[in] my_pool  pool index of this process
     * @param[in] nproc    total number of processes (MPI_COMM_WORLD size)
     * @param[in] nkstot   total number of k-points (without spin)
     * @param[in] nspin    number of spin components
     */
    ParaKmeshWorld(const MPI_Comm& comm, int kpar, int my_pool, int nproc, int nkstot, int nspin);
#endif

    /// Number of pools.
    int kpar() const { return kpar_; }

    /// Pool index of this process.
    int my_pool() const { return my_pool_; }

    /// Rank within the pool.
    int rank_in_pool() const { return rank_in_pool_; }

    /// Total number of processes.
    int nproc() const { return nproc_; }

    /// Number of spin components.
    int nspin() const { return nspin_; }

    /// Total number of k-points (without spin).
    int nkstot() const { return nkstot_; }

    /// Number of k-points in this pool.
    int nks_local() const { return nks_local_; }

    /// Global start index of this pool's k-points.
    int startk_global() const { return startk_global_; }

    /// Number of k-points in the given pool.
    int nks_pool(int pool) const;

    /// Global start index of the given pool's k-points.
    int startk_pool(int pool) const;

    /// Which pool owns the given global k-point index.
    int which_pool(int ik_global) const;

    /// First MPI_COMM_WORLD rank of the given pool.
    int startpro_pool(int pool) const;

    /// Maximum number of k-points across all pools.
    int max_nks_pool() const;

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

    std::vector<int> nks_pool_;      ///< k-points per pool
    std::vector<int> startk_pool_;   ///< global start index per pool
    std::vector<int> whichpool_;     ///< pool index per k-point
    std::vector<int> startpro_pool_; ///< first world rank per pool
};

} // namespace Parallel

#endif // PARA_KMESH_WORLD_H
