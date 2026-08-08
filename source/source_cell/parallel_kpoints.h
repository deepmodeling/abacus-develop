/**
 * @file parallel_kpoints.h
 * @brief Parallel_Kpoints class for k-point parallelization.
 */
#ifndef PARALLEL_KPOINTS_H
#define PARALLEL_KPOINTS_H

#include "source_base/complexarray.h"
#include "source_base/global_function.h"
#include "source_base/realarray.h"
#include "source_base/vector3.h"

/**
 * @brief Parallel_Kpoints class for k-point parallelization.
 */
class Parallel_Kpoints
{
  public:
    /**
     * @brief Default constructor.
     */
    Parallel_Kpoints(){};

    /**
     * @brief Destructor.
     */
    ~Parallel_Kpoints(){};

    /**
     * @brief Initialize k-point parallelization information.
     *
     * @param nkstot_in total number of k-points
     * @param kpar_in number of pools
     * @param my_pool_in pool index of current processor
     * @param rank_in_pool_in rank within pool
     * @param nproc_in number of processors
     * @param nspin_in number of spin components
     */
    void kinfo(int& nkstot_in,
               const int& kpar_in,
               const int& my_pool_in,
               const int& rank_in_pool_in,
               const int& nproc_in,
               const int& nspin_in);

    /**
     * @brief Collect value from each pool to wk.
     *
     * @param value output collected value
     * @param wk k-point weights
     * @param ik k-point index
     */
    void pool_collection(double& value, const double* wk, const int& ik);

    /**
     * @brief Collect value from each pool to overlap.
     *
     * @param valuea output collected value a
     * @param valueb output collected value b
     * @param a input array a
     * @param b input array b
     * @param ik k-point index
     */
    void pool_collection(double* valuea,
                         double* valueb,
                         const ModuleBase::realArray& a,
                         const ModuleBase::realArray& b,
                         const int& ik);

    /**
     * @brief Collect complex value from each pool.
     *
     * @param value output collected value
     * @param w input complex array
     * @param ik k-point index
     */
    void pool_collection(std::complex<double>* value, const ModuleBase::ComplexArray& w, const int& ik) const;

    /**
     * @brief Auxiliary template function for pool collection.
     *
     * @param value output collected value
     * @param w input array
     * @param dim dimension
     * @param ik k-point index
     */
    template <class T, class V>
    void pool_collection_aux(T* value, const V& w, const int& dim, const int& ik) const;

#ifdef __MPI
    /**
     * @brief Gather k-points from all processors.
     *
     * @param vec_local k-point vector in local processor
     * @param vec_global k-point vector in all processors
     */
    void gatherkvec(const std::vector<ModuleBase::Vector3<double>>& vec_local,
                    std::vector<ModuleBase::Vector3<double>>& vec_global) const;
#endif

    /// K-point information
    std::vector<int> nks_pool;    ///< number of k-points in each pool (without spin)
    std::vector<int> startk_pool; ///< the first k-point in each pool (without spin)
    std::vector<int> whichpool;   ///< whichpool[k]: the pool which k belongs to, dim: nkstot_np

    int nkstot_np = 0; ///< number of k-points without spin
    int nks_np = 0;    ///< number of k-points without spin in the present pool

    /**
     * @brief Get the first processor in the pool.
     *
     * @param pool pool index
     * @return first processor in the pool
     */
    int get_startpro_pool(const int& pool) const
    {
        return startpro_pool[pool];
    }

    /**
     * @brief Get the maximum number of k-points in all pools.
     * @return maximum number of k-points
     */
    int get_max_nks_pool() const
    {
        return *std::max_element(nks_pool.begin(), nks_pool.end());
    }

  public:
    int kpar = 0;         ///< number of pools
    int my_pool = 0;      ///< the pool index of the present processor
    int rank_in_pool = 0; ///< the rank in the present pool
    int nproc = 1;        ///< number of processors
    int nspin = 1;        ///< number of spins

  private:
    std::vector<int> startpro_pool; ///< the first processor in each pool

#ifdef __MPI
    /**
     * @brief Get number of k-points in each pool.
     *
     * @param nkstot total number of k-points
     */
    void get_nks_pool(const int& nkstot);

    /**
     * @brief Get start k-point index for each pool.
     *
     * @param nkstot total number of k-points
     */
    void get_startk_pool(const int& nkstot);

    /**
     * @brief Get which pool each k-point belongs to.
     *
     * @param nkstot total number of k-points
     */
    void get_whichpool(const int& nkstot);

    /**
     * @brief Set first processor for each pool.
     */
    void set_startpro_pool();
#endif
};

#endif
