#ifndef PARALLEL_KPOINTS_H
#define PARALLEL_KPOINTS_H

#include "module_base/complexarray.h"
#include "module_base/global_function.h"
#include "module_base/realarray.h"
#include "module_base/vector3.h"

class Parallel_Kpoints
{
  public:
    Parallel_Kpoints(){};
    ~Parallel_Kpoints(){};

    void kinfo(int& nkstot_in,
               const int& kpar_in,
               const int& my_pool_in,
               const int& rank_in_pool_in,
               const int& nproc_in,
               const int& nspin_in);

    // collect value from each pool to wk.
    void pool_collection(double& value, const double* wk, const int& ik);

    // collect value from each pool to overlap.
    void pool_collection(double* valuea,
                         double* valueb,
                         const ModuleBase::realArray& a,
                         const ModuleBase::realArray& b,
                         const int& ik);
    void pool_collection(std::complex<double>* value, const ModuleBase::ComplexArray& w, const int& ik) const;
    template <class T, class V>
    void pool_collection_aux(T* value, const V& w, const int& dim, const int& ik) const;
#ifdef __MPI
    /**
     * @brief gather kpoints from all processors
     *
     * @param vec_local kpoint vector in local processor
     * @param vec_global kpoint vector in all processors
     */
    void gatherkvec(const std::vector<ModuleBase::Vector3<double>>& vec_local,
                    std::vector<ModuleBase::Vector3<double>>& vec_global) const;
#endif

    // information about pool, dim: KPAR
    // int* nproc_pool = nullptr;    it is not used

    // inforamation about kpoints, dim: KPAR
    std::vector<int> nks_pool;    // number of k-points in each pool, here use k-points without spin
    std::vector<int> startk_pool; // the first k-point in each pool, here use k-points without spin

    // information about which pool each k-point belongs to,
    std::vector<int> whichpool; // whichpool[k] : the pool which k belongs to, dim: nkstot_np

    int nkstot_np = 0; // number of k-points without spin, kv.set_nkstot(nkstot_np) * nspin(1 or 2)
    int nks_np = 0;    // number of k-points without spin in the present pool

    // get the first processor in the pool
    int get_startpro_pool(const int& pool) const

    {
        return startpro_pool[pool];
    }

    // get the maximum number of k-points in all pools
    int get_max_nks_pool() const
    {
        return *std::max_element(nks_pool.begin(), nks_pool.end());
    }

  private:

    int kpar = 0;         // number of pools
    int my_pool = 0;      // the pool index of the present processor
    int rank_in_pool = 0; // the rank in the present pool
    int nproc = 1;        // number of processors
    int nspin = 1;        // number of spins

    std::vector<int> startpro_pool; // the first processor in each pool
#ifdef __MPI
    void get_nks_pool(const int& nkstot);
    void get_startk_pool(const int& nkstot);
    void get_whichpool(const int& nkstot);

    void set_startpro_pool();
#endif
};

#endif
