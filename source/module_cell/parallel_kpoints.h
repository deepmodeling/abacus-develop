#ifndef PARALLEL_KPOINTS_H
#define PARALLEL_KPOINTS_H

#include "module_base/complexarray.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/realarray.h"
#include "module_base/vector3.h"

class Parallel_Kpoints
{
  public:
    Parallel_Kpoints();
    ~Parallel_Kpoints();

    void kinfo(int& nkstot);

    // collect value from each pool to wk.
    void pool_collection(double& value, const double* wk, const int& ik);

    // collect value from each pool to overlap.
    void pool_collection(double* valuea,
                         double* valueb,
                         const ModuleBase::realArray& a,
                         const ModuleBase::realArray& b,
                         const int& ik);
    void pool_collection(std::complex<double>* value, const ModuleBase::ComplexArray& w, const int& ik);
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

    // information about pool, dim: GlobalV::KPAR
    int* nproc_pool = nullptr;
    int* startpro_pool = nullptr;

    // inforamation about kpoints, dim: GlobalV::KPAR
    int* nks_pool = nullptr;    // number of k-points in each pool
    int* startk_pool = nullptr; // the first k-point in each pool
    int kpar = 0;               // number of pools

    // information about which pool each k-point belongs to,
    int* whichpool = nullptr; // whichpool[k] : the pool which k belongs to, dim: nkstot_np
    int nkstot_np = 0;        // number of k-points without spin, kv.nkstot = nkstot_np * nspin(1 or 2)
    int nks_np = 0;           // number of k-points without spin in the present pool

  private:
#ifdef __MPI
    void get_nks_pool(const int& nkstot);
    void get_startk_pool(const int& nkstot);
    void get_whichpool(const int& nkstot);
#endif
};

#endif
