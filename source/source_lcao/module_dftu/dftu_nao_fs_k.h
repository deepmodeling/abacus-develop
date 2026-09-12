#ifndef DFTU_NAO_FS_K_H
#define DFTU_NAO_FS_K_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_base/matrix.h"
#include "source_cell/klist.h"
#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"

#include <cassert>
#include <complex>
#include <string>
#include <vector>


class ForceStressArrays;
class Plus_U_Base;

namespace DFTU_LCAO {

/// @brief Immutable shared environment for the DFT+U force/stress kernels.
///
/// Bundles the system objects (unit cell, neighbor grid, orbital
/// parallelization, folded-matrix arrays and DFT+U state) together with
/// the scalar configuration (orbital cutoffs, ks_solver name and spinor
/// count) that stays constant for one force/stress calculation. It holds
/// no per-k-point data: density matrices, k vectors and output matrices
/// remain explicit function arguments.
class DftuFsEnv
{
  public:
    DftuFsEnv(Plus_U_Base& dftu,
              const UnitCell& ucell,
              const Grid_Driver& gd,
              const Parallel_Orbitals& pv,
              ForceStressArrays& fsr,
              const std::vector<double>& orb_cutoff,
              const std::string& ks_solver,
              const int npol)
        : dftu_(&dftu),
          ucell_(&ucell),
          gd_(&gd),
          pv_(&pv),
          fsr_(&fsr),
          orb_cutoff_(&orb_cutoff),
          ks_solver_(ks_solver),
          npol_(npol)
    {
        assert(npol_ == 1 || npol_ == 2);
    }

    /// @brief DFT+U state; the onsite-potential builders take it non-const.
    Plus_U_Base& dftu() const
    {
        return *dftu_;
    }

    const UnitCell& ucell() const
    {
        return *ucell_;
    }

    const Grid_Driver& gd() const
    {
        return *gd_;
    }

    const Parallel_Orbitals& pv() const
    {
        return *pv_;
    }

    /// @brief Folded dS/dH matrix arrays; buffers are filled by the caller.
    ForceStressArrays& fsr() const
    {
        return *fsr_;
    }

    const std::vector<double>& orb_cutoff() const
    {
        return *orb_cutoff_;
    }

    const std::string& ks_solver() const
    {
        return ks_solver_;
    }

    int npol() const
    {
        return npol_;
    }

  private:
    Plus_U_Base* dftu_;
    const UnitCell* ucell_;
    const Grid_Driver* gd_;
    const Parallel_Orbitals* pv_;
    ForceStressArrays* fsr_;
    const std::vector<double>* orb_cutoff_;
    std::string ks_solver_;
    int npol_;
};

/// @brief Top-level entry: drives force/stress from DFT+U.
void force_stress(const DftuFsEnv& env,
                  const bool cal_force,
                  const bool cal_stress,
                  std::vector<std::vector<double>>* dmk_d,
                  std::vector<std::vector<std::complex<double>>>* dmk_c,
                  ModuleBase::matrix& force_dftu,
                  ModuleBase::matrix& stress_dftu,
                  const K_Vectors& kv,
                  const bool gamma_only_local);

} // namespace DFTU_LCAO

#endif // DFTU_NAO_FS_K_H
