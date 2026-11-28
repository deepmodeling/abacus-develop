#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_DM_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_DM_H

#include "source_base/matrix.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_psi/psi.h"

#include <vector>

class Grid_Driver;
class UnitCell;

namespace LR
{
/// Ground-state density matrix (DMR) from all-band KS orbitals. `psi_all` is distributed per
/// `pv_all`: the KS solver's `pv` on the ks-lr path, or the all-band descriptor when read from file.
template <typename T>
module_dm::DensityMatrix<T, double> make_gs_dm(const psi::Psi<T>& psi_all, const ModuleBase::matrix& wg_all,
                                               const Parallel_Orbitals& pv_all, const Parallel_Orbitals& pmat,
                                               const UnitCell& ucell, const Grid_Driver& gd,
                                               const std::vector<double>& orb_cutoff, int nspin,
                                               const std::vector<ModuleBase::Vector3<double>>& kvec_d, int nk);
} // namespace LR

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_DM_H
