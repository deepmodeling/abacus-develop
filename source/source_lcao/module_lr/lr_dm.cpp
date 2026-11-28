#include "source_lcao/module_lr/lr_dm.h"

#include "source_cell/module_neighbor/sltk_grid_driver.h"
#include "source_cell/unitcell.h"
#include "source_estate/module_dm/dm_from_psi.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"

#include <complex>

namespace LR
{

template <typename T>
module_dm::DensityMatrix<T, double> make_gs_dm(const psi::Psi<T>& psi_all, const ModuleBase::matrix& wg_all,
                                               const Parallel_Orbitals& pv_all, const Parallel_Orbitals& pmat,
                                               const UnitCell& ucell, const Grid_Driver& gd,
                                               const std::vector<double>& orb_cutoff, int nspin,
                                               const std::vector<ModuleBase::Vector3<double>>& kvec_d, int nk)
{
    module_dm::DensityMatrix<T, double> dm(&pmat, nspin, kvec_d, nk);
    // `psi_all` is distributed per `pv_all`: using the wrong descriptor reads its local block with
    // the wrong block size/process-grid assumption, which is invisible at nprocs=1.
    module_dm::dm_from_psi(&pv_all, wg_all, psi_all, dm);     // nbands is important here
    LR_Util::initialize_DMR(dm, pmat, ucell, gd, orb_cutoff); // nbands is not important here
    dm.cal_dmr(-1);
    return dm;
}

template module_dm::DensityMatrix<double, double> make_gs_dm<double>(
    const psi::Psi<double>&, const ModuleBase::matrix&, const Parallel_Orbitals&, const Parallel_Orbitals&,
    const UnitCell&, const Grid_Driver&, const std::vector<double>&, int,
    const std::vector<ModuleBase::Vector3<double>>&, int);
template module_dm::DensityMatrix<std::complex<double>, double> make_gs_dm<std::complex<double>>(
    const psi::Psi<std::complex<double>>&, const ModuleBase::matrix&, const Parallel_Orbitals&,
    const Parallel_Orbitals&, const UnitCell&, const Grid_Driver&, const std::vector<double>&, int,
    const std::vector<ModuleBase::Vector3<double>>&, int);

} // namespace LR
