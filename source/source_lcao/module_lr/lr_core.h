#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_LR_CORE_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_LR_CORE_H

#include "source_base/matrix.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_psi/psi.h"
#include "source_lcao/module_lr/lr_dim.h"
#include "source_lcao/module_lr/lr_gs.h"
#include "source_lcao/module_lr/lr_mpl.h"
#include "source_lcao/module_lr/lr_relax.h"
#include "source_lcao/module_lr/potentials/lr_gs_pot.h"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"

#include <memory>
#include <string>
#include <vector>

namespace LR
{
/// The ground state the LR solver works from.
template <typename T>
struct GsState
{
    std::unique_ptr<psi::Psi<T>> psi;
    psi::Psi<T>* psi_all = nullptr;
    std::unique_ptr<psi::Psi<T>> psi_all_own;
    ModuleBase::matrix eig;
    ModuleBase::matrix eig_all;
    ModuleBase::matrix wg;
    ModuleBase::matrix wg_all;
    Parallel_Orbitals para_mat;
    Parallel_Orbitals para_mat_all;
};

/// The excitation amplitudes and their distribution.
struct Amp
{
    std::vector<ct::Tensor> x;
    int nloc = 1;
    Parallel_Blkcyclic para_c;
    std::vector<Parallel_Blkcyclic> para_x;
    std::vector<std::string> spin_types;
};

/// The LR potentials, including the ground-state one needed by the force path.
struct PotSet
{
    std::vector<std::shared_ptr<PotHxcLR>> pot;
    std::shared_ptr<PotHxcLR> pot_hxc_gs;
    GsPot gs_pot;
};

/// State of an excited-state geometry relaxation.
template <typename T>
struct RelaxRun
{
    RelaxTarget<T> target;
    double etot_gs = 0.0;
    ModuleBase::matrix force_gs;
    ModuleBase::matrix lr_force;
    std::vector<MplLvc> mpl;
};

/// LR solver state that does not belong to the esolver interface.
template <typename T>
struct Core
{
    Dim dim;
    Amp amp;
    GsState<T> gs;
    ZWin<T> zw;
    PotSet pot;
    RelaxRun<T> relax;

    /// Read the ground-state bands this solver needs, from file.
    void read_wfc(const std::string& in_dir, const Input_para& inp, const K_Vectors& kv, std::ostream& log);
    /// Widen the ground-state window into the Z window used by the force path.
    void fill_z(const int* desc_src, const K_Vectors& kv, bool cal_force, bool first, std::ostream& log);
};

} // namespace LR

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_CORE_H
