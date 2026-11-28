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

#include <fstream>
#include <memory>
#include <string>
#include <vector>

class Grid_Driver;
class Structure_Factor;
class TwoCenterBundle;
class pseudopot_cell_vl;
class Charge;
class surchem;
#ifdef __EXX
template <typename Tdata>
class Exx_LRI;
#endif
namespace module_dm
{
template <typename TK, typename TR>
class DensityMatrix;
}

namespace LR
{
template <typename T>
struct GradientInputs;
template <typename T>
struct GsForceCtx;
template <typename T>
struct StateCtx;

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

/// All-band ground-state data owned by the ground-state solver, handed over for one refresh.
template <typename T>
struct KsBands
{
    psi::Psi<T>* psi = nullptr;              ///< all-band KS orbitals
    const ModuleBase::matrix* ekb = nullptr; ///< all-band eigenvalues
    const ModuleBase::matrix* wg = nullptr;  ///< all-band occupations
    const Parallel_Orbitals* pv = nullptr;   ///< distribution of `psi`
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

/// Borrowed handles the force and relaxation path needs. The esolver binds them once, after the
/// ground state has been read; every pointer outlives the solver.
template <typename T>
struct Refs
{
    const UnitCell* ucell = nullptr;
    const K_Vectors* kv = nullptr;
    Grid_Driver* gd = nullptr;
    const std::vector<double>* orb_cutoff = nullptr;
    ModulePW::PW_Basis* rho = nullptr;
    ModulePW::PW_Basis* rhod = nullptr;
    pseudopot_cell_vl* vloc = nullptr;
    Structure_Factor* sfac = nullptr;
    surchem* solvent = nullptr;
    const TwoCenterBundle* tcb = nullptr;
    const Input_para* inp = nullptr;
    std::ofstream* log = nullptr;
    Parallel_Grid* pgrid = nullptr;
    const std::string* in_dir = nullptr;
    const std::string* out_dir = nullptr;
    const std::string* xc_kernel = nullptr;
    /// The file (`lr`) path owns `sf`/`vloc`; on `ks-lr` the ground-state solver refreshed them.
    bool refresh_sfac = true;
    bool has_float_data = false;
    int my_rank = 0;
    /// The layout `gs.psi_all` is distributed by: the KS solver's `pv`, or `gs.para_mat_all`.
    const Parallel_Orbitals* pv_all = nullptr;
#ifdef __EXX
    /// Points at the solver's member so a rebuilt `Exx_LRI` is seen without rebinding.
    const std::shared_ptr<Exx_LRI<T>>* exx_lri = nullptr;
    double hybrid_alpha = 1.0;
#endif
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
    Refs<T> refs;

    void bind(const Refs<T>& r) { this->refs = r; }

    /// Read the ground-state bands this solver needs, from file.
    void read_wfc(const std::string& in_dir, const Input_para& inp, const K_Vectors& kv, std::ostream& log);
    /// Widen the ground-state window into the Z window used by the force path.
    void fill_z(const int* desc_src, const K_Vectors& kv, bool cal_force, bool first, std::ostream& log);
    /// Re-read the band window that follows the current geometry, then rebuild the Z window.
    /// `first` selects the one-time spin-channel setup; later calls only validate it.
    void refresh_gs(const KsBands<T>& bands, const K_Vectors& kv, const Input_para& inp, bool first,
                    std::ostream& log);
    /// Set up the AO/band distributions and the band window on the ground-state solver's layout.
    void setup_ks_layout(const Parallel_Orbitals& pv, const K_Vectors& kv, const UnitCell& ucell,
                         const Input_para& inp);
    /// Set up the same distributions for the standalone path, where the ground state comes from files.
    void setup_file_layout(const K_Vectors& kv, const UnitCell& ucell, const Input_para& inp);
    /// Solve the TDA eigenproblem, or read its states from file under `lr_solver = spectrum`.
    /// The `nspin x nstates` eigenvalues are written to `ekb`.
    void solve(double* ekb, int istep);
    /// Electron-hole densities, the absorption spectrum and the per-state gradients after the run.
    void output(const ModuleBase::matrix& ekb);
    /// Allocate the excitation amplitudes and fill their initial guess.
    void setup_x();
    /// Build the LR kernel potentials (and the ground-state one used by the force path) on `chg_gs`.
    void init_pot(const Charge& chg_gs);

    /// Excited-state gradients d(Omega)/dR, one matrix per state solved. `istate_only >= 0`
    /// restricts it to that one state; -1 does all `nstates`.
    std::vector<ModuleBase::matrix> cal_force(int ispin, const ModuleBase::matrix& ekb, int istate_only);
    /// Gradients for caller-supplied vectors already widened into the Z window.
    std::vector<ModuleBase::matrix> cal_force_xz(int ispin, const ct::Tensor& xz, const std::vector<double>& omega,
                                                 int label_begin);
    /// Open-shell counterpart of `cal_force_xz`: both channels live in one vector.
    std::vector<ModuleBase::matrix> cal_force_os(const ct::Tensor& xz, const std::vector<double>& omega,
                                                 int label_begin);
    /// Per-state gradients plus the gradient matrix of each degenerate multiplet.
    void cal_force_grad(int ispin, const ModuleBase::matrix& ekb, std::ofstream& ofs);
    /// Ground-state density matrix.
    module_dm::DensityMatrix<T, double> cal_dm_gs() const;
    /// Reproduce the ground-state force as a check of the force machinery.
    void test_force();

    /// Resolve and validate `lr_target_state` / `lr_target_spin`.
    void setup_relax();
    /// Re-choose the followed root by cross-geometry overlap, then refresh the references.
    void follow_target(const ModuleBase::matrix& ekb, std::ofstream& ofs);
    /// The LR half of the force for the current geometry, following the selected surface.
    ModuleBase::matrix cal_lr_force(const ModuleBase::matrix& ekb, std::ofstream& ofs);
    /// The force of the steepest-descending branch of the target multiplet (`lr_degen_mode = jt`).
    ModuleBase::matrix cal_jt_force(const std::vector<ModuleBase::matrix>& diag, const ModuleBase::matrix& ekb,
                                    std::ofstream& ofs);
    /// Widen a multiplet's eigenvectors into the Z window, one block each.
    ct::Tensor pad_group_z(int ispin, const std::vector<int>& group) const;

  private:
    GradientInputs<T> gradient_inputs() const;
    GsForceCtx<T> gs_force_ctx() const;
    StateCtx<T> state_ctx(const ModuleBase::matrix& ekb) const;
};

} // namespace LR

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_LR_CORE_H
