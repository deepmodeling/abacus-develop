#ifndef ABACUS_SOURCE_ESOLVER_ESOLVER_LR_H
#define ABACUS_SOURCE_ESOLVER_ESOLVER_LR_H

#include "source_cell/unitcell.h"
#include "source_esolver/esolver_fp.h"
#include "source_estate/elecstate.h"
#include "source_estate/elecstate_lcao.h"
#include "source_hamilt/hamilt.h"
#include "source_io/module_parameter/input_parameter.h"

#include <memory>
#include <string>
#include <vector> //future tensor

#include "source_esolver/esolver_ks_lcao.h" //for the move constructor
#include "source_estate/module_dm/density_matrix.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_hamilt/module_gint/gint_info.h"
#include "source_lcao/module_lr/hamilt_casida.h"
#include "source_lcao/module_lr/potentials/lr_gs_pot.h"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
#include "source_lcao/module_lr/lr_mpl.h"
#include "source_lcao/module_lr/lr_core.h"
#include "source_lcao/module_lr/lr_sforce.h"
#include "source_lcao/module_lr/lr_zeq.h"
#include "source_lcao/module_lr/lr_relax.h"
#ifdef __EXX
// #include <RI/physics/Exx.h>
#include "source_hamilt/module_xc/exx_info.h" // for Exx_Info value member
#include "source_lcao/module_ri/exx_lri.h"
#endif
namespace LR
{
template <typename T>
struct GradientInputs;
}
namespace ModuleESolver
{
struct LRRuntime
{
    std::ofstream &running_log;
    std::ofstream &warning_log;
    int rank;
};

/// Excited State Solver: Linear Response TDDFT (Tamm Dancoff Approximation)
template <typename T, typename TR = double>
class ESolver_LR : public ModuleESolver::ESolver_FP
{
  public:
    ESolver_LR(
        const Input_para &inp,
        const std::string &in_dir,
        const std::string &out_dir,
        const LRRuntime &runtime);
    ~ESolver_LR() {}

    /// input: input, call, basis(LCAO), psi(ground state), elecstate
    // initialize sth. independent of the ground state
    virtual void before_all_runners(BaseCell &basecell, const Input_para &inp) override;
    virtual void runner(BaseCell &basecell, int istep) override;
    virtual void after_all_runners(BaseCell &basecell) override;

    /// Total energy of the excited state being relaxed: E_gs + Omega. Zero outside a
    /// relaxation, where nothing consumes it and the old behaviour is kept.
    virtual double cal_energy() override;
    /// Force on the atoms in the relaxed excited state, F = -d(E_gs + Omega)/dR (Ry/Bohr).
    virtual void cal_force(BaseCell &basecell, ModuleBase::matrix &force) override;
    /// Not implemented: there is no excited-state stress. `cell-relax` is rejected at input.
    virtual void cal_stress(BaseCell &basecell, ModuleBase::matrix &stress) override;

  protected:
    const std::string in_dir;
    const std::string out_dir;
    const UnitCell *ucell_ = nullptr;
    std::vector<double> orb_cutoff_;

    /// cached aliases, read once here instead of at every log/print call site below
    std::ofstream &ofs_running_;
    std::ofstream &ofs_warning_;
    const int my_rank_;

    /// @brief the ground-state solver, kept alive across ionic steps (esolver_type = "ks-lr").
    /// Null on the `lr` path, where the ground state comes from files instead.
    std::unique_ptr<ModuleESolver::ESolver_KS_LCAO<T, TR>> ks_;

    // Geometry-dependent objects that the ground-state solver already builds for the current
    // structure. They are aliased rather than rebuilt: `ESolver_KS_LCAO::before_scf` refreshes
    // its own copies every ionic step, so recomputing them here would be both wasted work and
    // a chance for the two to disagree. On the `lr` path the pointers are bound to this
    // object's own members below, which `initialize_from_unitcell_` fills from file.
    Grid_Driver gd_own_;                    ///< only used when `ks_` is null
    TwoCenterBundle two_center_bundle_own_; ///< only used when `ks_` is null
    Grid_Driver *gd_ptr_ = nullptr;
    const TwoCenterBundle *tcb_ptr_ = nullptr;
    Parallel_Grid *pgrid_ptr_ = nullptr;
    Structure_Factor *sf_ptr_ = nullptr;
    pseudopot_cell_vl *locpp_ptr_ = nullptr;

    Grid_Driver &gd() const { return *this->gd_ptr_; }
    const TwoCenterBundle &tcb() const { return *this->tcb_ptr_; }
    Parallel_Grid &pgrid() const { return *this->pgrid_ptr_; }
    Structure_Factor &sfac() const { return *this->sf_ptr_; }
    pseudopot_cell_vl &vloc() const { return *this->locpp_ptr_; }
    /// bind the aliases above; `ks_` must already be set (or null for the `lr` path)
    void bind_ground_state_aliases_();

    // @brief only needed for force calculation

    std::string xc_kernel;

    void initialize_from_unitcell_(UnitCell &ucell, const Input_para &inp);
    /// one-time setup from the ground-state solver; ends by calling `refresh_from_ks_`
    void initialize_from_ks_(UnitCell &ucell, const Input_para &inp);
    /// re-read everything that depends on the atomic positions, once per ionic step
    void refresh_from_ks_(UnitCell &ucell);
    bool ks_initialized_ = false; ///< whether `initialize_from_ks_` has already run
    bool exx_owned_ = false;      ///< `exx_lri` was built here (so its Cs/Vs are ours to refresh)

    std::unique_ptr<ModuleGint::GintInfo> gint_info_ = nullptr;

    /// LR solver state that is not part of the esolver interface.
    LR::Core<T> lr_;


    // ---------------- the Z-vector (CPSCF) window ----------------
    // The Z-vector equation enforces the Brillouin condition in EVERY occupied-virtual
    // rotation, so it must not be confined to the `nvirt` window X lives in.
    // It should be the whole AO virtual space.
    //
    // These mirror `psi_ks` / `eig_ks` / `paraC_` / `paraX_` / `nvirt` / `nloc_per_state`
    // but span every virtual band the ground state produced (the input `nbands`), so the
    // window is widened by raising *nbands*, not `nvirt`. X keeps its own window, so Omega
    // -- and with it any finite-difference reference -- is untouched.
    /// Widen `nst` X blocks starting at `istate_begin` from the X window into the Z window,
    /// zero-filling the virtual rows X does not have.
    /// @brief variables for parallel distribution of matrix in AO representation
    LCAO_Orbitals orb_;                            ///< numerical atomic orbital data for single-point evaluation
    std::vector<std::complex<double>> velocity_mo; ///< store the velocity matrix elements in MO basis

    /// @brief read in the ground state wave function, band energy and occupation
    virtual void read_ks_wfc();
    /// @brief  read in the ground state charge density
    void read_ks_chg(Charge &chg);

    /// Build the LR kernel potentials (overridden by BSE).
    virtual void init_pot(const Charge &chg_gs);

    /// @brief check the legality of the input parameters
    /// @brief set nocc, nvirt, nbasis, npairs and nstates
    void set_dim();
    /// reset nocc, nvirt, npairs after read ground-state wavefunction when nspin=2
    void reset_dim_spin();

    ///========================== for gradient calculation =========================
    /// Borrow the geometry-dependent objects `LR::Core` needs.
    LR::Refs<T> make_refs();

#ifdef __EXX
    /// Tdata of Exx_LRI is same as T, for the reason, see operator_lr_exx.h
    std::shared_ptr<Exx_LRI<T>> exx_lri = nullptr;
    Exx_Info exx_info;
#endif
};
} // namespace ModuleESolver

#endif // ABACUS_SOURCE_ESOLVER_ESOLVER_LR_H
