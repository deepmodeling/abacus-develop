/**
 * @file spin_constrain.h
 * @brief Core header for the DeltaSpin (spin-constrained DFT) module.
 *
 * @par Purpose
 * Implements constrained local spin density (CLSD) functional calculations,
 * where atomic magnetic moments are constrained to target values via
 * Lagrange multipliers (lambda). The constrained energy functional is:
 *   E'[rho] = E[rho] - sum_i lambda_i . (M_i - M_target_i)
 * where lambda_i is the Lagrange multiplier (magnetic force) on atom i,
 * M_i is the computed magnetic moment, and M_target_i is the target moment.
 *
 * @par Algorithm
 * The lambda optimization uses a conjugate-gradient-like scheme (run_lambda_loop):
 *   1. Compute magnetic moments Mi from current wavefunction
 *   2. Calculate residual: delta_spin = Mi - M_target
 *   3. Build search direction (steepest descent or Polak-Ribiere CG)
 *   4. Apply lambda update: lambda += alpha * search_direction
 *   5. Re-diagonalize Hamiltonian with DeltaSpin correction
 *   6. Compute new Mi, find optimal alpha via linear interpolation
 *   7. Repeat until RMS(delta_spin) < sc_thr
 *
 * @par Basis Set Support
 * - LCAO: Uses real-space projection via DeltaSpin operator on density matrix
 * - PW (Plane Wave): Uses subspace diagonalization with OnsiteProjector becp coefficients
 *
 * @par Spin Types
 * - nspin=2 (collinear): Only z-component constrained, npol=1, uses spin_sign (+1/-1)
 * - nspin=4 (non-collinear): Full xyz components constrained, npol=2, full Pauli matrices
 *
 * @par Convergence Criteria
 * - RMS error: sqrt(mean(delta_spin^2)) < sc_thr (adaptive threshold)
 * - Gradient decay: max(dM/dlambda) per atom type < decay_grad[itype]
 * - Maximum steps: nsc (default 50), minimum steps: nsc_min
 *
 * @par Internal layout
 * All basis-set-independent constraint data (lambda, target_mag, Mi, constrain,
 * indexing maps, loop parameters) is owned by the ScState member `state_`
 * (see deltaspin_state.h). The public setters/getters below are thin
 * forwarding shells kept for backward compatibility with existing call sites.
 */
#ifndef SPIN_CONSTRAIN_H
#define SPIN_CONSTRAIN_H

#include <complex>
#include <fstream>
#include <map>
#include <vector>

#include "source_base/constants.h"
#include "source_base/complexmatrix.h"
#include "source_base/matrix.h"
#include "source_base/tool_quit.h"
#include "source_base/tool_title.h"
#include "source_base/vector3.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_cell/klist.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/operator.h"
#include "source_estate/elecstate.h"

#include "deltaspin_state.h"
#include "deltaspin_pw_cache.h"

#ifdef __LCAO
#include "source_estate/module_dm/density_matrix.h" // mohan add 2025-11-02
#endif

namespace spinconstrain
{

/**
 * @brief Singleton class implementing spin-constrained DFT (DeltaSpin).
 *
 * @par Template parameter TK
 * - std::complex<double>: Used for nspin=4 (non-collinear) and internally for nspin=2
 * - double: Stub specialization for nspin=2 collinear (all methods are no-ops)
 *
 * @par Design rationale
 * - Singleton pattern: Only one SpinConstrain instance per TK type is needed,
 *   shared across the SCF loop. Prevents duplicate state management.
 * - void* pointers (p_hamilt, psi): Type-erased to avoid template dependency cycles
 *   with the Hamiltonian and Psi classes. Cast to concrete types at call sites.
 * - subspace data caching (sub_h_save, sub_s_save, becp_save): For PW basis, the
 *   subspace Hamiltonian and becp are computed once per SCF iteration and reused
 *   across multiple lambda steps, avoiding expensive re-computation.
 *
 * @par Key workflow (PW basis):
 *   SCF iteration -> run_lambda_loop()
 *     -> cal_mw_from_lambda() [first call saves subspace data]
 *       -> calculate_delta_hcc() [H += becp^† * lambda * becp]
 *       -> diag_responce() [subspace diagonalization, update becp]
 *       -> accumulate_Mi_from_becp() [compute magnetic moments]
 *     -> BFGS optimizer updates lambda
 *     -> Repeat until RMS(Mi - M_target) < sc_thr
 *     -> update_psi_charge() [final full-space update]
 *
 * @par Error handling
 * - assert(sub_h_save != nullptr): Called before subspace operations;
 *   failure means cal_mw_from_lambda() was not called before update_psi_charge().
 *   Solution: Ensure cal_mw_from_lambda() is called at least once per SCF step.
 * - "atomCounts is not set": init_sc() was not called or UnitCell data is missing.
 * - "nspin must be 2 or 4": Invalid spin configuration. nspin=1 is not supported.
 */
template <typename TK>
class SpinConstrain
{
public:
    /**
     * =============================================================
     * PUBLIC INTERFACE - Main entry points for the ESolver layer
     * =============================================================
     */

    /**
     * @brief Master initialization: populate all SC parameters from UnitCell and input.
     *
     * @details Called once at the start of a DeltaSpin calculation. Performs:
     *   1. Set input parameters (convergence threshold, max steps, trial alpha)
     *   2. Get atom/orbital/lnchi counts from UnitCell for indexing
     *   3. Set nspin and npol (nspin=4 -> npol=2, nspin=2 -> npol=1)
     *   4. Load target_mag, lambda, constrain from UnitCell (parsed from STRU)
     *   5. For nspin=2: force x,y constraint flags to 0 (collinear: only z is constrained)
     *   6. Set solver parameters (k-point list, Hamiltonian, psi, electronic state)
     *
     * @param sc_thr_in Convergence threshold for RMS(Mi - M_target) in uB
     * @param nsc_in Maximum number of inner lambda optimization steps
     * @param nsc_min_in Minimum number of inner steps before early exit checks
     * @param alpha_trial_in Initial trial step size (eV/uB^2), converted to Ry internally
     * @param sccut_in Maximum lambda change per step (eV/uB), converted to Ry internally
     * @param sc_drop_thr_in Fraction of initial RMS for adaptive threshold
     * @param ucell Unit cell with atomic positions, STRU constraint data
     * @param direction_only_in If true, only optimize spin direction (|lambda| -> 0)
     * @param ParaV_in Parallel orbitals distribution info (LCAO only)
     * @param nspin_in Spin type: 2=collinear, 4=non-collinear
     * @param kv_in K-point vector list
     * @param p_hamilt_in Pointer to Hamiltonian (HamiltLCAO or HamiltPW)
     * @param psi_in Pointer to wavefunctions (Psi<TK>)
     * @param dm_in Pointer to density matrix (LCAO only)
     * @param pelec_in Pointer to electronic state (for charge, weights, ekb)
     * @param pw_wfc_in PW basis for wavefunction storage (PW only)
     */
  void init_sc(double sc_thr_in,
               int nsc_in,
               int nsc_min_in,
               double alpha_trial_in,
               double sccut_in,
               double sc_drop_thr_in,
               const UnitCell& ucell,
               bool direction_only_in,
               Parallel_Orbitals* ParaV_in,
               int nspin_in,
               const K_Vectors& kv_in,
               void* p_hamilt_in,
               void* psi_in,
#ifdef __LCAO
			   elecstate::DensityMatrix<TK, double> *dm_in, // mohan add 2025-11-02
#endif
			   elecstate::ElecState* pelec_in,
               ModulePW::PW_Basis_K* pw_wfc_in = nullptr);

  /**
   * @brief Calculate atomic magnetic moments using real-space projection (LCAO basis).
   *
   * @details Uses the DeltaSpin operator to compute magnetic moments from the density
   * matrix. For nspin=2, extracts only the z-component. For nspin=4, extracts
   * all three components from the interleaved 4-component spinor density matrix.
   * The moments are stored in state_.Mi_ (indexed by global atom index iat).
   *
   * @param step Current SCF iteration number (for logging)
   * @param print Whether to print moments to ofs_running
   */
  void cal_mi_lcao(const int& step, bool print = false);

  // The PW-basis magnetic-moment path (cal_mi_pw) has been lifted to the
  // free function spinconstrain::pw::cal_mi_pw() in deltaspin_pw_mi.h.

  /**
   * @brief Core workflow: apply lambda -> solve Hamiltonian -> compute magnetic moments.
   *
   * @details This is the central function called repeatedly during lambda optimization.
   *
   * @par LCAO path:
   *   1. Update lambda in DeltaSpin operator
   *   2. Solve HSolverLCAO (diagonalize without charge update)
   *   3. Calculate weights from eigenvalues
   *   4. Call cal_mi_lcao() to compute moments
   *
   * @par PW path:
   *   1. [First call only] Save subspace H, S, and becp from Hamiltonian
   *   2. Apply DeltaSpin correction via calculate_delta_hcc()
   *   3. Diagonalize in subspace via diag_responce(), update becp
   *   4. Calculate weights from new eigenvalues
   *   5. Call accumulate_Mi_from_becp() for each k-point
   *   6. MPI reduce Mi across k-pools
   *
   * @param i_step Current inner lambda step (-1 = initialization, 0+ = optimization)
   * @param delta_lambda Change in lambda from previous step (for incremental H correction)
   */
  void cal_mw_from_lambda(int i_step,
		  const ModuleBase::Vector3<double>* delta_lambda = nullptr);

  /**
   * @brief Calculate the spin constraint energy contribution: E_scon = -sum(lambda_i . Mi_i).
   *
   * @details Returns 0.0 if magnetic moments are not yet converged, because the
   * constraint energy is only physically meaningful when Mi ≈ M_target.
   * This energy is added to the total DFT energy in the SCF loop.
   *
   * @return Constraint energy in Ry (0.0 if not converged)
   */
  double cal_escon() { return state_.cal_escon(); }

  /// @brief Get the cached constraint energy from the last cal_escon() call (Ry)
  double get_escon() const { return state_.get_escon(); }

  /**
   * @brief Main lambda optimization loop using conjugate-gradient-like scheme.
   *
   * @details Iteratively adjusts Lagrange multipliers (lambda) to drive atomic
   * magnetic moments (Mi) toward target values. Uses:
   * - Polak-Ribiere formula for beta (conjugate direction)
   * - Linear interpolation for optimal step size (alpha_opt)
   * - Adaptive alpha_trial adjustment based on convergence behavior
   * - Gradient decay check for early termination
   *
   * @param outer_step Current SCF outer iteration number
   * @param rerun If true, use full PW solver for final charge update
   */
  void run_lambda_loop(int outer_step,
		  bool rerun,
		  std::ostream& ofs_running);

  /// @brief RMS error of the most recent lambda optimization loop (-1.0 if none has run).
  double get_last_rms_error() const { return last_rms_error_; }

  /**
   * @brief Alternative mode: sweep lambda values linearly for energy landscape mapping.
   *
   * @details Used for debugging or plotting E(lambda) curves. Scans from
   * sc_scan_lambda_start to sc_scan_lambda_end in sc_scan_steps steps.
   * Results written to lambda_scan_results.dat.
   *
   * @param outer_step Current SCF outer iteration number
   */
  void run_lambda_linear_scan(int outer_step, std::ostream& ofs_running);

  /// @brief Reset DeltaSpin operator initialization state when constraints change
  void reset_dspin_operator();

  /**
   * @brief Update wavefunctions and charge density after lambda optimization.
   *
   * @details Dispatcher to LCAO or PW (CPU/GPU) update paths.
   * For PW: performs subspace diagonalization + optional full-space refinement.
   *
   * @param delta_lambda Lambda change for incremental H correction
   * @param pw_solve If true, run full PW solver for refinement; if false, just update weights
   * @param full_update If true, apply full lambda (not delta) to H correction
   */
  void update_psi_charge(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve = true, bool full_update = false);

  // The PW-basis update implementation (update_psi_charge_pw_cpu/gpu) and the
  // subspace Hamiltonian correction (calculate_delta_hcc) have been lifted to
  // free functions spinconstrain::pw::update_psi_charge_pw_{cpu,gpu}() and
  // spinconstrain::pw::calculate_delta_hcc() in deltaspin_pw_mi.h.
  // (The old declaration update_psi_charge_pw() never had a definition.)

#ifdef __LCAO
  /// LCAO magnetic-moment helpers (orbital-matrix and mu*dm paths) have been
  /// lifted to free functions in deltaspin_lcao_mi.h (namespace spinconstrain::lcao).
#endif

  /// Lambda loop helpers (print_rms_stop, check_restriction, check_gradient_decay,
  /// cal_alpha_opt, print_header, print_termination) have been lifted to free
  /// functions in lambda_loop_helper.h. The class now only carries state and
  /// the core lambda-loop driver (run_lambda_loop).
  /// print_Mi and print_Mag_Force have also been lifted to lambda_loop_helper.h.

  /// @brief Use full PW solver (rerun) for higher precision in lambda loop
  bool higher_mag_prec = false;

public:
    /**
     * =============================================================
     * EXTERNAL POINTERS - Set by init_sc(), used throughout the module
     * =============================================================
     *
     * @par Design rationale for void* pointers
     * The Hamiltonian and Psi types differ between LCAO and PW bases.
     * Using void* avoids template coupling and allows the same SpinConstrain
     * code to work with both basis sets. Concrete types are recovered
     * via static_cast at call sites.
     */

    /// @brief Parallel orbitals distribution (row/col mapping for ScaLAPACK)
    Parallel_Orbitals *ParaV = nullptr;
    //--------------------------------------------------------------------------------
    // Pointers to external objects: Hamiltonian, wavefunctions, electronic state
    // These are type-erased void* to avoid coupling with specific Hamilt/Psi types
    void* p_hamilt = nullptr;     ///< Pointer to HamiltLCAO or HamiltPW
    void* psi = nullptr;          ///< Pointer to Psi<TK> wavefunction container
    elecstate::ElecState* pelec = nullptr;  ///< Electronic state: ekb, wg, charge, klist
    ModulePW::PW_Basis_K* pw_wfc_ = nullptr; ///< PW basis for wavefunction storage (PW only)
#ifdef __LCAO
    elecstate::DensityMatrix<TK, double>* dm_; ///< Density matrix pointer (LCAO only)
#endif
    const double meV_to_Ry = 7.349864435130999e-05; ///< Conversion factor
    K_Vectors kv_; ///< K-point vector list
    //--------------------------------------------------------------------------------

    /// Constraint parameters and runtime state (lambda, Mi, target_mag, indexing)
    ScState state_;

  public:
    /**
     * pubic methods for setting and getting spin-constrained DFT parameters
     * (thin forwarding shells to state_; kept for backward compatibility)
    */
    /// Public method to access the Singleton instance
    static SpinConstrain& getScInstance();
    /// Delete copy and move constructors and assign operators
    SpinConstrain(SpinConstrain const&) = delete;
    SpinConstrain(SpinConstrain&&) = delete;
    /// set element index to atom index map
    void set_atomCounts(const std::map<int, int>& atomCounts_in) { state_.set_atomCounts(atomCounts_in); }
    /// get element index to atom index map
    const std::map<int, int>& get_atomCounts() const { return state_.get_atomCounts(); }
    /// set element index to orbital index map
    void set_orbitalCounts(const std::map<int, int>& orbitalCounts_in) { state_.set_orbitalCounts(orbitalCounts_in); }
    /// get element index to orbital index map
    const std::map<int, int>& get_orbitalCounts() const { return state_.get_orbitalCounts(); }
    /// set lnchiCounts
    void set_lnchiCounts(const std::map<int, std::map<int, int>>& lnchiCounts_in) { state_.set_lnchiCounts(lnchiCounts_in); }
    /// get lnchiCounts
    const std::map<int, std::map<int, int>>& get_lnchiCounts() const { return state_.get_lnchiCounts(); }
    /// set sc_lambda
    void set_sc_lambda() { state_.set_sc_lambda(); }
    /// set sc_lambda from variable
    void set_sc_lambda(const ModuleBase::Vector3<double>* lambda_in, int nat_in) { state_.set_sc_lambda(lambda_in, nat_in); }
    /// set target_mag
    void set_target_mag() { state_.set_target_mag(); }
    /// set target_mag from variable
    void set_target_mag(const ModuleBase::Vector3<double>* target_mag_in, int nat_in) { state_.set_target_mag(target_mag_in, nat_in); }
    /// set target magnetic moment
    void set_target_mag(const std::vector<ModuleBase::Vector3<double>>& target_mag_in) { state_.set_target_mag(target_mag_in); }
    /// set constrain
    void set_constrain() { state_.set_constrain(); }
    /// set constrain from variable
    void set_constrain(const ModuleBase::Vector3<int>* constrain_in, int nat_in) { state_.set_constrain(constrain_in, nat_in); }
    /// get sc_lambda
    const std::vector<ModuleBase::Vector3<double>>& get_sc_lambda() const { return state_.get_sc_lambda(); }
    /// get target_mag
    const std::vector<ModuleBase::Vector3<double>>& get_target_mag() const { return state_.get_target_mag(); }
    /// get constrain
    const std::vector<ModuleBase::Vector3<int>>& get_constrain() const { return state_.get_constrain(); }
    /// get nat
    int get_nat() const { return state_.get_nat(); }
    /// get ntype
    int get_ntype() const { return state_.get_ntype(); }
    /// check atomCounts
    void check_atomCounts() const { state_.check_atomCounts(); }
    /// get iat
    int get_iat(int itype, int atom_index) { return state_.get_iat(itype, atom_index); }
    /// set nspin
    void set_nspin(int nspin) { state_.set_nspin(nspin); }
    /// get nspin
    int get_nspin() const { return state_.get_nspin(); }
    /// zero atomic magnetic moment
    void zero_Mi() { state_.zero_Mi(); }
    /// get decay_grad
    double get_decay_grad(int itype) const { return state_.get_decay_grad(itype); }
    /// set decay_grad
    void set_decay_grad() { state_.set_decay_grad(); }
    /// get decay_grad
    const std::vector<double>& get_decay_grad() const { return state_.get_decay_grad(); }
    /// set decay_grad from variable
    void set_decay_grad(const double* decay_grad_in, int ntype_in) { state_.set_decay_grad(decay_grad_in, ntype_in); }
    /// set decay grad switch
    void set_sc_drop_thr(double sc_drop_thr_in) { state_.set_sc_drop_thr(sc_drop_thr_in); }
    /// set input parameters
    void set_input_parameters(double sc_thr_in,
                              int nsc_in,
                              int nsc_min_in,
                              double alpha_trial_in,
                              double sccut_in,
                              double sc_drop_thr_in)
    { state_.set_input_parameters(sc_thr_in, nsc_in, nsc_min_in, alpha_trial_in, sccut_in, sc_drop_thr_in); }
    /// get sc_thr
    double get_sc_thr() const { return state_.get_sc_thr(); }
    /// get current adaptive sc threshold (max(initial_rms * sc_drop_thr_, sc_thr_))
    double get_current_sc_thr() const { return state_.get_current_sc_thr(); }
    /// get nsc
    int get_nsc() const { return state_.get_nsc(); }
    /// get nsc_min
    int get_nsc_min() const { return state_.get_nsc_min(); }
    /// get alpha_trial
    double get_alpha_trial() const { return state_.get_alpha_trial(); }
    /// get sccut
    double get_sccut() const { return state_.get_sccut(); }
    /// get sc_drop_thr
    double get_sc_drop_thr() const { return state_.get_sc_drop_thr(); }
    /// get computed magnetic moments Mi per atom
    const std::vector<ModuleBase::Vector3<double>>& get_Mi() const { return state_.get_Mi(); }
    /// get human-readable atom labels ("Fe_0", "Fe_1", ...) for table printing
    const std::vector<std::string>& get_atomLabels() const { return state_.get_atomLabels(); }
    /// @brief set orbital parallel info
    void set_ParaV(Parallel_Orbitals* ParaV_in);
    /// @brief set parameters for solver
    void set_solver_parameters(const K_Vectors& kv_in,
                               void* p_hamilt_in,
                               void* psi_in,
                               elecstate::ElecState* pelec_in);

  private:
    SpinConstrain(){};
    // Subspace buffers are owned by `pw_cache_` (RAII via release_cpu/gpu in the
    // PW update paths). The destructor is trivial; the singleton lives for the
    // whole program and the cache is released by the PW update functions.
    ~SpinConstrain() = default;
    SpinConstrain& operator=(SpinConstrain const&) = delete;  ///< Copy assignment deleted
    SpinConstrain& operator=(SpinConstrain &&) = delete;      ///< Move assignment deleted

  public:
    /// @brief Set DeltaSpin operator pointer for magnetic moment calculation (LCAO)
    /// @param op_in Base pointer, actual type is DeltaSpin<OperatorLCAO<TK, TR>>*
    void set_operator(hamilt::Operator<TK>* op_in);
    /// @brief Set magnetic moment convergence flag
    void set_mag_converged(bool is_Mi_converged_in) { state_.set_mag_converged(is_Mi_converged_in); }
    /// @brief Get magnetic moment convergence flag
    bool mag_converged() const { return state_.mag_converged(); }
    void set_npol(int npol) { state_.set_npol(npol); }
    int get_npol() const { return state_.get_npol(); }
    int get_nw() const { return state_.get_nw(); } ///< Total number of orbitals across all constrained atoms
    int get_iwt(int itype, int iat, int orbital_index) const { return state_.get_iwt(itype, iat, orbital_index); } ///< Convert (itype, iat, iw) to global orbital index
    /// @brief Get spin sign for k-point ik: +1 for spin-up, -1 for spin-down (nspin=2 only)
    int get_spin_sign(int ik) const;
  private:
    /// DeltaSpin operator pointer for LCAO magnetic moment calculation
    hamilt::Operator<TK>* p_operator = nullptr;

    /**
     * =============================================================
     * SUBSPACE DATA CACHING (PW basis only)
     * =============================================================
     *
     * @par Purpose
     * In the PW basis, the subspace Hamiltonian H_sub = <psi|H|psi> and
     * becp coefficients are expensive to compute. They are cached on the
     * first call to cal_mw_from_lambda() and reused across multiple lambda
     * steps within the same SCF iteration.
     *
     * @par Layout
     * - sub_h_save[ik * nbands * nbands + i * nbands + j]: H_sub for k-point ik
     * - sub_s_save: same layout for overlap matrix S_sub
     * - becp_save[ik * size_becp + ib * nkb * npol + ip]: becp coefficients
     * - lambda_in_sub_: lambda values at the time subspace data was saved,
     *   used to compute delta_lambda for incremental H corrections
     *
     * @par Memory management
     * Allocated with new[] on first cal_mw_from_lambda() call, freed in
     * update_psi_charge_pw_cpu/gpu() after final subspace diagonalization.
     */
  public:
    /// PW subspace data cache (H_sub/S_sub/becp + lambda snapshot). Owned object;
    /// buffers are device/host memory managed via allocate_cpu/gpu + release_cpu/gpu.
    pw::SubspaceCache pw_cache_;

  private:
    /// RMS error of the most recent lambda optimization loop; -1.0 if no loop has run.
    /// Used by ESolver to pass the current DeltaSpin RMS into the SCF iteration table.
    double last_rms_error_ = -1.0;

    /// Allow the lambda loop driver (still a member function) to record the RMS error.
    /// (kept private; run_lambda_loop is a member so it can write this directly)
};


} // namespace spinconstrain

#endif // SPIN_CONSTRAIN_H
