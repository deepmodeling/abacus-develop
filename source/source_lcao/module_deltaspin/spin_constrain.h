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
 */
#ifndef SPIN_CONSTRAIN_H
#define SPIN_CONSTRAIN_H

#include <complex>
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

#ifdef __LCAO
#include "source_estate/module_dm/density_matrix.h" // mohan add 2025-11-02
#endif

namespace spinconstrain
{

/**
 * @brief Convert spinor occupation matrix to magnetic moment vector using Pauli matrices.
 *
 * @details For a two-component spinor wavefunction, the spin density matrix is:
 *   rho = |a|^2    a*b  |   = | (1+Mz)/2    (Mx+iMy)/2 |
 *         |b*a    |b|^2  |     | (Mx-iMy)/2   (1-Mz)/2  |
 * The magnetic moment components are extracted via Pauli matrix traces:
 *   Mx = Tr(rho * sigma_x) = occ[1] + occ[2]           (real part)
 *   My = Tr(rho * sigma_y) = Im(occ[1] - occ[2])        (imaginary part)
 *   Mz = Tr(rho * sigma_z) = occ[0] - occ[3]            (real part)
 * where occ = {|a|^2, a*b, b*a, |b|^2} from becp coefficients.
 *
 * @param occ 4-element array of occupation matrix elements (complex)
 * @param weight k-point weight for integration
 * @return 3D magnetic moment vector (Mx, My, Mz) in Bohr magnetons
 */
inline ModuleBase::Vector3<double> pauli_to_moment(const std::complex<double> occ[4], double weight)
{
    return ModuleBase::Vector3<double>(
        weight * (occ[1] + occ[2]).real(),
        weight * (occ[1] - occ[2]).imag(),
        weight * (occ[0] - occ[3]).real()
    );
}

struct ScAtomData;

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
   * The moments are stored in Mi_ (indexed by global atom index iat).
   *
   * @param step Current SCF iteration number (for logging)
   * @param print Whether to print moments to ofs_running
   */
  void cal_mi_lcao(const int& step, bool print = false);

  /**
   * @brief Calculate atomic magnetic moments using projector overlap (PW basis).
   *
   * @details For each k-point:
   *   1. Call OnsiteProjector::tabulate_atomic() to set up atomic projectors
   *   2. Call OnsiteProjector::overlap_proj_psi() to compute becp = <alpha|psi>
   *   3. Call accumulate_Mi_from_becp() to decompose becp into magnetic moments
   * Finally, sum Mi across all MPI k-pool ranks via Parallel_Reduce.
   */
  void cal_mi_pw();

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
  double cal_escon();

  /// @brief Get the cached constraint energy from the last cal_escon() call (Ry)
  double get_escon() const;

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
		  bool rerun = true);

  /**
   * @brief Alternative mode: sweep lambda values linearly for energy landscape mapping.
   *
   * @details Used for debugging or plotting E(lambda) curves. Scans from
   * sc_scan_lambda_start to sc_scan_lambda_end in sc_scan_steps steps.
   * Results written to lambda_scan_results.dat.
   *
   * @param outer_step Current SCF outer iteration number
   */
  void run_lambda_linear_scan(int outer_step);

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

  /**
   * @brief Wavefunction and charge density update implementation for PW basis.
   * @details Two-stage process:
   *          1. Subspace diagonalization: apply DeltaSpin correction and solve for each k-point
   *          2. Charge update: full-space diagonalization or direct charge update based on pw_solve
   */
  void update_psi_charge_pw(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update = false);
  
  /// CPU implementation of PW basis update
  void update_psi_charge_pw_cpu(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update = false);
  
#if ((defined __CUDA) || (defined __ROCM))
  /// GPU implementation of PW basis update
  void update_psi_charge_pw_gpu(const ModuleBase::Vector3<double>* delta_lambda, bool pw_solve, bool full_update = false);
#endif

  /**
   * @brief Compute DeltaSpin correction to the subspace Hamiltonian.
   *
   * @details Adds the constraint term to the Hamiltonian in the subspace:
   *   H_corrected = H_original + becp^† * delta_lambda * becp
   * For npol=2 (nspin=4), uses full 2x2 Pauli matrix coefficients:
   *   coeff = | lambda_z      lambda_x + i*lambda_y |
   *           | lambda_x - i*lambda_y   -lambda_z   |
   * For npol=1 (nspin=2), only the z-component with spin_sign.
   *
   * @param h_tmp Subspace Hamiltonian (nbands x nbands, in/out)
   * @param becp_k Projector coefficients for k-point ik
   * @param delta_lambda Lambda change per atom (or full lambda if full_update)
   * @param nbands Number of bands
   * @param nkb Total number of projectors
   * @param nh_iat Number of projectors per atom
   * @param ik K-point index
   * @param full_update If true, compute delta = lambda_current - lambda_at_save
   */
  void calculate_delta_hcc(std::complex<double>* h_tmp,
		  const std::complex<double>* becp_k,
		  const ModuleBase::Vector3<double>* delta_lambda,
		  const int nbands, const int nkb, const int* nh_iat, const int ik,
		  bool full_update = false);

#ifdef __LCAO
  /// @brief Convert orbital matrix to nested vector format [nspin][iat][iw]
  std::vector<std::vector<std::vector<double>>> convert(const ModuleBase::matrix& orbMulP);
  /// @brief Calculate magnetic moment from orbital matrix (LCAO alternative path)
  void calculate_MW(const std::vector<std::vector<std::vector<double>>>& AorbMulP);
  /// @brief Collect magnetic moment contributions from complex matrix mu*dm
  void collect_MW(ModuleBase::matrix& MecMulP,
                  const ModuleBase::ComplexMatrix& mud,
                  int nw,
                  int isk);
#endif

  /// Lambda loop helper: check if RMS error below threshold or max steps reached
  bool check_rms_stop(int outer_step, int i_step, double rms_error, double duration, double total_duration);

  /// Lambda loop helper: cap step size via restrict_current_ to prevent overshooting
  void check_restriction(const std::vector<ModuleBase::Vector3<double>>& search, double& alpha_trial);

  /**
   * @brief Lambda loop helper: check if dM/dlambda gradient has decayed below threshold.
   *
   * @details Computes the diagonal of the susceptibility matrix dM/dlambda for each
   * atom type. If max gradient < decay_grad[itype], the lambda optimization has
   * reached diminishing returns and should stop.
   *
   * @return true if gradient decayed below threshold, false otherwise
   */
  bool check_gradient_decay(std::vector<ModuleBase::Vector3<double>> new_spin,
                            std::vector<ModuleBase::Vector3<double>> old_spin,
                            std::vector<ModuleBase::Vector3<double>> new_delta_lambda,
                            std::vector<ModuleBase::Vector3<double>> old_delta_lambda,
                            bool print = false);
  /// @brief Lambda loop helper: calculate optimal step size via linear interpolation
  double cal_alpha_opt(std::vector<ModuleBase::Vector3<double>> spin,
                       std::vector<ModuleBase::Vector3<double>> spin_plus,
                       const double alpha_trial);
  /// Print header at start of lambda loop
  void print_header();
  /// Print termination message with final spin and lambda values
  void print_termination();

  /// Print magnetic moments to output stream
  void print_Mi(std::ofstream& ofs_running);

  /// Print magnetic force (defined as dL/dMi = -lambda[iat]) in eV/uB
  void print_Mag_Force(std::ofstream& ofs_running);

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
    double tpiba = 0.0; /// @brief 2*pi/a lattice constant scaling factor, saved from UnitCell
    const double meV_to_Ry = 7.349864435130999e-05; ///< Conversion factor
    K_Vectors kv_; ///< K-point vector list
    //--------------------------------------------------------------------------------

  public:
    /**
     * pubic methods for setting and getting spin-constrained DFT parameters
    */
    /// Public method to access the Singleton instance
    static SpinConstrain& getScInstance();
    /// Delete copy and move constructors and assign operators
    SpinConstrain(SpinConstrain const&) = delete;
    SpinConstrain(SpinConstrain&&) = delete;
    /// set element index to atom index map
    void set_atomCounts(const std::map<int, int>& atomCounts_in);
    /// get element index to atom index map
    const std::map<int, int>& get_atomCounts() const;
    /// set element index to orbital index map
    void set_orbitalCounts(const std::map<int, int>& orbitalCounts_in);
    /// get element index to orbital index map
    const std::map<int, int>& get_orbitalCounts() const;
    /// set lnchiCounts
    void set_lnchiCounts(const std::map<int, std::map<int, int>>& lnchiCounts_in);
    /// get lnchiCounts
    const std::map<int, std::map<int, int>>& get_lnchiCounts() const;
    /// set sc_lambda
    void set_sc_lambda();
    /// set sc_lambda from variable
    void set_sc_lambda(const ModuleBase::Vector3<double>* lambda_in, int nat_in);
    /// set target_mag
    void set_target_mag();
    /// set target_mag from variable
    void set_target_mag(const ModuleBase::Vector3<double>* target_mag_in, int nat_in);
    /// set target magnetic moment
    void set_target_mag(const std::vector<ModuleBase::Vector3<double>>& target_mag_in);
    /// set constrain
    void set_constrain();
    /// set constrain from variable
    void set_constrain(const ModuleBase::Vector3<int>* constrain_in, int nat_in);
    /// get sc_lambda
    const std::vector<ModuleBase::Vector3<double>>& get_sc_lambda() const;
    /// get target_mag
    const std::vector<ModuleBase::Vector3<double>>& get_target_mag() const;
    /// get constrain
    const std::vector<ModuleBase::Vector3<int>>& get_constrain() const;
    /// get nat
    int get_nat();
    /// get ntype
    int get_ntype();
    /// check atomCounts
    void check_atomCounts();
    /// get iat
    int get_iat(int itype, int atom_index);
    /// set nspin
    void set_nspin(int nspin);
    /// get nspin
    int get_nspin() const;
    /// zero atomic magnetic moment
    void zero_Mi();
    /// get decay_grad
    double get_decay_grad(int itype);
    /// set decay_grad
    void set_decay_grad();
    /// get decay_grad
    const std::vector<double>& get_decay_grad();
    /// set decay_grad from variable
    void set_decay_grad(const double* decay_grad_in, int ntype_in);
    /// set decay grad switch
    void set_sc_drop_thr(double sc_drop_thr_in);
    /// set input parameters
    void set_input_parameters(double sc_thr_in,
                              int nsc_in,
                              int nsc_min_in,
                              double alpha_trial_in,
                              double sccut_in,
                              double sc_drop_thr_in);
    /// get sc_thr
    double get_sc_thr() const;
    /// get nsc
    int get_nsc() const;
    /// get nsc_min
    int get_nsc_min() const;
    /// get alpha_trial
    double get_alpha_trial() const;
    /// get sccut
    double get_sccut() const;
    /// get sc_drop_thr
    double get_sc_drop_thr() const;
    /// @brief set orbital parallel info
    void set_ParaV(Parallel_Orbitals* ParaV_in);
    /// @brief set parameters for solver
    void set_solver_parameters(const K_Vectors& kv_in,
                               void* p_hamilt_in,
                               void* psi_in,
                               elecstate::ElecState* pelec_in);

  private:
    /**
     * =============================================================
     * PRIVATE DATA MEMBERS - Internal state of SpinConstrain
     * =============================================================
     *
     * @par Unit conversion
     * - lambda_: Ry/uB internally, but meV/uB in input file (STRU)
     * - target_mag_, Mi_: uB (Bohr magnetons)
     * - alpha_trial_: Ry/uB^2 internally, but input is eV/uB^2
     * - restrict_current_: Ry/uB internally, but input is eV/uB
     * - decay_grad_: uB^2/Ry internally, but uB^2/eV in ScDecayGrad
     *
     * @par Indexing
     * All per-atom arrays (lambda_, target_mag_, Mi_, constrain_) are indexed
     * by GLOBAL atom index (iat), which runs from 0 to nat-1. The mapping
     * from (element_type, local_atom_index) to iat is handled by get_iat().
     */
    SpinConstrain(){};                               ///< Private constructor (Singleton)
    ~SpinConstrain()
    {
        delete[] sub_h_save;
        delete[] sub_s_save;
        delete[] becp_save;
        sub_h_save = nullptr;
        sub_s_save = nullptr;
        becp_save = nullptr;
    };
    SpinConstrain& operator=(SpinConstrain const&) = delete;  ///< Copy assignment deleted
    SpinConstrain& operator=(SpinConstrain &&) = delete;      ///< Move assignment deleted
    std::map<int, std::vector<ScAtomData>> ScData; ///< Raw constraint data indexed by element type (itype)
    std::map<int, double> ScDecayGrad; ///< Gradient decay thresholds (uB^2/eV) per element type
    std::vector<double> decay_grad_;   ///< Gradient decay thresholds converted to uB^2/Ry, per element type
    std::map<int, int> atomCounts;     ///< Number of atoms per element type: {itype -> nat_itype}
    std::map<int, int> orbitalCounts;  ///< Number of orbitals per element type: {itype -> nw_itype}
    std::map<int, std::map<int, int>> lnchiCounts; ///< {itype -> {L -> nchi}}: angular momentum channels
    std::vector<ModuleBase::Vector3<double>> lambda_; ///< Lagrange multipliers (Ry/uB) per atom, 3 components
    std::vector<ModuleBase::Vector3<double>> target_mag_; ///< Target magnetic moments (uB) per atom
    std::vector<ModuleBase::Vector3<double>> Mi_; ///< Current computed magnetic moments (uB) per atom
    std::vector<std::string> atomLabels_; ///< Human-readable labels: "Fe_0", "Fe_1", etc.
    double escon_ = 0.0; ///< Cached constraint energy from last cal_escon() call (Ry)
    int nspin_ = 0; ///< Spin type: 2=collinear, 4=non-collinear
    int npol_ = 1; ///< Number of spinor components: 1 for nspin=2, 2 for nspin=4
    /**
     * =============================================================
     * LAMBDA LOOP INPUT PARAMETERS
     * =============================================================
     */
    int nsc_; ///< Maximum number of inner lambda optimization steps
    int nsc_min_; ///< Minimum steps before early exit checks (gradient decay)
    double sc_drop_thr_ = 1e-3; ///< Fraction of initial RMS for adaptive threshold
    double sc_thr_; ///< Convergence threshold for RMS(Mi - M_target) in uB
    double current_sc_thr_; ///< Adaptive threshold: max(initial_rms * sc_drop_thr_, sc_thr_)
    std::vector<ModuleBase::Vector3<int>> constrain_; ///< Per-atom/component constraint flags: 0=free, 1=constrained
    bool debug = false; ///< Debug flag for verbose output
    double alpha_trial_; ///< Initial trial step size (Ry/uB^2), adaptively adjusted during loop
    double restrict_current_; ///< Maximum allowed lambda change per step (Ry/uB), prevents overshooting
    bool direction_only_ = false; ///< If true, only optimize spin direction (project out parallel lambda component)

  public:
    /// @brief Set DeltaSpin operator pointer for magnetic moment calculation (LCAO)
    /// @param op_in Base pointer, actual type is DeltaSpin<OperatorLCAO<TK, TR>>*
    void set_operator(hamilt::Operator<TK>* op_in);
    /// @brief Set magnetic moment convergence flag
    void set_mag_converged(bool is_Mi_converged_in){this->is_Mi_converged = is_Mi_converged_in;}
    /// @brief Get magnetic moment convergence flag
    bool mag_converged() const {return this->is_Mi_converged;}
    void set_npol(int npol);
    int get_npol() const;
    int get_nw() const; ///< Total number of orbitals across all constrained atoms
    int get_iwt(int itype, int iat, int orbital_index) const; ///< Convert (itype, iat, iw) to global orbital index
    /// @brief Get spin sign for k-point ik: +1 for spin-up, -1 for spin-down (nspin=2 only)
    int get_spin_sign(int ik) const;
    /**
     * @brief Accumulate magnetic moments from becp coefficients for a single k-point.
     *
     * @details For npol=2 (nspin=4), computes full Pauli decomposition:
     *   occ[0] = sum(becp_up^* * becp_up), occ[1] = sum(becp_up^* * becp_dn),
     *   occ[2] = sum(becp_dn^* * becp_up), occ[3] = sum(becp_dn^* * becp_dn)
     *   Mi = pauli_to_moment(occ, weight)
     * For npol=1 (nspin=2), only z-component:
     *   occ = sum(|becp|^2), Mi.z += weight * occ * spin_sign
     *
     * @param becp Projector coefficients <alpha_{l,m}|psi_{k,i}>
     * @param nkb Total number of projectors
     * @param nbands Number of bands
     * @param npol Number of spinor components
     * @param ik K-point index (for spin_sign lookup in nspin=2)
     * @param wg_ik Band weights for this k-point
     * @param nh_iat Number of projectors per atom
     */
    void accumulate_Mi_from_becp(const std::complex<double>* becp,
                                 int nkb,
                                 int nbands,
                                 int npol,
                                 int ik,
                                 const double* wg_ik,
                                 const int* nh_iat);
  private:
    /// DeltaSpin operator pointer for LCAO magnetic moment calculation
    hamilt::Operator<TK>* p_operator = nullptr;
    /// @brief Flag: has the magnetic moment converged in the current SCF iteration?
    bool is_Mi_converged = false;

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
    TK* sub_h_save = nullptr;       ///< Cached subspace Hamiltonian for all k-points
    TK* sub_s_save = nullptr;       ///< Cached subspace overlap matrix for all k-points
    TK* becp_save = nullptr;        ///< Cached becp coefficients for all k-points
    std::vector<ModuleBase::Vector3<double>> lambda_in_sub_; ///< Lambda values when subspace was saved
};


/**
 * @brief Per-atom spin constraint parameters parsed from STRU file.
 *
 * @details Stores the raw constraint data for a single atom before
 * it is distributed to the flat arrays (lambda_, target_mag_, constrain_).
 * The constraint data is organized by element type (itype) in the ScData map.
 *
 * @par Target moment specification (mag_type):
 * - mag_type=0: Direct Cartesian components (mx, my, mz) in uB
 * - mag_type=1: Spherical coordinates (magnitude, theta, phi)
 *   - target_mag_val: |M| in uB
 *   - target_mag_angle1: polar angle theta (degrees) from z-axis
 *   - target_mag_angle2: azimuthal angle phi (degrees) in xy-plane
 *   Conversion: Mx = |M|*sin(theta)*cos(phi), My = |M|*sin(theta)*sin(phi), Mz = |M|*cos(theta)
 */
struct ScAtomData {
    int index;                              ///< Local atom index within its element type
    std::vector<double> lambda;             ///< Initial lambda values (Ry/uB), 3 components (x,y,z)
    std::vector<double> target_mag;         ///< Target magnetic moment (uB), 3 components
    std::vector<int> constrain;             ///< Constraint flags: 0=free, 1=constrained, per component
    int mag_type;                           ///< 0=Cartesian (mx,my,mz), 1=spherical (|M|,theta,phi)
    double target_mag_val;                  ///< For mag_type=1: target moment magnitude (uB)
    double target_mag_angle1;               ///< For mag_type=1: polar angle theta (degrees)
    double target_mag_angle2;               ///< For mag_type=1: azimuthal angle phi (degrees)
};

} // namespace spinconstrain

#endif // SPIN_CONSTRAIN_H
