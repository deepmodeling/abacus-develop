/**
 * @file deltaspin_state.h
 * @brief Constraint parameters and runtime state for the DeltaSpin
 *        (spin-constrained DFT) module.
 *
 * @par Purpose
 * ScState owns all basis-set-independent data of the spin-constrained
 * calculation: per-atom Lagrange multipliers (lambda), target magnetic
 * moments, constraint flags, computed moments (Mi), atom/orbital indexing
 * maps, and the lambda-loop convergence parameters. It is deliberately
 * non-template: none of this data depends on the wavefunction type TK.
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
#ifndef DELTASPIN_STATE_H
#define DELTASPIN_STATE_H

#include <map>
#include <string>
#include <vector>

#include "source_base/vector3.h"

namespace spinconstrain
{

/**
 * @brief Per-atom spin constraint parameters parsed from STRU file.
 *
 * @details Stores the raw constraint data for a single atom before
 * it is distributed to the flat arrays (lambda_, target_mag_, constrain_).
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

class ScState
{
public:
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
    /// set sc_lambda from ScData (parsed from STRU file)
    void set_sc_lambda();
    /// set sc_lambda from variable
    void set_sc_lambda(const ModuleBase::Vector3<double>* lambda_in, int nat_in);
    /// set target_mag from ScData (parsed from STRU file)
    void set_target_mag();
    /// set target_mag from variable
    void set_target_mag(const ModuleBase::Vector3<double>* target_mag_in, int nat_in);
    /// set target magnetic moment
    void set_target_mag(const std::vector<ModuleBase::Vector3<double>>& target_mag_in);
    /// set constrain from ScData
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
    int get_nat() const;
    /// get ntype
    int get_ntype() const;
    /// check atomCounts
    void check_atomCounts() const;
    /// get iat
    int get_iat(int itype, int atom_index) const;
    /// set nspin
    void set_nspin(int nspin);
    /// get nspin
    int get_nspin() const;
    /// set npol
    void set_npol(int npol);
    /// get npol
    int get_npol() const;
    /// zero atomic magnetic moment
    void zero_Mi();
    /// get decay_grad (root process only: reads ScDecayGrad from json)
    double get_decay_grad(int itype) const;
    /// set decay_grad (zero-initialize per element type)
    void set_decay_grad();
    /// get decay_grad
    const std::vector<double>& get_decay_grad() const;
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
    /// get current adaptive sc threshold (max(initial_rms * sc_drop_thr_, sc_thr_))
    double get_current_sc_thr() const;
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
    /// get computed magnetic moments Mi per atom
    const std::vector<ModuleBase::Vector3<double>>& get_Mi() const;
    /// get human-readable atom labels ("Fe_0", "Fe_1", ...) for table printing
    const std::vector<std::string>& get_atomLabels() const;
    /// Total number of orbitals across all constrained atoms
    int get_nw() const;
    /// Convert (itype, iat, iw) to global orbital index
    int get_iwt(int itype, int iat, int orbital_index) const;
    /// Set magnetic moment convergence flag
    void set_mag_converged(bool is_Mi_converged_in) { this->is_Mi_converged = is_Mi_converged_in; }
    /// Get magnetic moment convergence flag
    bool mag_converged() const { return this->is_Mi_converged; }

    /**
     * @brief Calculate the spin constraint energy contribution: E_scon = -sum(lambda_i . Mi_i).
     * @return Constraint energy in Ry
     */
    double cal_escon();
    /// Get the cached constraint energy from the last cal_escon() call (Ry)
    double get_escon() const;

    /**
     * =============================================================
     * PUBLIC FIELDS - transitional, accessed directly by the lambda
     * loop and Mi accumulation code. Will be reduced to accessors in a
     * later refactoring step.
     * =============================================================
     */
    std::vector<ModuleBase::Vector3<double>> lambda_; ///< Lagrange multipliers (Ry/uB) per atom, 3 components
    std::vector<ModuleBase::Vector3<double>> target_mag_; ///< Target magnetic moments (uB) per atom
    std::vector<ModuleBase::Vector3<double>> Mi_; ///< Current computed magnetic moments (uB) per atom
    std::vector<ModuleBase::Vector3<int>> constrain_; ///< Per-atom/component constraint flags: 0=free, 1=constrained
    std::vector<std::string> atomLabels_; ///< Human-readable labels: "Fe_0", "Fe_1", etc.
    int nspin_ = 0; ///< Spin type: 2=collinear, 4=non-collinear
    int npol_ = 1; ///< Number of spinor components: 1 for nspin=2, 2 for nspin=4
    double sc_thr_ = 0.0; ///< Convergence threshold for RMS(Mi - M_target) in uB
    double sc_drop_thr_ = 1e-3; ///< Fraction of initial RMS for adaptive threshold
    double current_sc_thr_ = 0.0; ///< Adaptive threshold: max(initial_rms * sc_drop_thr_, sc_thr_)
    double alpha_trial_ = 0.0; ///< Initial trial step size (Ry/uB^2), adaptively adjusted during loop
    double restrict_current_ = 0.0; ///< Maximum allowed lambda change per step (Ry/uB)
    int nsc_ = 0; ///< Maximum number of inner lambda optimization steps
    int nsc_min_ = 0; ///< Minimum steps before early exit checks (gradient decay)
    double escon_ = 0.0; ///< Cached constraint energy from last cal_escon() call (Ry)
    bool is_Mi_converged = false; ///< Has the magnetic moment converged in the current SCF iteration?
    bool direction_only_ = false; ///< If true, only optimize spin direction
    double tpiba = 0.0; ///< 2*pi/a lattice constant scaling factor, saved from UnitCell

private:
    std::map<int, std::vector<ScAtomData>> ScData; ///< Raw constraint data indexed by element type (itype)
    std::map<int, double> ScDecayGrad; ///< Gradient decay thresholds (uB^2/eV) per element type
    std::vector<double> decay_grad_;   ///< Gradient decay thresholds converted to uB^2/Ry, per element type
    std::map<int, int> atomCounts;     ///< Number of atoms per element type: {itype -> nat_itype}
    std::map<int, int> orbitalCounts;  ///< Number of orbitals per element type: {itype -> nw_itype}
    std::map<int, std::map<int, int>> lnchiCounts; ///< {itype -> {L -> nchi}}: angular momentum channels
};

} // namespace spinconstrain

#endif // DELTASPIN_STATE_H
