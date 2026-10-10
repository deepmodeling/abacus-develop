#ifndef ABACUS_SOURCE_ESOLVER_ESOLVER_LR_LCAO_TDDFT_H
#define ABACUS_SOURCE_ESOLVER_ESOLVER_LR_LCAO_TDDFT_H

#include "source_esolver/esolver_fp.h"
#include "source_io/module_parameter/input_parameter.h"
#include "source_cell/unitcell.h"
#include "source_hamilt/hamilt.h"
#include "source_estate/elecstate.h"
#include "source_hamilt/hamilt.h"
#include "source_estate/elecstate_lcao.h"

#include <vector>   //future tensor
#include <memory>
#include <string>

#include "source_esolver/esolver_ks_lcao.h" //for the move constructor
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
#include "source_lcao/module_lr/hamilt_casida.h"
#include "source_lcao/module_lr/state_track.h"
#include "source_hamilt/module_gint/gint_info.h"
#include "source_estate/module_pot/potential_new.h"
#ifdef __EXX
// #include <RI/physics/Exx.h>
#include "source_lcao/module_ri/exx_lri.h"
#include "source_hamilt/module_xc/exx_info.h" // for Exx_Info value member
#endif
namespace LR { template <typename T> struct GradientInputs; }
namespace ModuleESolver
{
    ///Excited State Solver: Linear Response TDDFT (Tamm Dancoff Approximation) 
    template<typename T, typename TR = double>
    class ESolver_LR : public ModuleESolver::ESolver_FP
    {
    public:
        ESolver_LR(const Input_para& inp, const std::string& in_dir, const std::string& out_dir);
        ~ESolver_LR() {}

        ///input: input, call, basis(LCAO), psi(ground state), elecstate
        // initialize sth. independent of the ground state
        virtual void before_all_runners(BaseCell& basecell, const Input_para& inp) override;
        virtual void runner(BaseCell& basecell, int istep) override;
        virtual void after_all_runners(BaseCell& basecell) override;

        /// Total energy of the excited state being relaxed: E_gs + Omega. Zero outside a
        /// relaxation, where nothing consumes it and the old behaviour is kept.
        virtual double cal_energy() override;
        /// Force on the atoms in the relaxed excited state, F = -d(E_gs + Omega)/dR (Ry/Bohr).
        virtual void cal_force(BaseCell& basecell, ModuleBase::matrix& force) override;
        /// Not implemented: there is no excited-state stress. `cell-relax` is rejected at input.
        virtual void cal_stress(BaseCell& basecell, ModuleBase::matrix& stress) override;

      protected:
        const std::string in_dir;
        const std::string out_dir;
        const UnitCell* ucell_ = nullptr;
        std::vector<double> orb_cutoff_;

        /// cached aliases, read once here instead of at every log/print call site below
        std::ofstream& ofs_running_ = GlobalV::ofs_running;
        std::ofstream& ofs_warning_ = GlobalV::ofs_warning;
        const int my_rank_ = GlobalV::MY_RANK;

        /// @brief the ground-state solver, kept alive across ionic steps (esolver_type = "ks-lr").
        /// Null on the `lr` path, where the ground state comes from files instead.
        std::unique_ptr<ModuleESolver::ESolver_KS_LCAO<T, TR>> ks_;

        // Geometry-dependent objects that the ground-state solver already builds for the current
        // structure. They are aliased rather than rebuilt: `ESolver_KS_LCAO::before_scf` refreshes
        // its own copies every ionic step, so recomputing them here would be both wasted work and
        // a chance for the two to disagree. On the `lr` path the pointers are bound to this
        // object's own members below, which `initialize_from_unitcell_` fills from file.
        Grid_Driver gd_own_;                  ///< only used when `ks_` is null
        TwoCenterBundle two_center_bundle_own_;  ///< only used when `ks_` is null
        Grid_Driver* gd_ptr_ = nullptr;
        const TwoCenterBundle* tcb_ptr_ = nullptr;
        Parallel_Grid* pgrid_ptr_ = nullptr;
        Structure_Factor* sf_ptr_ = nullptr;
        pseudopot_cell_vl* locpp_ptr_ = nullptr;

        Grid_Driver& gd() const { return *this->gd_ptr_; }
        const TwoCenterBundle& tcb() const { return *this->tcb_ptr_; }
        Parallel_Grid& pgrid() const { return *this->pgrid_ptr_; }
        Structure_Factor& sfac() const { return *this->sf_ptr_; }
        pseudopot_cell_vl& vloc() const { return *this->locpp_ptr_; }
        /// bind the aliases above; `ks_` must already be set (or null for the `lr` path)
        void bind_ground_state_aliases_();

        // not to use ElecState because 2-particle state is quite different from 1-particle state.
        // implement a independent one (ExcitedState) to pack physical properties if needed.
        // put the components of ElecState here: 
        std::vector<std::shared_ptr<LR::PotHxcLR>> pot;

        // ground state info 

        /// @brief ground state wave function
        std::unique_ptr<psi::Psi<T>> psi_ks;  ///< KS orbitals used in the [nocc+nvirt] window
        /// @brief all KS orbitals. On the `ks-lr` path this aliases the ground-state solver's
        /// `psi` (nk x nbands x nbasis -- far too big to copy every ionic step, and only read
        /// here); on the `lr` path it points at `psi_ks_all_own_`, filled from file.
        psi::Psi<T>* psi_ks_all_ = nullptr;
        std::unique_ptr<psi::Psi<T>> psi_ks_all_own_;

        /// @brief ground state bands, read from the file, or moved from ESolver_FP::pelec.ekb
        ModuleBase::matrix eig_ks;///< ground state eigenvalues in the [nocc+nvirt] window
        ModuleBase::matrix eig_ks_all; ///< all eigenvalues of ground state, read from the file, or moved from ESolver_FP::pelec.ekb
        ModuleBase::matrix wg_ks;   /// occupation numbers of ground state in the [nocc+nvirt] window
        ModuleBase::matrix wg_ks_all;   /// occupation number of all bands of ground state


        // @brief only needed for force calculation 
        std::unique_ptr<elecstate::Potential> pot_gs;
        std::unique_ptr<elecstate::Potential> pot_gs_hartree;   /// ground-state Hartree potential, only used for test_force
        double etxc_gs = 0.;
        double vtxc_gs = 0.;

        std::shared_ptr<LR::PotHxcLR> pot_hxc_gs; /// used in lr-grad, in the ground-state Hxc gradient term coming from dF/dC

        /// @brief Excited state wavefunction (locc, lvirt are local size of nocc and nvirt in each process)
        /// size of X: [neq][{nstate, nloc_per_state}], namely:
        /// - [nspin][{nstates, nk* (locc* lvirt}] for close- shell,
        /// -  [1][{nstates, nk * (locc[0] * lvirt[0]) + nk * (locc[1] * lvirt[1])}] for open-shell
        std::vector<ct::Tensor> X;
        int nloc_per_state = 1;

        std::vector<int> nocc;   ///< number of occupied orbitals for each spin used in the calculation
        int nocc_in = 1;    ///< nocc read from input (adjusted by nelec): max(spin-up, spindown)
        int nocc_max = 1;   ///< full occupied count in the largest spin channel
        std::vector<int> nvirt;   ///< number of virtual orbitals for each spin used in the calculation
        int nvirt_in = 1;   ///< nvirt read from input (adjusted by nelec): min(spin-up, spindown)
        int nbands = 2;
        int nbasis = 2;
        /// n_occ*nvirt, the basis size of electron-hole pair representation
        std::vector<int> npairs;
        /// how many 2-particle states to be solved
        int nstates = 1;
        int nspin = 1;
        int nk = 1;
        int nupdown = 0;
        bool openshell = false;
        std::string xc_kernel;

        void initialize_from_unitcell_(UnitCell& ucell, const Input_para& inp);
        /// one-time setup from the ground-state solver; ends by calling `refresh_from_ks_`
        void initialize_from_ks_(UnitCell& ucell, const Input_para& inp);
        /// re-read everything that depends on the atomic positions, once per ionic step
        void refresh_from_ks_(UnitCell& ucell);
        bool ks_initialized_ = false;   ///< whether `initialize_from_ks_` has already run
        bool exx_owned_ = false;        ///< `exx_lri` was built here (so its Cs/Vs are ours to refresh)

        // ---------- geometry relaxation on an excited state ----------
        /// resolve and validate `lr_target_state` / `lr_target_spin`; call once the dimensions
        /// (and therefore `openshell`) are final
        void setup_relax_target_();
        bool excited_relax_ = false;   ///< driving `calculation = relax` from an excited state
        int target_is_ = 0;            ///< spin block of the relaxed state in `X` and `pelec->ekb`
        double etot_gs_ = 0.0;         ///< ground-state total energy of the current step (Ry)
        ModuleBase::matrix force_gs_;  ///< ground-state force of the current step (Ry/Bohr, F = -dE/dR)
        /// The LR part of the excited-state force, -d(Omega)/dR (Ry/Bohr). 
        ModuleBase::matrix lr_force_;
        /// The state currently being followed, as an index into `X` / `pelec->ekb`.
        /// Seeded from `lr_target_state` on the first ionic step, then re-chosen at every
        /// later step by cross-geometry group overlap, then the single/JT reference.
        /// Following a fixed INDEX instead is what makes a relaxation fail near a
        /// degeneracy: the index always names the n-th lowest root, so as soon as two
        /// surfaces cross, "the target" jumps to a different diabatic state and the force
        /// is discontinuous. CG assumes a conservative field and cannot recover from that.
        int target_state_ = -1;
        /// Previous ionic step's amplitude for the followed state (local part), the
        /// reference the overlap is taken against. Empty on the first step.
        std::vector<T> target_X_prev_;
        LR::RootBasis<T> target_basis_prev_;
        /// Match the old subspace to current energy groups and refresh the single/group references.
        /// `ofs` receives the note when the followed root changes index, and the warning when
        /// no current root resembles the previous one.
        void follow_target_state_(std::ofstream& ofs);

        /// index of the relaxed state inside `pelec->ekb`
        int target_ekb_offset_() const
        { return this->openshell ? this->target_state_
                                 : this->target_is_ * this->nstates + this->target_state_; }

        std::vector<std::string> spin_types;

        std::unique_ptr<ModuleGint::GintInfo> gint_info_ = nullptr;
        void set_gint();

        /// @brief variables for parallel distribution of KS orbitals
        Parallel_2D paraC_;
        /// @brief variables for parallel distribution of excited states
        std::vector<Parallel_2D> paraX_;

        // ---------------- the Z-vector (CPSCF) window ----------------
        // The Z-vector equation enforces the Brillouin condition in EVERY occupied-virtual
        // rotation, so it must not be confined to the `nvirt` window X lives in. 
        // It should be the whole AO virtual space.
        //
        // These mirror `psi_ks` / `eig_ks` / `paraC_` / `paraX_` / `nvirt` / `nloc_per_state`
        // but span every virtual band the ground state produced (the input `nbands`), so the
        // window is widened by raising *nbands*, not `nvirt`. X keeps its own window, so Omega
        // -- and with it any finite-difference reference -- is untouched.
        std::unique_ptr<psi::Psi<T>> psi_ks_z_;
        ModuleBase::matrix eig_ks_z_;
        Parallel_2D paraC_z_;
        std::vector<Parallel_2D> paraX_z_;
        std::vector<int> nvirt_z_;
        int nbands_z_ = 0;
        int nloc_per_state_z_ = 0;
        /// (re)build the Z window from `psi_ks_all_` / `eig_ks_all`. `desc_src` describes the
        /// source wavefunction's 2D layout; unused (and may be null) in a serial build.
        void fill_z_window_(const int* desc_src);
        /// Widen `nst` X blocks starting at `istate_begin` from the X window into the Z window,
        /// zero-filling the virtual rows X does not have.
        ct::Tensor pad_X_to_z_(const int ispin, const int istate_begin, const int nst) const;
        /// @brief variables for parallel distribution of matrix in AO representation
        Parallel_Orbitals paraMat_;
        Parallel_Orbitals paraMat_all_; // for the parallelized size of the KS orbitals


        LCAO_Orbitals orb_; ///< numerical atomic orbital data for single-point evaluation
        std::vector<std::complex<double>> velocity_mo; ///< store the velocity matrix elements in MO basis
        int cal_nupdown_form_occ(const ModuleBase::matrix& wg);
        void setup_2center_table(TwoCenterBundle& two_center_bundle, LCAO_Orbitals& orb, UnitCell& ucell);

        /// @brief allocate and set the inital value of X
        void setup_eigenvectors_X();
        void set_X_initial_guess();

        /// @brief read in the ground state wave function, band energy and occupation
        virtual void read_ks_wfc();
        /// @brief  read in the ground state charge density
        void read_ks_chg(Charge& chg);

        virtual void init_pot(const Charge& chg_gs);

        /// @brief check the legality of the input parameters
        void parameter_check() const;
        /// @brief set nocc, nvirt, nbasis, npairs and nstates
        void set_dimension();
        /// reset nocc, nvirt, npairs after read ground-state wavefunction when nspin=2
        void reset_dim_spin2();

        /// setup Parallel_Orbitals info. beyond Parallel_2D
        void set_parallel_orbitals_band(Parallel_Orbitals& p, const int nbands_in);

        ///========================== for gradient calculation =========================
        void init_pot_groundstate(const Charge& chg_gs);
        /// Solve the Z-vector equation for `nst` independent blocks of `Xz`.
        ct::Tensor solve_zvector_eqation(const int ispin, const int nst, const ct::Tensor& Xz);
        /// Excited-state gradients d(Omega)/dR, one matrix per state solved. `istate_only >= 0`
        /// restricts it to that one state: geometry relaxation follows a single state, and the
        /// Z-vector solve dominates the cost. -1 does all `nstates`.
        std::vector<ModuleBase::matrix> cal_force(const int ispin, const int istate_only = -1);
        /// open-shell (spin-unrestricted) excited-state force: X holds [up | down] and every
        /// density matrix has two independent channels
        std::vector<ModuleBase::matrix> cal_force_openshell(const int istate_only = -1);
        /// @brief Gradients for excitation vectors supplied by the caller, already widened into
        ///        the Z window -- the two functions above are thin wrappers that widen the stored
        ///        eigenvectors and look up their `omega`.
        ///
        /// The blocks of `Xz` need not be the eigenvectors the Casida diagonalizer returned. Any
        /// normalized vector inside a degenerate multiplet is an eigenvector with the same
        /// `omega`, so passing a linear combination is what turns the per-state gradient into the
        /// full degenerate-subspace gradient matrix; see `cal_grad_matrix_degenerate` and
        /// `grad_degen.h`.
        ///
        /// @param omega        excitation energy of each block (Ry); its size sets the block count
        /// @param label_begin  state index the first block is reported under (labels only)
        std::vector<ModuleBase::matrix> cal_force_Xz(const int ispin, const ct::Tensor& Xz,
            const std::vector<double>& omega, const int label_begin);
        /// open-shell counterpart of `cal_force_Xz`
        std::vector<ModuleBase::matrix> cal_force_openshell_Xz(const ct::Tensor& Xz,
            const std::vector<double>& omega, const int label_begin);
        /// @brief The linear vibronic coupling (LVC) data of one degenerate multiplet.
        ///
        /// $H(\delta R)=\Omega_0\mathbb{1}+\sum_{A\alpha}\delta R_{A\alpha}G^{(A\alpha)}$ is the
        /// complete first-order description of a degeneracy, and $G$ is its parameter set. Kept as
        /// one object rather than as loose force matrices because the excited-state relaxation and
        /// (later) non-adiabatic dynamics need exactly the same data: near a degeneracy the correct
        /// propagation is on this coupled $d\times d$ model, not on an adiabatic gradient.
        struct MultipletLVC
        {
            int ispin = 0;                  ///< which spin channel (index into `spin_types`)
            std::vector<int> states;        ///< the multiplet's state indices, ascending
            double omega0 = 0.0;            ///< the common excitation energy (Ry)
            double omega_spread = 0.0;      ///< max - min over the members (Ry); 0 if exact
            /// G[k][l], symmetric, each a (nat, 3) force matrix -- i.e. the 3N matrices of d x d
            std::vector<std::vector<ModuleBase::matrix>> g;
            int dim() const { return static_cast<int>(states.size()); }
            /// $\operatorname{Tr}G/d$: the one smooth, basis-independent 3N vector field the
            /// multiplet has. Following it preserves the symmetric configuration.
            ModuleBase::matrix average_force() const;
        };
        /// LVC data of every multiplet found at this geometry, rebuilt on each call of
        /// `cal_force_and_grad_matrix_`. Empty unless `lr_degen_thr > 0`.
        std::vector<MultipletLVC> multiplet_lvc_;
        /// The multiplet `lr_target_state` belongs to, or empty when the target is non-degenerate
        /// or `lr_degen_mode = state`. Refreshed every ionic step by
        /// `resolve_target_multiplet_`, and it is what makes `cal_energy` and the reported gradient
        /// describe the same surface.
        std::vector<int> target_group_;
        /// Fill `target_group_` from the excitation energies of the current geometry.
        void resolve_target_multiplet_();
        /// The excitation energy the relaxation is minimising: the target state's own, or the
        /// multiplet average when `target_group_` is set.
        double target_omega_() const;
        /// The LR half of the force for the current geometry, following whichever surface
        /// `lr_degen_mode` selects. `ofs` receives the note when that is not a single state.
        ModuleBase::matrix cal_lr_force_relax_(std::ofstream& ofs);
        /// @brief The force of the steepest-descending branch of the target multiplet, i.e.
        ///        `lr_degen_mode = jt`.
        ///
        /// Assembles the off-diagonal part of the gradient matrix (which `average` does not need),
        /// solves the joint direction/mixing optimization in `LR::find_jt_direction`, and returns
        /// that branch's own force. Once a step has split the multiplet there is no group left and
        /// the ordinary single-state path takes over, so the mode is self-limiting.
        ///
        /// @param diag  the multiplet's per-state forces, already computed
        ModuleBase::matrix cal_jt_force_(const std::vector<ModuleBase::matrix>& diag,
            std::ofstream& ofs);
        LR::GradientInputs<T> gradient_inputs_() const;
        /// Widen a multiplet's eigenvectors into the Z window, one block each. Members need not be
        /// contiguous, so they are padded one at a time.
        ct::Tensor pad_group_to_z_(const int ispin, const std::vector<int>& group) const;
        /// @brief Per-state gradients of every state, plus the gradient matrix of each degenerate
        ///        multiplet when `lr_degen_thr` asks for it. The single-point entry point.
        ///
        /// Fills `multiplet_lvc_`.
        void cal_force_and_grad_matrix_(const int ispin, std::ofstream& ofs);
        /// @brief The gradient matrix of one degenerate multiplet,
        ///        $G^{(A\alpha)}_{kl}=\langle X_k|\partial A/\partial R_{A\alpha}|X_l\rangle$.
        ///
        /// At a degeneracy no single state has a gradient vector -- the branch slopes along a
        /// displacement $u$ are the eigenvalues of $\sum_{A\alpha}u_{A\alpha}G^{(A\alpha)}$, whose
        /// eigenvectors depend on $u$ -- so this whole matrix, not its diagonal, is the first-order
        /// information. It is obtained from the polarization identity
        /// $G_{kl}=\mathcal F[(X_k{+}X_l)/\sqrt2]-\tfrac12(G_{kk}+G_{ll})$, which needs no new
        /// physics. `grad_degen.h` derives why that is exact.
        ///
        /// @param group  state indices of the multiplet, from `LR::group_degenerate_states`
        /// @param diag   their per-state gradients, i.e. $G_{kk}$, already computed by `cal_force`
        /// @param ofs    log stream the matrix and the precondition diagnostics are written to
        /// @return  G[k][l], symmetric, each entry a (nat, 3) force matrix
        std::vector<std::vector<ModuleBase::matrix>> cal_grad_matrix_degenerate(const int ispin,
            const std::vector<int>& group, const std::vector<ModuleBase::matrix>& diag,
            std::ofstream& ofs);
        void test_force();   // test: reproduce the force of ground state
        module_dm::DensityMatrix<T, double> cal_dm_gs();  ///< ground-state density matrix

#ifdef __EXX
        /// Tdata of Exx_LRI is same as T, for the reason, see operator_lr_exx.h
        std::shared_ptr<Exx_LRI<T>> exx_lri = nullptr;
        /// share the ground-state solver's Exx_LRI. It is shared, not stolen: the KS solver
        /// keeps using it on the next ionic step.
        void share_exx_lri(std::shared_ptr<Exx_LRI<double>>&);
        void share_exx_lri(std::shared_ptr<Exx_LRI<std::complex<double>>>&);
        Exx_Info exx_info;
#endif
    };
}

#endif // ABACUS_SOURCE_ESOLVER_ESOLVER_LR_LCAO_TDDFT_H
