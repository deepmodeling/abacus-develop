#pragma once
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
#include "source_hamilt/module_gint/gint_info.h"
#include "source_estate/module_pot/potential_new.h"
#ifdef __EXX
// #include <RI/physics/Exx.h>
#include "source_lcao/module_ri/exx_lri.h"
#include "source_hamilt/module_xc/exx_info.h" // for Exx_Info value member
#endif
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
        int nocc_max = 1;   ///< nelec/2
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
        /// index of the relaxed state inside `pelec->ekb`
        int target_ekb_offset_() const
        { return this->openshell ? this->inp_->lr_target_state
                                 : this->target_is_ * this->nstates + this->inp_->lr_target_state; }

        std::vector<std::string> spin_types;

        std::unique_ptr<ModuleGint::GintInfo> gint_info_ = nullptr;
        void set_gint();

        /// @brief variables for parallel distribution of KS orbitals
        Parallel_2D paraC_;
        /// @brief variables for parallel distribution of excited states
        std::vector<Parallel_2D> paraX_;
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
        /// Solve the Z-vector equation. `istate_only >= 0` restricts it to that one excited state
        /// (the returned tensor then holds a single block); -1 solves all `nstates`.
        ct::Tensor solve_zvector_eqation(const int ispin, const int istate_only = -1);
        /// Excited-state gradients d(Omega)/dR, one matrix per state solved. `istate_only` as above:
        /// geometry relaxation follows a single state, and the Z-vector solve dominates the cost.
        std::vector<ModuleBase::matrix> cal_force(const int ispin, const int istate_only = -1);
        /// open-shell (spin-unrestricted) excited-state force: X holds [up | down] and every
        /// density matrix has two independent channels
        std::vector<ModuleBase::matrix> cal_force_openshell(const int istate_only = -1);
        void test_force();   // test: reproduce the force of ground state
        elecstate::DensityMatrix<T, double> cal_dm_gs();  ///< ground-state density matrix

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
