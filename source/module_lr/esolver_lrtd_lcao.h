#pragma once
#include "module_esolver/esolver_fp.h"
#include "module_parameter/input_parameter.h"
#include "module_cell/unitcell.h"
#include "module_hamilt_general/hamilt.h"
#include "module_hsolver/hsolver.h"
#include "module_elecstate/elecstate_lcao.h"

#include <vector>   //future tensor
#include <memory>

#include "module_esolver/esolver_ks_lcao.h" //for the move constructor
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_lr/potentials/pot_hxc_lrtd.h"
#include "module_lr/hamilt_casida.h"
#ifdef __EXX
// #include <RI/physics/Exx.h>
#include "module_ri/Exx_LRI.h"
#endif
namespace LR
{
    ///Excited State Solver: Linear Response TDDFT (Tamm Dancoff Approximation) 
    template<typename T, typename TR = double>
    class ESolver_LR : public ModuleESolver::ESolver_FP
    {
    public:
        /// @brief  a move constructor from ESolver_KS_LCAO
        ESolver_LR(ModuleESolver::ESolver_KS_LCAO<T, TR>&& ks_sol, const Input_para& inp, UnitCell& ucell);
        /// @brief a from-scratch constructor
        ESolver_LR(const Input_para& inp, UnitCell& ucell);
        ~ESolver_LR() {
            delete this->psi_ks;
        }

        ///input: input, call, basis(LCAO), psi(ground state), elecstate
        // initialize sth. independent of the ground state
        virtual void before_all_runners(const Input_para& inp, UnitCell& cell) override {};
        virtual void runner(int istep, UnitCell& ucell) override;
        virtual void after_all_runners() override;

        virtual double cal_energy()  override { return 0.0; };
        virtual void cal_force(ModuleBase::matrix& force) override {};
        virtual void cal_stress(ModuleBase::matrix& stress) override {};

    protected:
        const Input_para& input;
        const UnitCell& ucell;
        std::vector<double> orb_cutoff_;

        // not to use ElecState because 2-particle state is quite different from 1-particle state.
        // implement a independent one (ExcitedState) to pack physical properties if needed.
        // put the components of ElecState here: 
        std::vector<std::shared_ptr<PotHxcLR>> pot;

        // ground state info 

        /// @brief ground state wave function
        psi::Psi<T>* psi_ks = nullptr;

        /// @brief ground state bands, read from the file, or moved from ESolver_FP::pelec.ekb
        ModuleBase::matrix eig_ks;///< energy of ground state

        /// @brief Excited state wavefunction (locc, lvirt are local size of nocc and nvirt in each process)
        /// size of X: [neq][{nstate, nloc_per_band}], namely:
        /// - [nspin][{nstates, nk* (locc* lvirt}] for close- shell,
        /// -  [1][{nstates, nk * (locc[0] * lvirt[0]) + nk * (locc[1] * lvirt[1])}] for open-shell
        std::vector<ct::Tensor> X;
        int nloc_per_band = 1;

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

        Grid_Technique gt_;
        Gint_Gamma gint_g_;
        Gint_k gint_k_;
        typename TGint<T>::type* gint_ = nullptr;
        void set_gint();

        /// @brief variables for parallel distribution of KS orbitals
        Parallel_2D paraC_;
        /// @brief variables for parallel distribution of excited states
        std::vector<Parallel_2D> paraX_;
        /// @brief variables for parallel distribution of matrix in AO representation
        Parallel_Orbitals paraMat_;

        TwoCenterBundle two_center_bundle_;

        /// @brief allocate and set the inital value of X
        void setup_eigenvectors_X();
        void set_X_initial_guess();

        /// @brief read in the ground state wave function, band energy and occupation
        void read_ks_wfc();
        /// @brief  read in the ground state charge density
        void read_ks_chg(Charge& chg);

        void init_pot(const Charge& chg_gs);

        /// @brief check the legality of the input parameters
        void parameter_check() const;
        /// @brief set nocc, nvirt, nbasis, npairs and nstates
        void set_dimension();
        /// reset nocc, nvirt, npairs after read ground-state wavefunction when nspin=2
        void reset_dim_spin2();

#ifdef __EXX
        /// Tdata of Exx_LRI is same as T, for the reason, see operator_lr_exx.h
        std::shared_ptr<Exx_LRI<T>> exx_lri = nullptr;
        void move_exx_lri(std::shared_ptr<Exx_LRI<double>>&);
        void move_exx_lri(std::shared_ptr<Exx_LRI<std::complex<double>>>&);
        Exx_Info& exx_info;
#endif
    };
}
