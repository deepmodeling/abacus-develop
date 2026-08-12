#pragma once
#include "esolver_lr_lcao_tddft.h"
#include "source_lcao/module_bse/molecular_lri.h"

namespace ModuleESolver
{
    template<typename T> using Real = typename GetTypeReal<T>::type;

    template<typename T, typename TR = double>
    class ESolver_BSE : public ESolver_LR<T, TR> {
    public:
        ESolver_BSE(const Input_para& inp, const std::string& in_dir, const std::string& out_dir)
            : ESolver_LR<T, TR>(inp, in_dir, out_dir), rpa_dir(inp.rpa_outdir)
        {
        }

        ~ESolver_BSE() override {
            //delete this->psi_ks; // already deleted in ESolver_LR
            delete this->psi_ks_global;
        }

        const std::string rpa_dir;
        LR_IO::RI_kRlist kRlist;
        psi::Psi<T>* psi_ks_global; ///< global version of psi_ks
        ModuleBase::matrix eig_gw; ///< GW energy
        double cbm_energy, vbm_energy; //conduction band minimum and valence band maximum
        double direct_gap;
        std::vector<double> tda_ene, full_ene; // in Rydberg

        /// @brief [nspin_types][{nstates, nk* (locc* lvirt}]
        std::vector<ct::Tensor> full_X, full_Y;

        std::unique_ptr<BSE::MolecularLRI<T>> mo_lri;

        virtual void before_all_runners(BaseCell& basecell, const Input_para& inp) override;
        virtual void runner(BaseCell& basecell, int istep) override;
        virtual void after_all_runners(BaseCell& basecell) override;

        void lri_init();

        /// @brief solve IPA without construct and diagonalize matrix
        void ipa_solver();

        /// @brief read in the ground state wave function, gw band energy and occupation
        virtual void read_ks_wfc() override;

        /// @brief init Hartree potential for BSE,
        /// @attention here we don't multiply 2 for singlet or 0 for triplet, we will do it in HamiltBSE
        virtual void init_pot(const Charge& chg_gs) override;

        /// @brief X for tda and also Y for full BSE excitation
        void allocate_eigen_infos();
    };

}
