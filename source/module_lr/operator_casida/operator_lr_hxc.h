#pragma once
#include "module_cell/klist.h"
#include "module_hamilt_general/operator.h"
#include "module_lr/utils/gint_template.h"
#include "module_hamilt_lcao/module_gint/grid_technique.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_lr/potentials/pot_hxc_lrtd.h"
#include "module_lr/utils/lr_util.h"
#include "module_lr/utils/lr_util_hcontainer.h"
namespace LR
{
    /// @brief  Hxc part of A operator for LR-TDDFT
    template<typename T = double, typename Device = base_device::DEVICE_CPU>
    class OperatorLRHxc : public hamilt::Operator<T, Device>
    {
    public:
        //when nspin=2, nks is 2 times of real number of k-points. else (nspin=1 or 4), nks is the real number of k-points
        OperatorLRHxc(const int& nspin,
            const int& naos,
            const int& nocc,
            const int& nvirt,
            const psi::Psi<T, Device>* psi_ks_in,
            std::vector<std::unique_ptr<elecstate::DensityMatrix<T, T>>>& DM_trans_in,
            typename TGint<T>::type* gint_in,
            std::weak_ptr<PotHxcLR> pot_in,
            const UnitCell& ucell_in,
            const std::vector<double>& orb_cutoff,
            Grid_Driver& gd_in,
            const K_Vectors& kv_in,
            Parallel_2D* pX_in,
            Parallel_2D* pc_in,
            Parallel_Orbitals* pmat_in)
            : nspin(nspin), naos(naos), nocc(nocc), nvirt(nvirt),
            psi_ks(psi_ks_in), DM_trans(DM_trans_in), gint(gint_in), pot(pot_in),
            ucell(ucell_in), orb_cutoff_(orb_cutoff), gd(gd_in), kv(kv_in),
            pX(pX_in), pc(pc_in), pmat(pmat_in)
        {
            ModuleBase::TITLE("OperatorLRHxc", "OperatorLRHxc");
            this->cal_type = hamilt::calculation_type::lcao_gint;
            this->act_type = 2;
            this->is_first_node = true;
            this->hR = std::unique_ptr<hamilt::HContainer<T>>(new hamilt::HContainer<T>(pmat_in));
            this->initialize_HR(*this->hR, ucell_in, gd_in, pmat_in);
            this->DM_trans[0]->init_DMR(*this->hR);
        };
        ~OperatorLRHxc() { };

        void init(const int ik_in) override {};

        // virtual psi::Psi<T> act(const psi::Psi<T>& psi_in) const override;
        virtual void act(const psi::Psi<T>& psi_in, psi::Psi<T>& psi_out, const int nbands) const override;
    private:
        template<typename TR>   //T=double, TR=double; T=std::complex<double>, TR=std::complex<double>/double
        void initialize_HR(hamilt::HContainer<TR>& hR, const UnitCell& ucell, Grid_Driver& gd, const Parallel_Orbitals* pmat) const
        {
            for (int iat1 = 0; iat1 < ucell.nat; iat1++)
            {
                auto tau1 = ucell.get_tau(iat1);
                int T1, I1;
                ucell.iat2iait(iat1, &I1, &T1);
                AdjacentAtomInfo adjs;
                gd.Find_atom(ucell, tau1, T1, I1, &adjs);
                for (int ad = 0; ad < adjs.adj_num + 1; ++ad)
                {
                    const int T2 = adjs.ntype[ad];
                    const int I2 = adjs.natom[ad];
                    int iat2 = this->ucell.itia2iat(T2, I2);
                    if (pmat->get_row_size(iat1) <= 0 || pmat->get_col_size(iat2) <= 0) { continue; }
                    const ModuleBase::Vector3<int>& R_index = adjs.box[ad];
                    if (ucell.cal_dtau(iat1, iat2, R_index).norm() * this->ucell.lat0 >= orb_cutoff_[T1] + orb_cutoff_[T2]) { continue; }
                    hamilt::AtomPair<TR> tmp(iat1, iat2, R_index.x, R_index.y, R_index.z, pmat);
                    hR.insert_pair(tmp);
                }
            }
            hR.allocate(nullptr, true);
            hR.set_paraV(pmat);
            if (std::is_same<T, double>::value) { hR.fix_gamma(); }
        }
        template<typename TR>
        void init_DM_trans(const int& nbands, std::vector<std::unique_ptr<elecstate::DensityMatrix<T, TR>>>& DM_trans)const
        {
            // LR_Util::print_DMR(*this->DM_trans[0], ucell.nat, "DMR[ib=" + std::to_string(0) + "]");
            if (this->next_op != nullptr)
            {
                int prev_size = DM_trans.size();
                if (prev_size > nbands) { for (int ib = nbands;ib < prev_size;++ib) { DM_trans[ib].reset(); } }
                DM_trans.resize(nbands);
                for (int ib = prev_size;ib < nbands;++ib)
                {
                    // the first dimenstion of DensityMatrix is nk=nks/nspin 
                    DM_trans[ib] = LR_Util::make_unique<elecstate::DensityMatrix<T, TR>>(this->pmat, this->nspin, this->kv.kvec_d, this->kv.get_nks() / nspin);
                    DM_trans[ib]->init_DMR(*this->hR);
                }
            }
        }
        void grid_calculation(const int& nbands, const int& iband_dm)const;

        //global sizes
        const int& nspin;
        const int nspin_solve = 1;    ///< in singlet-triplet calculation, the Casida equation is solved respectively so nspin_solve in a single problem is 1
        const int& naos;
        const int& nocc;
        const int& nvirt;
        const K_Vectors& kv;
        /// ground state wavefunction
        const psi::Psi<T, Device>* psi_ks = nullptr;

        /// transition density matrix
        std::vector<std::unique_ptr<elecstate::DensityMatrix<T, T>>>& DM_trans;

        /// transition hamiltonian in AO representation
        std::unique_ptr<hamilt::HContainer<T>> hR = nullptr;

        /// parallel info
        Parallel_2D* pc = nullptr;
        Parallel_2D* pX = nullptr;
        Parallel_Orbitals* pmat = nullptr;

        std::weak_ptr<PotHxcLR> pot;

        typename TGint<T>::type* gint = nullptr;

        const UnitCell& ucell;
        std::vector<double> orb_cutoff_;
        Grid_Driver& gd;

        bool tdm_sym = false; ///< whether transition density matrix is symmetric

        /// test
        mutable bool first_print = true;
    };
}
