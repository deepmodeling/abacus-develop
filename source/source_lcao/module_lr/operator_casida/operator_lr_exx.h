#ifdef __EXX
#ifndef ABACUS_SOURCE_LCAO_MODULE_LR_OPERATOR_CASIDA_OPERATOR_LR_EXX_H
#define ABACUS_SOURCE_LCAO_MODULE_LR_OPERATOR_CASIDA_OPERATOR_LR_EXX_H

#include "source_hamilt/operator.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_ri/exx_lri.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/dm_trans/dm_diff.h"
namespace LR
{

    /// @brief  Exx part of A operator
    /// The single source of truth is `LR_Util::hybrid_xc_list` -- see there for what the list
    /// does and does not decide.
    inline const std::set<std::string>& exx_kernel_list() { return LR_Util::hybrid_xc_list(); };

    /// @brief does the GROUND STATE carry exact exchange, i.e. is `dft_functional` a hybrid?
    inline bool gs_is_hybrid(const std::string& dft_functional) { return exx_kernel_list().count(LR_Util::tolower(dft_functional)) > 0; };

    template<typename T = double>
    class OperatorLREXX : public hamilt::Operator<T, base_device::DEVICE_CPU>
    {
    public:
        /// @brief type of molecular orbital to atomic orbital transformation:
        /// CC_vo: MO = C_v^* AO  C_o;
        /// CC_oo: MO = C_o^* AO  C_o;
        /// CXC: MO = C_v^* AO X C_v- C_o^* X^* AO C_o
        /// CXC_o: MO = C_o^* X^* AO C_o
        enum class MO_TO_AO_TYPE { CC_vo, CC_oo, CXC, CXC_o };

        OperatorLREXX(const int& nspin,
            const int& naos,
            const int& nocc,
            const int& nvirt,
            const UnitCell& ucell_in,
            const psi::Psi<T>& psi_ks_in,
            const module_dm::DensityMatrix<T, T>& DM_trans_in,
            // HContainer<double>* hR_in,
            std::weak_ptr<Exx_LRI<T>> exx_lri_in,
            const K_Vectors& kv_in,
            const Parallel_2D& pX_in,
            const Parallel_2D& pc_in,
            const Parallel_Orbitals& pmat_in,
            const bool cal_force,
            const double& alpha = 1.0,
            const MO_TO_AO_TYPE dm_pq_in = MO_TO_AO_TYPE::CC_vo,
            const std::vector<int>& aims_nbasis = {},
            const hamilt::calculation_type cal_type_in = hamilt::calculation_type::lr_dmtrans_exx)
            : nspin(nspin), naos(naos), nocc(nocc), nvirt(nvirt), nk(kv_in.get_nks() / nspin),
            psi_ks(psi_ks_in), DM_trans(DM_trans_in), exx_lri(exx_lri_in), kv(kv_in),
            pX(pX_in), pc(pc_in), pmat(pmat_in), ucell(ucell_in), alpha(alpha), dm_pq_(dm_pq_in),
            aims_nbasis(aims_nbasis)
        {
            ModuleBase::TITLE("OperatorLREXX", "OperatorLREXX");
            this->cal_type = cal_type_in;
            this->is_first_node = false;

            // reduce psi_ks for later use
            this->psi_ks_full.resize(this->nk, nocc + nvirt, this->naos);
            this->psi_ks_full.zero_out();
            for (int ik = 0;ik < nk;++ik)
            {
                LR_Util::gather_2d_to_full(this->pc, &this->psi_ks(ik, 0, 0), &this->psi_ks_full(ik, 0, 0), false, this->naos, nocc + nvirt);
            }
            if (!this->exx_lri.expired())
            {
                this->exx_lri.lock()->Hexxs.resize(1);
            }
        };

        void init(const int ik_in) override {};

        virtual void act(const int nbands,
                         const int nbasis,
                         const int npol,
                         const T* psi_in,
                         T* hpsi,
                         const int ngk_ik = 0,
                         const bool is_first_node = false) const override;

    private:
        //global sizes
        const int nspin = 1;
        const int naos = 1;
        const int nocc = 1;
        const int nvirt = 1;
        const int nk = 1;  ///< number of k-points
        const double alpha = 1.0;   //(allow non-ref constant)
        MO_TO_AO_TYPE dm_pq_ = MO_TO_AO_TYPE::CC_vo;
        const bool cal_dm_trans = false;
        const bool tdm_sym = false; ///< whether transition density matrix is symmetric
        const K_Vectors& kv;
        /// ground state wavefunction
        const psi::Psi<T>& psi_ks = nullptr;
        psi::Psi<T> psi_ks_full;

        /// transition density matrix 
        const module_dm::DensityMatrix<T, T>& DM_trans;

        /// C, V tensors of RI, and LibRI interfaces
        /// gamma_only: T=double, Tpara of exx (equal to Tpara of Ds(R) ) is also double 
        /// multi-k: both the transition density and the exchange response are complex
        /// so TR in DensityMatrix and Tdata in Exx_LRI are all equal to T
        std::weak_ptr<Exx_LRI<T>> exx_lri;

        const UnitCell& ucell;

        ///parallel info
        const Parallel_2D& pc;
        const Parallel_2D& pX;
        const Parallel_Orbitals& pmat;

        /// number of basis functions per type in the aims benchmark (empty for the normal case)
        const std::vector<int> aims_nbasis;

        /// only for gradient calculation


        /// Build and communicate the exchange response once per input density.
        void cal_Hs() const;

        /// Gamma and complex Bloch matrix projection, including benchmark indices.
        void project_k(const T* psi_in, T* hpsi) const;

        void cal_coxt_cvx(const T* x_istate, psi::Psi<T>& coxt_full, psi::Psi<T>& cvx_full) const   // C_o X^T, C_v X (only for gradients)
        {
            ModuleBase::TITLE("OperatorLREXX", "cal_coxt_cvx");
            const auto& c = this->psi_ks;
            // allocate local cvx
            Parallel_2D pcx;
            LR_Util::setup_2d_division(pcx, pX.get_block_size(), naos, nocc, pX.blacs_ctxt);
            ct::Tensor cvx(ct::DataTypeToEnum<T>::value, DEV::CpuDevice, { pcx.get_col_size(), pcx.get_row_size() });

            // allocate local coxt
            Parallel_2D pcxt;
            LR_Util::setup_2d_division(pcxt, pX.get_block_size(), naos, nvirt, pX.blacs_ctxt);
            ct::Tensor coxt(ct::DataTypeToEnum<T>::value, DEV::CpuDevice, { pcxt.get_col_size(), pcxt.get_row_size() });

            // calculate global coxt_full, cvx_full
            cvx_full.zero_out();
            coxt_full.zero_out();
            for (int ik = 0;ik < nk;++ik)
            {
                c.fix_k(ik);
                const int start = ik * pX.get_local_size();
                CvX(c.get_pointer(), pc, x_istate + start, pX, naos, nocc, nvirt, cvx.data<T>(), pcx);
                LR_Util::gather_2d_to_full(pcx, cvx.data<T>(), &cvx_full(ik, 0, 0), false, naos, nocc);
                CoXT(c.get_pointer(), pc, x_istate + start, pX, naos, nocc, nvirt, coxt.data<T>(), pcxt);
                LR_Util::gather_2d_to_full(pcxt, coxt.data<T>(), &coxt_full(ik, 0, 0), false, naos, nvirt);
            }
        }

    };
}
#endif

#endif // ABACUS_SOURCE_LCAO_MODULE_LR_OPERATOR_CASIDA_OPERATOR_LR_EXX_H
