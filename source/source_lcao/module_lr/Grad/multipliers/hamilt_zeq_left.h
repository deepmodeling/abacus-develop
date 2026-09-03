#pragma once
#include "source_lcao/module_lr/hamilt_casida.h"
#include "hamilt_zeq_ulr.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_lr/Grad/xc/pot_grad_xc.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_hxc.h"
#include "source_lcao/module_lr/Grad/dm_diff/dm_diff.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#ifdef __EXX
#include "source_lcao/module_lr/operator_casida/operator_lr_exx.h"
#endif
namespace LR
{
    template<typename T>
    class Z_vector_L : public HamiltLR<T>
    {
        using ATYPE = typename OperatorLRHxc<T>::MO_TO_AO_TYPE;
#ifdef __EXX
        using ATYPE_EXX = typename OperatorLREXX<T>::MO_TO_AO_TYPE;
#endif
    public:
        Z_vector_L(const std::string& xc_kernel,
            const int& nspin,
            const int& naos,
            const std::vector<int>& nocc,
            const std::vector<int>& nvirt,
            const UnitCell& ucell,
            const std::vector<double>& orb_cutoff,
            const Grid_Driver& gd,
            const psi::Psi<T>& psi_ks,
            const ModuleBase::matrix& eig_ks,
#ifdef __EXX
            std::weak_ptr<Exx_LRI<T>> exx_lri,
            const double& exx_alpha,
#endif 
            std::weak_ptr<PotLRBase> pot_hxc_gs,
            const K_Vectors& kv,
            const std::vector<Parallel_2D>& pX,
            const Parallel_2D& pc,
            const Parallel_Orbitals& pmat,
            const std::string& spin_type)
            : HamiltLR<T>(xc_kernel, nspin, naos, nocc, nvirt, ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot_hxc_gs, kv, pX, pc, pmat, spin_type, PARAM.globalv.global_out_dir)
        {
            ModuleBase::TITLE("Z_vector_L", "Z_vector_L");
            this->DM_trans = LR_Util::make_unique<elecstate::DensityMatrix<T, T>>(&pmat, 1, kv.kvec_d, this->nk);
            LR_Util::initialize_DMR(*this->DM_trans, pmat, ucell, gd, orb_cutoff);
            // Hessian (A+B) with GS XC kernel 
            // 1. diag term in A
            this->ops = new OperatorLRDiag<T>(eig_ks.c, pX[0], kv.get_nks() / nspin, nocc[0], nvirt[0]);
            // 2. $H_{ia}[D^Z]$, equals to $2K_{ab}[D^Z]$ when $D^Z$ is symmetrized.
            // Factor 4 (not 2): the singlet kernel is $K^S_\text{Hxc}=2$`pot_hxc_gs` (not doubled), 
            // while `pot` is the already-doubled singlet potential, so $H^S=2K^S$ here needs 4.
            // The EXX line below is already $2\alpha$ and is consistent.
            hamilt::Operator<T>* op_hz = new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks,
                *this->DM_trans, pot_hxc_gs, ucell, orb_cutoff, gd, kv, pX, pc, pmat,
                { 0 }, 4.0, ATYPE::CC_vo);
            this->ops->add(op_hz);
#ifdef __EXX
            if (gs_is_hybrid())
            {
                hamilt::Operator<T>* op_hz_exx = new OperatorLREXX<T>(nspin, naos, nocc[0], nvirt[0], ucell, psi_ks,
                    *this->DM_trans, exx_lri, kv, pX[0], pc, pmat,
                    2.0 * exx_alpha, //alpha; H=2K when D is symmetrized
                    ATYPE_EXX::CC_vo);
                this->ops->add(op_hz_exx);
            }
#endif
            this->cal_dm_trans = [&, this](const int& is, const T* X)->void
                {
                    const auto psi_ks_is = LR_Util::get_psi_spin(psi_ks, is, this->nk);
#ifdef __MPI
                    std::vector<ct::Tensor>  dm_trans_2d = cal_dm_trans_pblas(X, this->pX[is], psi_ks_is, pc, naos, nocc[is], nvirt[is], pmat);
                    for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat);
#else
                    std::vector<ct::Tensor>  dm_trans_2d = cal_dm_trans_blas(X, psi_ks_is, nocc[is], nvirt[is]);
                    for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
                    // LR_Util::print_tensor<T>(dm_trans_2d[0], "dm_trans_2d[0]", &pmat);
                    // tensor to vector, then set DMK
                    for (int ik = 0;ik < this->nk;++ik) { this->DM_trans->set_DMK_pointer(ik, dm_trans_2d[ik].data<T>()); }
                };
        }
    };

    /// @brief Open-shell (spin-unrestricted) counterpart of `Z_vector_L`.
    ///
    /// Left-hand side of the Z-vector equation: the orbital Hessian built with the
    /// GROUND-STATE kernel,
    ///     $H_{ai\sigma}[D^Z]+Z_{ai\sigma}(\epsilon_{a\sigma}-\epsilon_{i\sigma})$.
    ///
    /// Factor 2.0 on the Hxc blocks (the closed-shell `Z_vector_L` uses 4.0): there
    /// `pot_hxc_gs` is `S2_gs = S2_singlet/2`, so recovering $H^S=2K^S=4\cdot$`pot` needs 4.
    /// Here `pot_hxc_gs` is `S2_updown`, which IS one $K_{\sigma\sigma'}$ component with no
    /// halving, so only the $H=2K$ factor remains. The spin sum $\sum_{\sigma'}$ that the
    /// singlet combination had baked in is now done by the block loop in `ZeqULR::hPsi`.
    template<typename T>
    class Z_vector_UL : public ZeqULR<T>
    {
        using ATYPE = typename OperatorLRHxc<T>::MO_TO_AO_TYPE;
#ifdef __EXX
        using ATYPE_EXX = typename OperatorLREXX<T>::MO_TO_AO_TYPE;
#endif
    public:
        Z_vector_UL(const std::string& xc_kernel,
            const int& nspin,
            const int& naos,
            const std::vector<int>& nocc,
            const std::vector<int>& nvirt,
            const UnitCell& ucell,
            const std::vector<double>& orb_cutoff,
            const Grid_Driver& gd,
            const psi::Psi<T>& psi_ks,
            const ModuleBase::matrix& eig_ks,
#ifdef __EXX
            std::weak_ptr<Exx_LRI<T>> exx_lri,
            const double& exx_alpha,
#endif
            std::weak_ptr<PotLRBase> pot_hxc_gs,
            const K_Vectors& kv,
            const std::vector<Parallel_2D>& pX,
            const Parallel_2D& pc,
            const Parallel_Orbitals& pmat)
            : ZeqULR<T>(nocc, nvirt, pX, kv.get_nks() / nspin),
            naos_(naos), pc_(pc), pmat_(pmat), psi_ks_(psi_ks)
        {
            ModuleBase::TITLE("Z_vector_UL", "Z_vector_UL");
            this->DM_trans = LR_Util::make_unique<elecstate::DensityMatrix<T, T>>(&pmat, 1, kv.kvec_d, this->nk);
            LR_Util::initialize_DMR(*this->DM_trans, pmat, ucell, gd, orb_cutoff);

            // 1. the orbital-energy difference, diagonal blocks only
            this->ops[0] = new OperatorLRDiag<T>(eig_ks.c, pX[0], this->nk, nocc[0], nvirt[0]);
            this->ops[3] = new OperatorLRDiag<T>(eig_ks.c + this->nk * (nocc[0] + nvirt[0]),
                pX[1], this->nk, nocc[1], nvirt[1]);

            // 2. $H_{ai\sigma}[D^Z]=2\sum_{\sigma'}K_{ai\sigma}[D^Z_{\sigma'}]$
            auto newHxc = [&](const int sl, const int sr)
                {
                    return new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks,
                        *this->DM_trans, pot_hxc_gs, ucell, orb_cutoff, gd, kv, pX, pc, pmat,
                        { sl, sr }, T(2.0), ATYPE::CC_vo);
                };
            this->ops[0]->add(newHxc(0, 0));
            this->ops[1] = newHxc(0, 1);
            this->ops[2] = newHxc(1, 0);
            this->ops[3]->add(newHxc(1, 1));

#ifdef __EXX
            // exchange is spin-diagonal ($\delta_{\sigma\sigma'}$), so only blocks 0 and 3.
            // Factor 2*alpha is unchanged from the closed-shell version: the EXX part of the
            // kernel carries no singlet/triplet combination, only $H=2K$.
            if (gs_is_hybrid())
            {
                for (int is : {0, 1})
                {
                    this->psi_ks_spin_.push_back(LR_Util::get_psi_spin(psi_ks, is, this->nk));
                }
                for (int is : {0, 1})
                {
                    this->ops[(is << 1) + is]->add(new OperatorLREXX<T>(nspin, naos, nocc[is], nvirt[is],
                        ucell, this->psi_ks_spin_[is], *this->DM_trans, exx_lri, kv, pX[is], pc, pmat,
                        2.0 * exx_alpha, ATYPE_EXX::CC_vo));
                }
            }
#endif
        }

    protected:
        void set_dm(const int is, const T* const X) const override
        {
            const auto psi_ks_is = LR_Util::get_psi_spin(this->psi_ks_, is, this->nk);
#ifdef __MPI
            this->dm_buf_ = cal_dm_trans_pblas(X, this->pX[is], psi_ks_is, this->pc_,
                this->naos_, this->nocc[is], this->nvirt[is], this->pmat_);
            for (auto& t : this->dm_buf_) { LR_Util::matsym(t.template data<T>(), this->naos_, this->pmat_); }
#else
            this->dm_buf_ = cal_dm_trans_blas(X, psi_ks_is, this->nocc[is], this->nvirt[is]);
            for (auto& t : this->dm_buf_) { LR_Util::matsym(t.template data<T>(), this->naos_); }
#endif
            for (int ik = 0;ik < this->nk;++ik)
            {
                this->DM_trans->set_DMK_pointer(ik, this->dm_buf_[ik].template data<T>());
            }
        }

    private:
        const int naos_ = 1;
        const Parallel_2D& pc_;
        const Parallel_Orbitals& pmat_;
        const psi::Psi<T>& psi_ks_;
        std::vector<psi::Psi<T>> psi_ks_spin_;
        std::unique_ptr<elecstate::DensityMatrix<T, T>> DM_trans;
        /// the tensors `DM_trans` points into; kept alive for the whole `act` chain
        mutable std::vector<ct::Tensor> dm_buf_;
    };
}
