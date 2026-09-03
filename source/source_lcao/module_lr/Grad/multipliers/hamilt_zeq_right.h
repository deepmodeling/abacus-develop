#pragma once
#include "source_hamilt/hamilt.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_lr/Grad/xc/pot_grad_xc.h"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_hxc.h"
#include "hamilt_zeq_ulr.h"
#include "source_lcao/module_lr/Grad/xc/operator_gxc_ulr.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#ifdef __EXX
#include "source_lcao/module_lr/operator_casida/operator_lr_exx.h"
#endif
namespace LR
{
    template<typename T>
    class Z_vector_R : public HamiltLR<T>
    {
        using ATYPE = typename OperatorLRHxc<T>::MO_TO_AO_TYPE;
#ifdef __EXX
        using ATYPE_EXX = typename OperatorLREXX<T>::MO_TO_AO_TYPE;
#endif
    public:
        Z_vector_R(const std::string& xc_kernel,
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
            std::weak_ptr<PotHxcLR> pot,
            std::weak_ptr<PotHxcLR> pot_hxc_gs,
            const K_Vectors& kv,
            const std::vector<Parallel_2D>& pX,
            const Parallel_2D& pc,
            const Parallel_Orbitals& pmat,
            const std::string& spin_type = "singlet")
            : HamiltLR<T>(xc_kernel, nspin, naos, nocc, nvirt, ucell, orb_cutoff, gd, psi_ks, eig_ks,
#ifdef __EXX
                exx_lri, exx_alpha,
#endif
                pot, kv, pX, pc, pmat, spin_type, PARAM.globalv.global_out_dir)
        {
            ModuleBase::TITLE("Z_vector_R", "Z_vector_R");

            this->DM_trans = LR_Util::make_unique<elecstate::DensityMatrix<T, T>>(&pmat, 1, kv.kvec_d, this->nk);
            LR_Util::initialize_DMR(*this->DM_trans, pmat, ucell, gd, orb_cutoff);
            this->DM_diff = LR_Util::make_unique<elecstate::DensityMatrix<T, T>>(&pmat, 1, kv.kvec_d, this->nk);
            LR_Util::initialize_DMR(*this->DM_diff, pmat, ucell, gd, orb_cutoff);

            // note: calculation_type cannot repeated, or it will be ignored in ops->add()
            // 1. $2\sum_bX_{ib}K_{ab}[D^X]-2\sum_jX_{ja}K_{ij}[D^X]$
            // kernel: excited state
            this->ops = new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks,
                *this->DM_trans, pot, ucell, orb_cutoff, gd, kv, pX, pc, pmat,
                { 0 }, T(-2.0), ATYPE::CXC);
#ifdef __EXX
            if (exx_kernel_list().count(xc_kernel))
            {
                hamilt::Operator<T>* op_hz_exx = new OperatorLREXX<T>(nspin, naos, nocc[0], nvirt[0], ucell, psi_ks,
                    *this->DM_trans, exx_lri, kv, pX[0], pc, pmat,
                    -2.0 * exx_alpha, //alpha; H=2K when D is symmetrized
                    ATYPE_EXX::CXC, {}, hamilt::calculation_type::lr_dmtrans_exx);
                this->ops->add(op_hz_exx);
            }
#endif
            // 2. $H_{ia}[T]$, equals to $2K_{ab}[T]$ when $T$ is symmetrized
            // kernel: ground state
            // Factor -4 (not -2), for the same reason as in `hamilt_zeq_left.h`:
            // $K^S_\text{Hxc}=2$`pot_hxc_gs`, so $H^S=2K^S$ needs 4. 
            hamilt::Operator<T>* op_ht = new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks,
                *this->DM_diff, pot_hxc_gs, ucell, orb_cutoff, gd, kv, pX, pc, pmat,
                { 0 }, T(-4.0), ATYPE::CC_vo, hamilt::calculation_type::lr_dmdiff_hxc);
            this->ops->add(op_ht);
#ifdef __EXX
            if (gs_is_hybrid())
            {
                hamilt::Operator<T>* op_ht_exx = new OperatorLREXX<T>(nspin, naos, nocc[0], nvirt[0], ucell, psi_ks,
                    *this->DM_diff, exx_lri, kv, pX[0], pc, pmat,
                    -2.0 * exx_alpha, //alpha; H=2K when D is symmetrized
                    ATYPE_EXX::CC_vo, {}, hamilt::calculation_type::lr_dmdiff_exx);
                this->ops->add(op_ht_exx);
            }
#endif

            // 3. $2\sum_{jb,kc} g^{xc}_{ia, jb, kc}X_{jb}X_{kc}$
            // NOT singlet-only for a local functional where the triplet kernel is
            // $K^T_{xc}=f_{uu}-f_{ud}\ne0$ -- exactly what `PotHxcLR`'s `S2_triplet` branch
            // evaluates -- so $\partial K^T$ carries a $g^{xc}$ term too, with the "-" spin
            // combination. `PotGradXCLR` picks it via the `triplet` flag.
            // (LR-Grad-formulas/GGA-kxc-to-v积分公式.md section 5.)
            if (LR_Util::has_local_xc(xc_kernel))
            {
                this->pot_grad = std::make_shared<PotGradXCLR>(pot.lock()->xc_kernel_components(), pot.lock()->get_rho_basis(), ucell, pot.lock()->nrxx, spin_type == "triplet");
                hamilt::Operator<T>* op_gxc = new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks,
                    *this->DM_trans, this->pot_grad, ucell, orb_cutoff, gd, kv, pX, pc, pmat,
                    { 0 }, T(-2.0), ATYPE::CC_vo, hamilt::calculation_type::lr_dmtrans_gxc);
                assert(op_gxc != nullptr);
                this->ops->add(op_gxc);
            }
            // // test: op_ht only 
            // delete this->ops;
            // this->ops = new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks,
            //     *this->DM_diff, pot_hxc_gs, ucell, orb_cutoff, gd, kv, pX, pc, pmat,
            //     { 0 }, T(-2.0), ATYPE::CC_vo);

            // $D^X$ is fed to the CXC operators TRANSPOSED.
            //
            // The RHS needs $K^{S/T}_{ba}[D^X]$ -- the same kernel the Casida equation was solved
            // with, which `HamiltLR` builds from the un-symmetrized $D^X$ (`tdm_sym = false`).
            // `CVCX_virt`/`CVCX_occ` give the kernel matrix with its two MO indices in
            // the opposite order to what this term needs, and since
            //     $(K[D])^T = K[D^T]$,
            // transposing the density matrix on the way in restores it.
            this->cal_dm_trans = [&, this](const int& is, const T* X)->void
                {
                    const auto psi_ks_is = LR_Util::get_psi_spin(psi_ks, is, this->nk);
#ifdef __MPI
                    std::vector<ct::Tensor> dm_trans_2d = cal_dm_trans_pblas(X, this->pX[is], psi_ks_is, pc, naos, nocc[is], nvirt[is], pmat);
                    for (auto& t : dm_trans_2d) LR_Util::mattrans(t.data<T>(), naos, pmat);
#else
                    std::vector<ct::Tensor> dm_trans_2d = cal_dm_trans_blas(X, psi_ks_is, nocc[is], nvirt[is]);
                    for (auto& t : dm_trans_2d)
                    {
                        T* d = t.data<T>();
                        for (int u = 0;u < naos;++u) { for (int v = u + 1;v < naos;++v) { std::swap(d[u * naos + v], d[v * naos + u]); } }
                    }
#endif
                    for (int ik = 0;ik < this->nk;++ik) { this->DM_trans->set_DMK_pointer(ik, dm_trans_2d[ik].data<T>()); }
                };

            this->cal_dm_diff = [&, this](const int& is, const T* const X)->void
                {
                    const auto psi_ks_is = LR_Util::get_psi_spin(psi_ks, is, this->nk);
#ifdef __MPI
                    std::vector<ct::Tensor> dm_diff_2d = cal_dm_diff_pblas(X, this->pX[is], psi_ks_is, pc, naos, nocc[is], nvirt[is], pmat);
                    for (auto& t : dm_diff_2d) LR_Util::matsym(t.data<T>(), naos, pmat);
#else
                    std::vector<ct::Tensor> dm_diff_2d = cal_dm_diff_blas(X, psi_ks_is, naos, nocc[is], nvirt[is]);
                    for (auto& t : dm_diff_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
                    for (int ik = 0;ik < this->nk;++ik) { this->DM_diff->set_DMK_pointer(ik, dm_diff_2d[ik].data<T>()); }
                    // std::cout << "difference density matrix" << std::endl;
                    // for (int ik = 0;ik < this->nk;++ik) { LR_Util::print_value(dm_diff_2d[ik].data<T>(), naos, naos); }
                    // std::cout << "test: set dm_diff to zero" << std::endl;
                    // for (int ik = 0;ik < this->nk;++ik) { dm_diff_2d[ik].zero(); }
                };
        }
        virtual void hPsi(const T* const psi, T* const hpsi, const int ld_psi, const int nband) const override
        {
            assert(ld_psi == this->nk * this->pX[0].get_local_size());
            for (int ib = 0;ib < nband;++ib)
            {
                const int offset = ib * ld_psi;
                this->cal_dm_trans(0, psi + offset);  // transition density matrix, only for test
                this->cal_dm_diff(0, psi + offset);  // difference density matrix
                hamilt::Operator<T>* node(this->ops);
                while (node != nullptr)
                {
                    node->act(/*nband=*/1, ld_psi, /*npol=*/1, psi + offset, hpsi + offset);
                    node = (hamilt::Operator<T>*)(node->next_op);
                }
            }
        }

    private:
        std::unique_ptr<elecstate::DensityMatrix<T, T>> DM_diff;
        std::function<void(const int&, const T* const)> cal_dm_diff;
        std::shared_ptr<PotLRBase> pot_grad;
    };

    /// @brief Open-shell (spin-unrestricted) counterpart of `Z_vector_R`.
    ///
    /// Right-hand side of the Z-vector equation (the operators produce $-R$):
    /// $$R_{ai\sigma}=H_{ai\sigma}[T]+\sum_bX_{bi\sigma}H_{ba\sigma}[D^X]
    ///   -\sum_kX_{ak\sigma}H_{ik\sigma}[D^X]
    ///   +2\sum_{\kappa\lambda\sigma'}\sum_{\alpha\beta\sigma''}D^X D^X g^{xc}.$$
    ///
    /// Factors relative to the closed-shell `Z_vector_R`:
    /// - the $K[D^X]$ term keeps -2.0. There `pot` is `S2_singlet`, which already IS $K^S$,
    ///   so the factor is just $H=2K$; here `pot` is `S2_updown`, which is one
    ///   $K_{\sigma\sigma'}$ component, and the $\sum_{\sigma'}$ is done by the block loop.
    /// - the $H[T]$ term goes -4.0 -> -2.0, because `pot_hxc_gs` is no longer the halved
    ///   `S2_gs` but `S2_updown` (same argument as in `Z_vector_UL`).
    template<typename T>
    class Z_vector_UR : public ZeqULR<T>
    {
        using ATYPE = typename OperatorLRHxc<T>::MO_TO_AO_TYPE;
#ifdef __EXX
        using ATYPE_EXX = typename OperatorLREXX<T>::MO_TO_AO_TYPE;
#endif
    public:
        Z_vector_UR(const std::string& xc_kernel,
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
            std::weak_ptr<PotHxcLR> pot,
            std::weak_ptr<PotHxcLR> pot_hxc_gs,
            const K_Vectors& kv,
            const std::vector<Parallel_2D>& pX,
            const Parallel_2D& pc,
            const Parallel_Orbitals& pmat)
            : ZeqULR<T>(nocc, nvirt, pX, kv.get_nks() / nspin),
            naos_(naos), pc_(pc), pmat_(pmat), psi_ks_(psi_ks)
        {
            ModuleBase::TITLE("Z_vector_UR", "Z_vector_UR");
            this->DM_trans = LR_Util::make_unique<elecstate::DensityMatrix<T, T>>(&pmat, 1, kv.kvec_d, this->nk);
            LR_Util::initialize_DMR(*this->DM_trans, pmat, ucell, gd, orb_cutoff);
            this->DM_diff = LR_Util::make_unique<elecstate::DensityMatrix<T, T>>(&pmat, 1, kv.kvec_d, this->nk);
            LR_Util::initialize_DMR(*this->DM_diff, pmat, ucell, gd, orb_cutoff);

            // 1. $\sum_bX_{bi\sigma}H_{ba\sigma}[D^X]-\sum_kX_{ak\sigma}H_{ik\sigma}[D^X]$,
            //    excited-state kernel, $D^X$ fed TRANSPOSED (see `set_dm` below)
            auto newCXC = [&](const int sl, const int sr)
                {
                    return new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks,
                        *this->DM_trans, pot, ucell, orb_cutoff, gd, kv, pX, pc, pmat,
                        { sl, sr }, T(-2.0), ATYPE::CXC);
                };
            this->ops[0] = newCXC(0, 0);
            this->ops[1] = newCXC(0, 1);
            this->ops[2] = newCXC(1, 0);
            this->ops[3] = newCXC(1, 1);

            // 2. $H_{ai\sigma}[T]$, ground-state kernel
            auto newHT = [&](const int sl, const int sr)
                {
                    return new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks,
                        *this->DM_diff, pot_hxc_gs, ucell, orb_cutoff, gd, kv, pX, pc, pmat,
                        { sl, sr }, T(-2.0), ATYPE::CC_vo, hamilt::calculation_type::lr_dmdiff_hxc);
                };
            this->ops[0]->add(newHT(0, 0));
            this->ops[1]->add(newHT(0, 1));
            this->ops[2]->add(newHT(1, 0));
            this->ops[3]->add(newHT(1, 1));

            // 3. $2\sum_{\sigma'\sigma''}D^X_{\sigma'}D^X_{\sigma''}g^{xc}_{\dots,ai\tau}$.
            // Quadratic in $D^X$, so it lives outside the 2x2 block structure (see `band_extra`).
            // Factor -2.0, straight from the spin-orbital formula -- the closed-shell code uses the
            // same number only because its `PotGradXCLR` carries the doubled S/T combination.
            if (LR_Util::has_local_xc(xc_kernel))
            {
                this->gxc_ = LR_Util::make_unique<OperatorGxcULR<T>>(pot.lock()->xc_kernel_components(),
                    pot.lock()->get_rho_basis(), ucell, orb_cutoff, gd, kv, pmat, pc, psi_ks,
                    nocc, nvirt, naos, pX, pX, LR_Util::MO_TYPE::VO, T(-2.0));
            }

#ifdef __EXX
            for (int is : {0, 1}) { this->psi_ks_spin_.push_back(LR_Util::get_psi_spin(psi_ks, is, this->nk)); }
            if (exx_kernel_list().count(xc_kernel))
            {
                for (int is : {0, 1})
                {
                    this->ops[(is << 1) + is]->add(new OperatorLREXX<T>(nspin, naos, nocc[is], nvirt[is],
                        ucell, this->psi_ks_spin_[is], *this->DM_trans, exx_lri, kv, pX[is], pc, pmat,
                        -2.0 * exx_alpha, ATYPE_EXX::CXC, {}, hamilt::calculation_type::lr_dmtrans_exx));
                }
            }
            if (gs_is_hybrid())
            {
                for (int is : {0, 1})
                {
                    this->ops[(is << 1) + is]->add(new OperatorLREXX<T>(nspin, naos, nocc[is], nvirt[is],
                        ucell, this->psi_ks_spin_[is], *this->DM_diff, exx_lri, kv, pX[is], pc, pmat,
                        -2.0 * exx_alpha, ATYPE_EXX::CC_vo, {}, hamilt::calculation_type::lr_dmdiff_exx));
                }
            }
#endif
        }

    protected:
        void band_extra(const T* const X, T* const out) const override
        {
            if (this->gxc_) { this->gxc_->act(X, out); }
        }

        /// Rebuild BOTH density matrices from the `is` block of X: the CXC operators read
        /// $D^X$ (transposed, un-symmetrized -- see the closed-shell `Z_vector_R` for why),
        /// the CC_vo operators read the difference density matrix $T$ (symmetrized).
        void set_dm(const int is, const T* const X) const override
        {
            const auto psi_ks_is = LR_Util::get_psi_spin(this->psi_ks_, is, this->nk);
#ifdef __MPI
            this->dmx_buf_ = cal_dm_trans_pblas(X, this->pX[is], psi_ks_is, this->pc_,
                this->naos_, this->nocc[is], this->nvirt[is], this->pmat_);
            for (auto& t : this->dmx_buf_) { LR_Util::mattrans(t.template data<T>(), this->naos_, this->pmat_); }
            this->dmd_buf_ = cal_dm_diff_pblas(X, this->pX[is], psi_ks_is, this->pc_,
                this->naos_, this->nocc[is], this->nvirt[is], this->pmat_);
            for (auto& t : this->dmd_buf_) { LR_Util::matsym(t.template data<T>(), this->naos_, this->pmat_); }
#else
            this->dmx_buf_ = cal_dm_trans_blas(X, psi_ks_is, this->nocc[is], this->nvirt[is]);
            for (auto& t : this->dmx_buf_)
            {
                T* d = t.template data<T>();
                for (int u = 0;u < this->naos_;++u)
                {
                    for (int v = u + 1;v < this->naos_;++v) { std::swap(d[u * this->naos_ + v], d[v * this->naos_ + u]); }
                }
            }
            this->dmd_buf_ = cal_dm_diff_blas(X, psi_ks_is, this->naos_, this->nocc[is], this->nvirt[is]);
            for (auto& t : this->dmd_buf_) { LR_Util::matsym(t.template data<T>(), this->naos_); }
#endif
            for (int ik = 0;ik < this->nk;++ik)
            {
                this->DM_trans->set_DMK_pointer(ik, this->dmx_buf_[ik].template data<T>());
                this->DM_diff->set_DMK_pointer(ik, this->dmd_buf_[ik].template data<T>());
            }
        }

    private:
        const int naos_ = 1;
        const Parallel_2D& pc_;
        const Parallel_Orbitals& pmat_;
        const psi::Psi<T>& psi_ks_;
        std::vector<psi::Psi<T>> psi_ks_spin_;
        std::unique_ptr<elecstate::DensityMatrix<T, T>> DM_trans;
        std::unique_ptr<elecstate::DensityMatrix<T, T>> DM_diff;
        mutable std::vector<ct::Tensor> dmx_buf_;
        mutable std::vector<ct::Tensor> dmd_buf_;
        std::unique_ptr<OperatorGxcULR<T>> gxc_;
    };
}
