#pragma once
#include "pot_grad_xc.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_lr/dm_trans/dm_trans.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_hcontainer.h"
#include "source_lcao/module_lr/ao_to_mo_transformer/ao_to_mo.h"
#include "source_hamilt/module_gint/gint_interface.h"
#include "source_hamilt/module_hcontainer/hcontainer_funcs.h"

namespace LR
{
    /// @brief The open-shell $g^{xc}$ term, i.e.
    ///     $\sum_{\kappa\lambda\sigma'}\sum_{\alpha\beta\sigma''}D^X_{\sigma'}D^X_{\sigma''}
    ///      g^{xc}_{\kappa\lambda\sigma',\alpha\beta\sigma'',pq\tau}$
    /// projected onto the MO block selected by `mo_type` (VO for the Z-vector right-hand side,
    /// OO for $W^c$).
    ///
    /// This cannot be an `OperatorLRHxc` block like the other kernel terms: it is quadratic in
    /// $D^X$ rather than bilinear in a (out, in) spin pair, so BOTH transition-density channels
    /// have to be on the grid at the same time while the free spin $\tau$ is the output index.
    /// Hence the density is built once for the whole X vector and the potential is evaluated
    /// twice, once per $\tau$.
    ///
    /// Only the real (gamma-only) path is implemented, matching the rest of the LR gradient --
    /// `solve_Z_CG` has no complex version either.
    template<typename T>
    class OperatorGxcULR
    {
    public:
        OperatorGxcULR(const KernelXC& kxc,
            const ModulePW::PW_Basis& rho_basis,
            const UnitCell& ucell,
            const std::vector<double>& orb_cutoff,
            const Grid_Driver& gd,
            const K_Vectors& kv,
            const Parallel_Orbitals& pmat,
            const Parallel_2D& pc,
            const psi::Psi<T>& psi_ks,
            const std::vector<int>& nocc,
            const std::vector<int>& nvirt,
            const int naos,
            const std::vector<Parallel_2D>& pX,       ///< layout of the X blocks (always VO)
            const std::vector<Parallel_2D>& pout,     ///< layout of the output blocks (VO or OO)
            const LR_Util::MO_TYPE mo_type,
            const T factor)
            : pot_grad_(kxc, rho_basis, ucell, rho_basis.nrxx, /*triplet=*/false),
            ucell_(ucell), gd_(gd), kv_(kv), pmat_(pmat), pc_(pc), psi_ks_(psi_ks),
            nocc_(nocc), nvirt_(nvirt), naos_(naos), pX_(pX), pout_(pout),
            orb_cutoff_(orb_cutoff), mo_type_(mo_type), factor_(factor),
            nk_(kv.get_nks() / PARAM.inp.nspin), nrxx_(rho_basis.nrxx)
        {
            for (int is : {0, 1}) { this->psi_spin_.push_back(LR_Util::get_psi_spin(psi_ks, is, this->nk_)); }
            this->hR_ = LR_Util::make_unique<hamilt::HContainer<double>>(&pmat);
            LR_Util::initialize_HR<T, double>(*this->hR_, ucell, gd, orb_cutoff);
        }

        /// @brief `out += factor * <p tau| v^{(2)}_tau[X, X] |q tau>` for both spin channels
        void act(const T* const X_full, T* const out_full) const
        {
            ModuleBase::TITLE("OperatorGxcULR", "act");
            ModuleBase::timer::start("OperatorGxcULR", "act");
            const std::vector<int> off_x = { 0, nk_ * pX_[0].get_local_size() };
            const std::vector<int> off_out = { 0, nk_ * pout_[0].get_local_size() };

            // 1. the two transition-density channels, on the grid together
            std::vector<std::vector<ct::Tensor>> dmk(2);
            for (int is : {0, 1})
            {
#ifdef __MPI
                dmk[is] = cal_dm_trans_pblas(X_full + off_x[is], pX_[is], psi_spin_[is], pc_,
                    naos_, nocc_[is], nvirt_[is], pmat_);
                for (auto& t : dmk[is]) { LR_Util::matsym(t.template data<T>(), naos_, pmat_); }
#else
                dmk[is] = cal_dm_trans_blas(X_full + off_x[is], psi_spin_[is], nocc_[is], nvirt_[is]);
                for (auto& t : dmk[is]) { LR_Util::matsym(t.template data<T>(), naos_); }
#endif
            }
            elecstate::DensityMatrix<T, double> dm = LR_Util::build_dm_from_dmk_spin<T, double>(dmk,
                pmat_, nk_, kv_.kvec_d, ucell_, gd_, orb_cutoff_);

            double** rho1 = nullptr;
            LR_Util::_allocate_2order_nested_ptr(rho1, 2, nrxx_);
            for (int is : {0, 1}) { ModuleBase::GlobalFunc::ZEROS(rho1[is], nrxx_); }
            ModuleGint::cal_gint_rho(dm.get_DMR_vector(), 2, rho1, false);

            // 2. one potential per output spin, then AO -> MO.
            // The V(R) container is real, so for complex T only the gamma-only case is right --
            // the multi-k real/imag split that `OperatorLRHxc` does is not reproduced here. That
            // is not a new restriction: `solve_Z_CG` has no complex version either.
            for (int tau : {0, 1})
            {
                ModuleBase::matrix v2(1, nrxx_);   // zero-initialized
                pot_grad_.cal_v_eff_openshell(rho1, ucell_, v2, tau);

                this->hR_->set_zero();
                ModuleGint::cal_gint_vl(v2.c, this->hR_.get());
                std::vector<ct::Tensor> v_2d(nk_, LR_Util::newTensor<T>({ pmat_.get_col_size(), pmat_.get_row_size() }));
                for (auto& v : v_2d) { v.zero(); }
                const int nrow = ModuleBase::GlobalFunc::IS_COLUMN_MAJOR_KS_SOLVER(PARAM.inp.ks_solver)
                    ? pmat_.get_row_size() : pmat_.get_col_size();
                for (int ik = 0;ik < nk_;++ik)
                {
                    folding_HR(*this->hR_, v_2d[ik].template data<T>(), kv_.kvec_d[ik], nrow, 1);
                }
#ifdef __MPI
                ao_to_mo_pblas(v_2d, pmat_, psi_spin_[tau], pc_, naos_, nocc_[tau], nvirt_[tau],
                    pout_[tau], out_full + off_out[tau], /*add_on=*/true, mo_type_, factor_);
#else
                ao_to_mo_blas(v_2d, psi_spin_[tau], nocc_[tau], nvirt_[tau],
                    out_full + off_out[tau], /*add_on=*/true, mo_type_, factor_);
#endif
            }
            LR_Util::_deallocate_2order_nested_ptr(rho1, 2);
            ModuleBase::timer::end("OperatorGxcULR", "act");
        }

    private:
        PotGradXCLR pot_grad_;
        const UnitCell& ucell_;
        const Grid_Driver& gd_;
        const K_Vectors& kv_;
        const Parallel_Orbitals& pmat_;
        const Parallel_2D& pc_;
        const psi::Psi<T>& psi_ks_;
        std::vector<psi::Psi<T>> psi_spin_;
        const std::vector<int>& nocc_;
        const std::vector<int>& nvirt_;
        const int naos_ = 1;
        const std::vector<Parallel_2D>& pX_;
        const std::vector<Parallel_2D>& pout_;
        std::vector<double> orb_cutoff_;
        const LR_Util::MO_TYPE mo_type_ = LR_Util::VO;
        const T factor_ = T(1);
        const int nk_ = 1;
        const int nrxx_ = 1;
        std::unique_ptr<hamilt::HContainer<double>> hR_;
    };
}
