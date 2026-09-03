#pragma once
#include "source_hamilt/hamilt.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_lr/Grad/xc/pot_grad_xc.h"
#include "source_lcao/module_lr/Grad/xc/operator_gxc_ulr.h"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_hxc.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include <ATen/ops/linalg_op.h>
#ifdef __EXX
#include "source_lcao/module_lr/operator_casida/operator_lr_exx.h"
#endif
#include "source_base/module_external/scalapack_connector.h"

namespace LR
{
    // $\sum_a X^*_{ai} X_{aj}(\Omega - \epsilon_a)$
    template<typename T>
    void add_ediff_term(T* const inout,
        const T* const X,
        const double eig_ext_istate,
        const double* const eig_ks,
        const int nk,
        const int nocc,
        const int nvirt,
        const Parallel_2D& px,
        const Parallel_2D& p_occ_occ)
    {
        for (int ik = 0;ik < nk;++ik)
        {
            // 1. calculate (\Omega-\epsilon_a)X_{aj}
            std::vector<T> wX(px.get_local_size());
            const int eigks_start_k = ik * (nocc + nvirt);
            const int x_start_k = ik * px.get_local_size();
            const int inout_start_k = ik * p_occ_occ.get_local_size();
            for (int la = 0;la < px.get_row_size();++la)
            {
                const int ga = px.local2global_row(la);
                const double weight = eig_ext_istate - eig_ks[eigks_start_k + nocc + ga];
                std::cout << "la=" << la << ", ks-eig=" << eig_ks[eigks_start_k + nocc + ga] << ", weight=" << weight << ", X2=" << X[x_start_k + la] << std::endl;
                for (int li = 0;li < px.get_col_size();++li)
                {
                    const int idx = li * px.get_row_size() + la;
                    wX[idx] = weight * X[x_start_k + idx];
                }
            }
            // 2. matrix multiplication (parallel)
            const int i1 = 1;
            ScalapackConnector::gemm('C', 'N', nocc, nocc, nvirt,
                T(1.0), X + x_start_k, i1, i1, px.desc,
                wX.data(), i1, i1, px.desc,
                T(1.0)/*add-on*/, inout + inout_start_k, i1, i1, p_occ_occ.desc);
        }
    }

    template<typename T>
    void cal_W_from_Z(T* const W,
        const T* const Z,
        const T* const X,
        const double eig,
        const double* const eig_ks,
        const int& nspin,
        const int& naos,
        const std::vector<int>& nocc,
        const std::vector<int>& nvirt,
        const UnitCell& ucell,
        const std::vector<double>& orb_cutoff,
        const Grid_Driver& gd,
        const psi::Psi<T>& psi_ks,

#ifdef __EXX
        std::weak_ptr<Exx_LRI<T>> exx_lri,
        const double& exx_alpha,
#endif 
        std::weak_ptr<PotHxcLR> pot_hxc_gs,
        const K_Vectors& kv,
        const std::vector<Parallel_2D>& px,
        const Parallel_2D& pc,
        const std::vector<Parallel_2D>& p_occ_occ,   // < for W
        const Parallel_Orbitals& pmat,
        const std::string xc_kernel,
        const std::string& spin_type = "singlet")
    {
        ModuleBase::TITLE("cal_W_from_Z", "cal_W_from_Z");
        using ATYPE = typename OperatorLRHxc<T>::MO_TO_AO_TYPE;
#ifdef __EXX
        using ATYPE_EXX = typename OperatorLREXX<T>::MO_TO_AO_TYPE;
#endif
        const int nk = kv.get_nks() / nspin;
        // allocate memory for DMs
        elecstate::DensityMatrix<T, T> DM_trans(&pmat, 1, kv.kvec_d, nk);   //DX
        LR_Util::initialize_DMR(DM_trans, pmat, ucell, gd, orb_cutoff);
        elecstate::DensityMatrix<T, T> DM_diff_relaxed(&pmat, 1, kv.kvec_d, nk);    //T+DZ
        LR_Util::initialize_DMR(DM_diff_relaxed, pmat, ucell, gd, orb_cutoff);
        /// operators
        // 1. 0.5$H_ij[T+Z]$, equals to $K_ij[T+Z]$ when $(T+Z)$ is symmetrized
        // Note that K_Hxc(singlet) = 2* pot_hxc_gs, that's why there's factor 2 here.
        OperatorLRHxc<T> op_ht(nspin, naos, nocc, nvirt, psi_ks,
            DM_diff_relaxed, pot_hxc_gs, ucell, orb_cutoff, gd, kv, p_occ_occ, pc, pmat,
            { 0 }, T(2.0), ATYPE::CC_oo);
#ifdef __EXX
        OperatorLREXX<T> op_ht_exx(nspin, naos, nocc[0], nvirt[0], ucell, psi_ks,
            DM_diff_relaxed, exx_lri, kv, p_occ_occ[0], pc, pmat,
            exx_alpha, ATYPE_EXX::CC_oo);
#endif
        // 2. $2\sum_{jb,kc} g^{xc}_{ia, jb, kc}X_{jb}X_{kc}$
        // use pointer here for polymorphism
        // but `weak_ptr = make_shared()` will cause a segment fault because the shared_ptr is a temporary object
        // correct way is to use `shared_ptr = make_shared()` and then assign it to weak_ptr
        // `weak_ptr=shared_ptr` is automatically called in the constructor of OperatorLRHxc, so we don't need to do it manually
        // if `pot_grad` is passed into a function rather than a class, we need to write `weak_ptr=shared_ptr` explicitly
        std::shared_ptr<PotGradXCLR> pot_grad =
            std::make_shared<PotGradXCLR>(pot_hxc_gs.lock()->xc_kernel_components(), pot_hxc_gs.lock()->get_rho_basis(), 
            ucell, pot_hxc_gs.lock()->nrxx, spin_type == "triplet");
        OperatorLRHxc<T> op_gxc(nspin, naos, nocc, nvirt, psi_ks,
            DM_trans, pot_grad, ucell, orb_cutoff, gd, kv, p_occ_occ, pc, pmat,
            // Factor 1.0 according to the $W^c$ formula (`pot_grad` carries $2*g^{xc}$: uu+ud or uu-ud).
            // NOTE this factor is NOT shared with the Z-vector RHS `op_gxc` in `hamilt_zeq_right.h`:
            // that one is a different object and its original -2.0 is correct, as the formula and H2-DZP scan confirms.
            { 0 }, T(1.0), ATYPE::CC_oo);

        std::vector<ct::Tensor> dm_trans_2d, dm_diff_2d;
        auto cal_dm_trans = [&](const int is, const T* const x_ptr)->void //DX
            {
                const auto psi_ks_is = LR_Util::get_psi_spin(psi_ks, is, nk);
#ifdef __MPI
                dm_trans_2d = cal_dm_trans_pblas(x_ptr, px[is], psi_ks_is, pc, naos, nocc[is], nvirt[is], pmat);
                for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat);
#else
                dm_trans_2d = cal_dm_trans_blas(x_ptr, psi_ks_is, nocc[is], nvirt[is]);
                for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
                for (int ik = 0;ik < nk;++ik) { DM_trans.set_DMK_pointer(ik, dm_trans_2d[ik].data<T>()); }
            };
        auto cal_dm_diff_relaxed = [&](const int& is, const T* const x_ptr, const T* const z_ptr)->void  // T+DZ
            {
                const auto psi_ks_is = LR_Util::get_psi_spin(psi_ks, is, nk);
#ifdef __MPI
                std::vector<ct::Tensor> z_2d = cal_dm_trans_pblas(z_ptr, px[is], psi_ks_is, pc, naos, nocc[is], nvirt[is], pmat);
                for (auto& t : z_2d) LR_Util::matsym(t.data<T>(), naos, pmat);
#else
                std::vector<ct::Tensor> z_2d = cal_dm_trans_blas(z_ptr, psi_ks_is, nocc[is], nvirt[is]);
                for (auto& t : z_2d) LR_Util::matsym(t.data<T>(), naos);
#endif

#ifdef __MPI
                dm_diff_2d = cal_dm_diff_pblas(x_ptr, px[is], psi_ks_is, pc, naos, nocc[is], nvirt[is], pmat);
#else
                dm_diff_2d = cal_dm_diff_blas(x_ptr, psi_ks_is, naos, nocc[is], nvirt[is]);
#endif
                for (int ik = 0;ik < nk;++ik)
                {
                    dm_diff_2d[ik] = dm_diff_2d[ik] + z_2d[ik];
                    DM_diff_relaxed.set_DMK_pointer(ik, dm_diff_2d[ik].data<T>());
                }
            };

        // act the operators onto current state
        const int ld_vo = nk * px[0].get_local_size();
        const int ld_oo = nk * p_occ_occ[0].get_local_size();
        cal_dm_trans(0, X);  // transition density matrix DX
        cal_dm_diff_relaxed(0, X, Z);  // relaxed difference density matrix T+DZ
        // the 3 terms
        op_ht.act(/*nband=*/1, ld_oo, /*npol=*/1, X, W);    //comment out this line to test H[T+Z]=0
        // std::cout << "W (H[T+Z])) local terms: " << std::endl;
        // LR_Util::print_value(W, nk, p_occ_occ[0].get_col_size(), p_occ_occ[0].get_row_size());
        if (LR::gs_is_hybrid())  // H[T+Z] term depends on ground-state kernel (dft_functional)
            op_ht_exx.act(/*nband=*/1, ld_oo, /*npol=*/1, X, W);
        // std::cout << "W (H[T+Z])) local +exx terms: " << std::endl;
        // LR_Util::print_value(W, nk, p_occ_occ[0].get_col_size(), p_occ_occ[0].get_row_size());
        // Not singlet-only: $K^T_{xc}=f_{uu}-f_{ud}\ne0$ for a local functional, so $W^{c,T}$ has a
        // $g^{xc}$ term as well, built from the "-" spin combination (the `triplet` flag above).
        if (LR_Util::has_local_xc(xc_kernel))
            op_gxc.act(/*nband=*/1, ld_oo, /*npol=*/1, X, W);

        std::cout << "W (H[T+Z]) + W(gxc) terms: " << std::endl;
        LR_Util::print_value(W, nk, p_occ_occ[0].get_col_size(), p_occ_occ[0].get_row_size());
        add_ediff_term(W, X, eig, eig_ks, nk, nocc[0], nvirt[0], px[0], p_occ_occ[0]);
    }

    /// @brief Open-shell (spin-unrestricted) counterpart of `cal_W_from_Z`.
    ///
    /// $$W^c_{ij\sigma}=\tfrac12 H_{ij\sigma}[T+D^Z]
    ///   +\sum_a X_{ai\sigma}(\Omega-\epsilon_{a\sigma})X_{aj\sigma}
    ///   +\sum_{\kappa\lambda\sigma'}\sum_{\alpha\beta\sigma''}D^X D^X g^{xc}_{\dots,ij\sigma}$$
    ///
    /// `W[is]` is the occ-occ block of spin channel `is`, distributed over `p_occ_occ[is]`.
    /// Factor 1.0 on the $H[T+D^Z]$ operator (the closed-shell version uses 2.0): there
    /// $\tfrac12 H^S = K^S = 2\cdot$`pot_hxc_gs` because `pot_hxc_gs` is the halved `S2_gs`;
    /// here it is `S2_updown`, one $K_{\sigma\sigma'}$ component, and the $\sum_{\sigma'}$
    /// is done by the block loop below.
    template<typename T>
    void cal_W_from_Z_openshell(std::vector<std::vector<T>>& W,
        const T* const Z,
        const T* const X,
        const double eig,
        const double* const eig_ks,
        const int& nspin,
        const int& naos,
        const std::vector<int>& nocc,
        const std::vector<int>& nvirt,
        const UnitCell& ucell,
        const std::vector<double>& orb_cutoff,
        const Grid_Driver& gd,
        const psi::Psi<T>& psi_ks,
#ifdef __EXX
        std::weak_ptr<Exx_LRI<T>> exx_lri,
        const double& exx_alpha,
#endif
        std::weak_ptr<PotHxcLR> pot_hxc_gs,
        const K_Vectors& kv,
        const std::vector<Parallel_2D>& px,
        const Parallel_2D& pc,
        const std::vector<Parallel_2D>& p_occ_occ,
        const Parallel_Orbitals& pmat,
        const std::string xc_kernel)
    {
        ModuleBase::TITLE("cal_W_from_Z_openshell", "cal_W_from_Z_openshell");
        using ATYPE = typename OperatorLRHxc<T>::MO_TO_AO_TYPE;
#ifdef __EXX
        using ATYPE_EXX = typename OperatorLREXX<T>::MO_TO_AO_TYPE;
#endif
        const int nk = kv.get_nks() / nspin;
        const std::vector<int> ld_x = { nk * px[0].get_local_size(), nk * px[1].get_local_size() };
        const std::vector<int> ld_oo = { nk * p_occ_occ[0].get_local_size(), nk * p_occ_occ[1].get_local_size() };
        const std::vector<int> off_x = { 0, ld_x[0] };
        const int nband_window = nocc[0] + nvirt[0];   // common KS window, see `set_dimension`

        W.assign(2, {});
        for (int is : {0, 1}) { W[is].assign(ld_oo[is], T(0.0)); }

        elecstate::DensityMatrix<T, T> DM_diff_relaxed(&pmat, 1, kv.kvec_d, nk);   // T+D^Z of one channel
        LR_Util::initialize_DMR(DM_diff_relaxed, pmat, ucell, gd, orb_cutoff);

        // $\tfrac12 H_{ij\sigma}[T+D^Z]=K_{ij\sigma}[T+D^Z]$, one operator per (out, in) spin pair
        std::vector<std::unique_ptr<OperatorLRHxc<T>>> op_ht(4);
        for (int sl : {0, 1})
        {
            for (int sr : {0, 1})
            {
                op_ht[(sl << 1) + sr] = LR_Util::make_unique<OperatorLRHxc<T>>(nspin, naos, nocc, nvirt, psi_ks,
                    DM_diff_relaxed, pot_hxc_gs, ucell, orb_cutoff, gd, kv, p_occ_occ, pc, pmat,
                    std::vector<int>({ sl, sr }), T(1.0), ATYPE::CC_oo);
            }
        }
#ifdef __EXX
        std::vector<psi::Psi<T>> psi_ks_spin;
        for (int is : {0, 1}) { psi_ks_spin.push_back(LR_Util::get_psi_spin(psi_ks, is, nk)); }
        std::vector<std::unique_ptr<OperatorLREXX<T>>> op_ht_exx(2);
        const bool with_exx = LR::gs_is_hybrid();
        if (with_exx)
        {   // exchange is spin-diagonal
            for (int is : {0, 1})
            {
                op_ht_exx[is] = LR_Util::make_unique<OperatorLREXX<T>>(nspin, naos, nocc[is], nvirt[is],
                    ucell, psi_ks_spin[is], DM_diff_relaxed, exx_lri, kv, p_occ_occ[is], pc, pmat,
                    exx_alpha, ATYPE_EXX::CC_oo);
            }
        }
#endif
        // the relaxed difference density matrix of one spin channel
        std::vector<ct::Tensor> dm_buf;
        auto set_dm_diff_relaxed = [&](const int is)->void
            {
                const auto psi_ks_is = LR_Util::get_psi_spin(psi_ks, is, nk);
                const T* const x_ptr = X + off_x[is];
                const T* const z_ptr = Z + off_x[is];
#ifdef __MPI
                std::vector<ct::Tensor> z_2d = cal_dm_trans_pblas(z_ptr, px[is], psi_ks_is, pc, naos, nocc[is], nvirt[is], pmat);
                for (auto& t : z_2d) { LR_Util::matsym(t.data<T>(), naos, pmat); }
                dm_buf = cal_dm_diff_pblas(x_ptr, px[is], psi_ks_is, pc, naos, nocc[is], nvirt[is], pmat);
#else
                std::vector<ct::Tensor> z_2d = cal_dm_trans_blas(z_ptr, psi_ks_is, nocc[is], nvirt[is]);
                for (auto& t : z_2d) { LR_Util::matsym(t.data<T>(), naos); }
                dm_buf = cal_dm_diff_blas(x_ptr, psi_ks_is, naos, nocc[is], nvirt[is]);
#endif
                for (int ik = 0;ik < nk;++ik)
                {
                    dm_buf[ik] = dm_buf[ik] + z_2d[ik];
                    DM_diff_relaxed.set_DMK_pointer(ik, dm_buf[ik].data<T>());
                }
            };

        for (int is_in : {0, 1})
        {
            set_dm_diff_relaxed(is_in);
            for (int is_out : {0, 1})
            {
                op_ht[(is_out << 1) + is_in]->act(/*nband=*/1, ld_oo[is_out], /*npol=*/1,
                    X + off_x[is_in], W[is_out].data());
            }
#ifdef __EXX
            if (with_exx)
            {   // $\delta_{\sigma\sigma'}$: only the diagonal block contributes
                op_ht_exx[is_in]->act(/*nband=*/1, ld_oo[is_in], /*npol=*/1,
                    X + off_x[is_in], W[is_in].data());
            }
#endif
        }

        // $\sum_{\sigma'\sigma''}D^X_{\sigma'}D^X_{\sigma''}g^{xc}_{\dots,ij\tau}$, coefficient 1
        // in the spin-orbital formula. Quadratic in $D^X$, so it does not fit the block loop above.
        // The kernel comes from `pot_hxc_gs`, mirroring the closed-shell `cal_W_from_Z`; it only
        // differs from the LR kernel used on the Z-vector right-hand side when
        // `xc_kernel != dft_functional`, in which case both branches share the same ambiguity.
        if (LR_Util::has_local_xc(xc_kernel))
        {
            OperatorGxcULR<T> gxc(pot_hxc_gs.lock()->xc_kernel_components(), pot_hxc_gs.lock()->get_rho_basis(),
                ucell, orb_cutoff, gd, kv, pmat, pc, psi_ks, nocc, nvirt, naos,
                px, p_occ_occ, LR_Util::MO_TYPE::OO, T(1.0));
            std::vector<T> w_flat(ld_oo[0] + ld_oo[1], T(0.0));
            gxc.act(X, w_flat.data());
            for (int is : {0, 1})
            {
                const int off = is * ld_oo[0];
                for (int i = 0;i < ld_oo[is];++i) { W[is][i] += w_flat[off + i]; }
            }
        }

        // $\sum_a X_{ai\sigma}(\Omega-\epsilon_{a\sigma})X_{aj\sigma}$ -- spin-diagonal
        for (int is : {0, 1})
        {
            add_ediff_term(W[is].data(), X + off_x[is], eig, eig_ks + is * nk * nband_window,
                nk, nocc[is], nvirt[is], px[is], p_occ_occ[is]);
        }
    }
}
