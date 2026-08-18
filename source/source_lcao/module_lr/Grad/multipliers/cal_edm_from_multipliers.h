#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/module_lr/dm_trans/dm_trans.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "source_lcao/module_lr/utils/lr_util_print.h"
#include "cal_multiplier_w_from_z.h"
#include <ATen/ops/linalg_op.h>
#ifdef __EXX
#include "source_lcao/module_lr/operator_casida/operator_lr_exx.h"
#endif
namespace LR
{

    // $X_{\mu i}=\sum_a c_{\mu a} X_{ai}$
    template <typename T>
    void cal_X_ao_occ(const T* const X, const Parallel_2D& px,
        const T* const c, const Parallel_2D& pc,
        T* const X_ao_occ, const Parallel_2D& px_ao_occ);

    // D=X1*X2^T
    // $D_{\mu\nu} = \sum_i X1_{\mu i}X2_{\nu i}$
    template <typename T>
    void matdot(const T* const vec1, const T* const vec2, const Parallel_2D& pvec,
        T* const dm, const Parallel_2D& pmat);

    template<typename T>
    void multiply_eig_onto_vec(const T* const vec, const double* const eig, const Parallel_2D& pvec, T* const evec)
    {
        for (int i = 0;i < pvec.get_col_size();++i)
        {
            const int gi = pvec.local2global_col(i);
            for (int j = 0;j < pvec.get_row_size();++j)
            {
                const int idx = i * pvec.get_row_size() + j;
                evec[idx] = vec[idx] * eig[gi];
            }
        }
    }

    template<typename T>
    ct::Tensor cal_edm_single_kpoint(const T* const vec, const Parallel_2D& pvec, const double* eig, const Parallel_2D& pmat)
    {
        const int nocc = pvec.get_global_col_size();
        const int naos = pvec.get_global_row_size();
        std::vector<T> eig_times_vec(pvec.get_local_size());
        multiply_eig_onto_vec(vec, eig, pvec, eig_times_vec.data());
        ct::Tensor edm_result = LR_Util::newTensor<T>({ pmat.get_col_size(), pmat.get_row_size() });
        matdot(eig_times_vec.data(), vec, pvec, edm_result.data<T>(), pmat);
        return edm_result;
    }

    template <typename T>
    std::vector<ct::Tensor> cal_edm_term4(const T* X,
        const double eig_ext_istate,    //1, the excitation energy of one state
        const double* const eig_ks,     // gocc+gvirt
        const psi::Psi<T>& c,
        const Parallel_2D& px,
        const Parallel_2D& pc, 
        const Parallel_Orbitals& pmat)
    {
        const int& naos = pmat.get_global_row_size();
        const int& nocc = px.get_global_col_size();
        const int& nvirt = px.get_global_row_size();
        // 4. $\sum_i (\Omega + \epsilon_i) \sum_{ab} C_{\mu a} X_{ia} C_{\nu b} X_{ib}$
        std::vector<ct::Tensor> edm(c.get_nk());
        Parallel_2D px_ao_occ;
        LR_Util::setup_2d_division(px_ao_occ, px.get_block_size(), naos, nocc, px.blacs_ctxt);
        for (int ik = 0;ik < c.get_nk();++ik)
        {
            const int idx_X = ik * px.get_local_size();
            std::vector<T> X_ao_occ(px_ao_occ.get_local_size());
            cal_X_ao_occ(X + idx_X, px, &c(ik, 0, 0), pc, X_ao_occ.data(), px_ao_occ);
            std::vector<double> eig_ks_plus_ext(nocc, 0.0);
            const int idx_eig_ks = ik * (nocc + nvirt);
            // for (int i = 0;i < nocc;++i) { eig_ks_plus_ext[i] = eig_ks[idx_eig_ks + i] + eig_ext_istate; }
            std::transform(eig_ks + idx_eig_ks, eig_ks + idx_eig_ks + nocc, eig_ks_plus_ext.begin(), [eig_ext_istate](double x) {return x + eig_ext_istate;});
            edm[ik] = cal_edm_single_kpoint(X_ao_occ.data(), px_ao_occ, eig_ks_plus_ext.data(), pmat);
        }
        return edm;
    }

    // calculate the excited state energy density matrix (for multiplying the overlap gradient in the gradient of lagrangian)
    // multi-k has not been supported yet
    template<typename T>
    std::vector<ct::Tensor> cal_edm_terms_from_XZWK(
        const T* const X,    //lvirt*locc
        const T* const Z,   //lvirt*locc
        const T* const W,   //locc*locc
        const T* const K_cvcx,   //lvirt*locc
        const double eig_ext_istate,    //1, the excitation energy of one state
        const double* const eig_ks,     // gocc+gvirt
        const psi::Psi<T>& c,
        const int nspin,
        const Parallel_2D& p_occ_occ,
        const Parallel_2D& p_virt_occ,
        const Parallel_2D& pc,
        const Parallel_Orbitals& pmat)
    {
        const int& naos = pmat.get_global_col_size();
        const int& nocc = p_occ_occ.get_global_col_size();
        const int& nvirt = p_virt_occ.get_global_row_size();
        const Parallel_2D& px = p_virt_occ;

        // 1. c * W * c
        const std::vector<ct::Tensor> cWc = cal_dm_trans_pblas(W, p_occ_occ, c, pc, naos, nocc, nvirt, pmat, (T)1., LR_Util::MO_TYPE::OO);

        // 2. edm of Z : $\sum_i \sum_a c_a epsilon_i Z_{ai} c_i
        std::vector<T> epsi_Z(px.get_local_size());
        multiply_eig_onto_vec(Z, eig_ks, px, epsi_Z.data());
        std::vector<ct::Tensor> cZc = cal_dm_trans_pblas(Z, px, c, pc, naos, nocc, nvirt, pmat);
        std::for_each(cZc.begin(), cZc.end(), [&](ct::Tensor& s) { LR_Util::matsym(s.data<T>(), naos, pmat); });

        //3. c * K_cvcx * c
        std::vector<ct::Tensor> cKc = cal_dm_trans_pblas(K_cvcx, px, c, pc, naos, nocc, nvirt, pmat, (T)2.0);
        std::for_each(cKc.begin(), cKc.end(), [&](ct::Tensor& s) { LR_Util::matsym(s.data<T>(), naos, pmat); });

        // 4. $\sum_i (\Omega + \epsilon_i) \sum_{ab} C_{\mu a} X_{ia} C_{\nu b} X_{ib}$
        const std::vector<ct::Tensor> edm = cal_edm_term4(X, eig_ext_istate, eig_ks, c, px, pc, pmat);

        if (PARAM.inp.test_force)
        {
            std::cout << "cWc: " << std::endl;
            LR_Util::print_value(cWc[0].data<T>(), pmat.get_col_size(), pmat.get_row_size());
            std::cout << "cZc: " << std::endl;
            LR_Util::print_value(cZc[0].data<T>(), pmat.get_col_size(), pmat.get_row_size());
            std::cout << "cKc: " << std::endl;
            LR_Util::print_value(cKc[0].data<T>(), pmat.get_col_size(), pmat.get_row_size());
            std::cout << "edm term 4: " << std::endl;
            LR_Util::print_value(edm[0].data<T>(), pmat.get_col_size(), pmat.get_row_size());
        }
        return edm + cWc + cZc + cKc;
    }

    template<typename T>
    std::vector<ct::Tensor> cal_edm_from_XZ_istate( //for one excited state
        const T* const X,   //lvirt*locc
        const T* const Z,   //lvirt*locc
        const double eig_ext_istate,    //1, the excitation energy of one state
        const double* const eig_ks,     // gocc+gvirt
        const elecstate::DensityMatrix<T, T>& dm_trans, // D_X
        const psi::Psi<T>& c,
        const int& nspin,
        const int& naos,
        const std::vector<int>& nocc,
        const std::vector<int>& nvirt,
        const UnitCell& ucell,
        const std::vector<double>& orb_cutoff,
#ifdef __EXX
        std::weak_ptr<Exx_LRI<T>> exx_lri,
        const double& exx_alpha,
#endif 
        std::weak_ptr<PotHxcLR> pot,
        std::weak_ptr<PotHxcLR> pot_hxc_gs,
        const K_Vectors& kv,
        const Grid_Driver& gd,
        const std::vector<Parallel_2D>& px,
        const Parallel_2D& pc,
        const Parallel_Orbitals& pmat,
        const std::string xc_kernel,
        const std::string& spin_type = "singlet")
    {
        const int nk = kv.get_nks() / nspin;
        // 1. calculate W multiplier 
        std::vector<Parallel_2D> p_occ_occ(nspin);
        for (int is = 0;is < nspin;++is) { LR_Util::setup_2d_division(p_occ_occ[is], 1, nocc[is], nocc[is], px[is].blacs_ctxt); }
        std::vector<T> W(p_occ_occ[0].get_local_size() * nk, 0.0);
        cal_W_from_Z(W.data(), Z, X, eig_ext_istate, eig_ks, nspin, naos, nocc, nvirt,
            ucell, orb_cutoff, gd, c,
#ifdef __EXX
            exx_lri, exx_alpha,
#endif
            pot_hxc_gs, kv, px, pc, p_occ_occ, pmat, xc_kernel, spin_type);
        std::cout << "W: " << std::endl;
        LR_Util::print_value(W.data(), nk, p_occ_occ[0].get_col_size(), p_occ_occ[0].get_row_size());

        // 2. build K_cvcx (nvirt*nocc) = \sum_i X_{ia} K_{ij} = \sum_i X_{ia} \sum_{\mu\nu} c_{\mu i} c_{\nu j} K_{\mu\nu}[D^X]
        // $2\sum_i X_{ai} K_{ij}[D_X]$ (D_X is symmetrized)
        OperatorLRHxc<T> op_K_cvcx(nspin, naos, nocc, nvirt, c,
            dm_trans, pot, ucell, orb_cutoff, gd, kv, px, pc, pmat,
            { 0 }, T(2.0), OperatorLRHxc<T>::MO_TO_AO_TYPE::CXC_o);
#ifdef __EXX
        OperatorLREXX<T> op_K_exx(nspin, naos, nocc[0], nvirt[0], ucell, c,
            dm_trans, exx_lri, kv, px[0], pc, pmat,
            2.0 * exx_alpha, OperatorLREXX<T>::MO_TO_AO_TYPE::CXC_o);
#endif
        const int ld_vo = nk * px[0].get_local_size();
        std::vector<T> K_cvcx(ld_vo, 0.0);
        op_K_cvcx.act(/*nbands=*/1, ld_vo, /*npol=*/1, X, K_cvcx.data());
        if (LR::exx_kernel_list().count(xc_kernel))
            op_K_exx.act(/*nbands=*/1, ld_vo, /*npol=*/1, X, K_cvcx.data());

        return cal_edm_terms_from_XZWK(X, Z, W.data(), K_cvcx.data(), eig_ext_istate, eig_ks, c, nspin, p_occ_occ[0], px[0], pc, pmat);
    }
}