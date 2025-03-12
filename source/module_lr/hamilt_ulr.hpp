#pragma once
#include "module_hamilt_general/hamilt.h"
#include "module_elecstate/module_dm/density_matrix.h"
#include "module_lr/operator_casida/operator_lr_diag.h"
#include "module_lr/operator_casida/operator_lr_hxc.h"
#include "module_lr/utils/lr_util_print.h"
#ifdef __EXX
#include "module_lr/operator_casida/operator_lr_exx.h"
#endif
namespace LR
{
    /// Unristricted TDDFT (TDA) for open-shell systems
    /// The A matrix is diveded by 4 blocks: uu, ud, du, dd
    template<typename T>
    class HamiltULR
    {
    public:
        template<typename TGint>
        HamiltULR(std::string& xc_kernel,
            const int& nspin,
            const int& naos,
            const std::vector<int>& nocc,   ///< {up, down}
            const std::vector<int>& nvirt,   ///< {up, down}
            const UnitCell& ucell_in,
            const std::vector<double>& orb_cutoff,
            Grid_Driver& gd_in,
            const psi::Psi<T>& psi_ks_in,
            const ModuleBase::matrix& eig_ks,
#ifdef __EXX
            std::weak_ptr<Exx_LRI<T>> exx_lri_in,
            const double& exx_alpha,
#endif 
            TGint* gint_in,
            std::vector<std::shared_ptr<PotHxcLR>>& pot_in,
            const K_Vectors& kv_in,
            const std::vector<Parallel_2D>& pX_in,   ///< {up, down}
            const Parallel_2D& pc_in,
            const Parallel_Orbitals& pmat_in) :nocc(nocc), nvirt(nvirt), pX(pX_in), nk(kv_in.get_nks() / nspin),
            ldim(nk* pX[0].get_local_size() + nk * pX[1].get_local_size()),
            gdim(nk* std::inner_product(nocc.begin(), nocc.end(), nvirt.begin(), 0))
        {
            ModuleBase::TITLE("HamiltULR", "HamiltULR");
            this->DM_trans = LR_Util::make_unique<elecstate::DensityMatrix<T, T>>(&pmat_in, 1, kv_in.kvec_d, nk);
            LR_Util::initialize_DMR(*this->DM_trans, pmat_in, ucell_in, gd_in, orb_cutoff);
            // this->DM_trans->init_DMR(&gd_in, &ucell_in); // too large due to not restricted by orb_cutoff
            this->ops.resize(4);

            this->ops[0] = new OperatorLRDiag<T>(eig_ks.c, pX_in[0], nk, nocc[0], nvirt[0]);
            this->ops[3] = new OperatorLRDiag<T>(eig_ks.c + nk * (nocc[0] + nvirt[0]), pX_in[1], nk, nocc[1], nvirt[1]);

            auto newHxc = [&](const int& sl, const int& sr) { return new OperatorLRHxc<T>(nspin, naos, nocc, nvirt, psi_ks_in,
                this->DM_trans, gint_in, pot_in[sl], ucell_in, orb_cutoff, gd_in, kv_in, pX_in, pc_in, pmat_in, { sl,sr }); };
            this->ops[0]->add(newHxc(0, 0));
            this->ops[1] = newHxc(0, 1);
            this->ops[2] = newHxc(1, 0);
            this->ops[3]->add(newHxc(1, 1));

#ifdef __EXX
            if (xc_kernel == "hf" || xc_kernel == "hse")
            {
                std::vector<psi::Psi<T>> psi_ks_spin = { LR_Util::get_psi_spin(psi_ks_in, 0, nk), LR_Util::get_psi_spin(psi_ks_in, 1, nk) };
                for (int is : {0, 1})
                {
                    this->ops[(is << 1) + is]->add(new OperatorLREXX<T>(nspin, naos, nocc[is], nvirt[is], ucell_in, psi_ks_spin[is],
                        this->DM_trans, exx_lri_in, kv_in, pX_in[is], pc_in, pmat_in,
                        xc_kernel == "hf" ? 1.0 : exx_alpha));
                }
            }
#endif

            this->cal_dm_trans = [&, this](const int& is, const T* const X)->void
                {
                    const auto psi_ks_is = LR_Util::get_psi_spin(psi_ks_in, is, nk);
                    // LR_Util::print_value(X, pX_in[is].get_local_size());
#ifdef __MPI
                    std::vector<ct::Tensor>  dm_trans_2d = cal_dm_trans_pblas(X, pX[is], psi_ks_is, pc_in, naos, nocc[is], nvirt[is], pmat_in, (T)1.0 / (T)nk);
                    if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat_in);
#else
                    std::vector<ct::Tensor>  dm_trans_2d = cal_dm_trans_blas(X, psi_ks_is, nocc[is], nvirt[is], (T)1.0 / (T)nk);
                    if (this->tdm_sym) for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
#endif
                    // LR_Util::print_tensor<T>(dm_trans_2d[0], "DMtrans(k=0)", &pmat_in);
                    // tensor to vector, then set DMK
                    for (int ik = 0;ik < nk;++ik) { this->DM_trans->set_DMK_pointer(ik, dm_trans_2d[ik].data<T>()); }
                };
        }
        ~HamiltULR()
        {
            for (auto& op : ops) { delete op; }
        }
        void hPsi(const T* const psi_in, T* const hpsi, const int ld_psi, const int& nband) const
        {
            ModuleBase::TITLE("HamiltULR", "hPsi");
            assert(ld_psi == this->ldim);
            const std::vector<int64_t> xdim_is = { nk * pX[0].get_local_size(), nk * pX[1].get_local_size() };
            /// band-wise act (also works for close-shell, but not efficient)
            for (int ib = 0;ib < nband;++ib)
            {
                const int offset_band = ib * ld_psi;
                for (int is_bj : {0, 1})
                {
                    const int offset_bj = offset_band + is_bj * xdim_is[0];
                    cal_dm_trans(is_bj, psi_in + offset_bj);   // calculate transition density matrix here
                    for (int is_ai : {0, 1})
                    {
                        const int offset_ai = offset_band + is_ai * xdim_is[0];
                        hamilt::Operator<T>* node(this->ops[(is_ai << 1) + is_bj]);
                        while (node != nullptr)
                        {
                            node->act(/*nband=*/1, xdim_is[is_bj], /*npol=*/1, psi_in + offset_bj, hpsi + offset_ai);
                            node = (hamilt::Operator<T>*)(node->next_op);
                        }
                    }
                }
            }
        }
        std::vector<T> matrix()const
        {
            ModuleBase::TITLE("HamiltULR", "matrix");
            const std::vector<int> npairs = { this->nocc[0] * this->nvirt[0], this->nocc[1] * this->nvirt[1] };
            const std::vector<int64_t> ldim_is = { nk * pX[0].get_local_size(), nk * pX[1].get_local_size() };
            const std::vector<int> gdim_is = { nk * npairs[0], nk * npairs[1] };
            std::vector<T> Amat_full(gdim * gdim);
            for (int is_bj : {0, 1})
            {
                const int no = this->nocc[is_bj];
                const int nv = this->nvirt[is_bj];
                const auto& px = this->pX[is_bj];
                const int loffset_bj = is_bj * ldim_is[0];
                const int goffset_bj = is_bj * gdim_is[0];
                for (int ik_bj = 0;ik_bj < nk;++ik_bj)
                {
                    for (int j = 0;j < no;++j)
                    {
                        for (int b = 0;b < nv;++b)
                        {
                            const int gcol = goffset_bj + ik_bj * npairs[is_bj] + j * nv + b;//global
                            std::vector<T> X_bj(this->ldim, T(0));
                            const int lj = px.global2local_col(j);
                            const int lb = px.global2local_row(b);
                            const int lcol = loffset_bj + ik_bj * px.get_local_size() + lj * px.get_row_size() + lb;//local
                            if (px.in_this_processor(b, j)) { X_bj[lcol] = T(1); }
                            this->cal_dm_trans(is_bj, X_bj.data() + loffset_bj);
                            std::vector<T> Aloc_col(this->ldim, T(0)); // a col of A matrix (local)
                            for (int is_ai : {0, 1})
                            {
                                const int goffset_ai = is_ai * gdim_is[0];
                                const int loffset_ai = is_ai * ldim_is[0];
                                const auto& pax = this->pX[is_ai];
                                hamilt::Operator<T>* node(this->ops[(is_ai << 1) + is_bj]);
                                while (node != nullptr)
                                {
                                    node->act(1, ldim_is[is_bj], /*npol=*/1, X_bj.data() + loffset_bj, Aloc_col.data() + loffset_ai);
                                    node = (hamilt::Operator<T>*)(node->next_op);
                                }
#ifdef __MPI
                                for (int ik_ai = 0;ik_ai < this->nk;++ik_ai)
                                {
                                    LR_Util::gather_2d_to_full(pax, Aloc_col.data() + loffset_ai + ik_ai * pax.get_local_size(),
                                        Amat_full.data() + gcol * gdim /*col, bj*/ + goffset_ai + ik_ai * npairs[is_ai]/*row, ai*/,
                                        false, nv, no);
                                }
#else
                                std::memcpy(Amat_full.data() + gcol * gdim + goffset_ai, Aloc_col.data() + goffset_ai, gdim_is[is_ai] * sizeof(T));
#endif
                            }
                        }
                    }
                }
            }
            std::cout << "Full A matrix:" << std::endl;
            LR_Util::print_value(Amat_full.data(), gdim, gdim);
            return Amat_full;
        }

        /// copy global data (eigenvectors) to local memory
        void global2local(T* lvec, const T* gvec, const int& nband) const
        {
            const std::vector<int> npairs = { this->nocc[0] * this->nvirt[0], this->nocc[1] * this->nvirt[1] };
            const std::vector<int64_t> ldim_is = { nk * pX[0].get_local_size(), nk * pX[1].get_local_size() };
            const std::vector<int> gdim_is = { nk * npairs[0], nk * npairs[1] };
            for (int ib = 0;ib < nband;++ib)
            {
                const int loffset_b = ib * this->ldim;
                const int goffset_b = ib * this->gdim;
                for (int is : {0, 1})
                {
                    const int loffset_bs = loffset_b + is * ldim_is[0];
                    const int goffset_bs = goffset_b + is * gdim_is[0];
                    for (int ik = 0;ik < nk;++ik)
                    {
                        const int loffset = loffset_bs + ik * pX[is].get_local_size();
                        const int goffset = goffset_bs + ik * npairs[is];
                        for (int lo = 0;lo < pX[is].get_col_size();++lo)
                        {
                            const int go = pX[is].local2global_col(lo);
                            for (int lv = 0;lv < pX[is].get_row_size();++lv)
                            {
                                const int gv = pX[is].local2global_row(lv);
                                lvec[loffset + lo * pX[is].get_row_size() + lv] = gvec[goffset + go * nvirt[is] + gv];
                            }
                        }
                    }
                }
            }
        }

    private:
        const std::vector<int>& nocc;
        const std::vector<int>& nvirt;

        const std::vector<Parallel_2D>& pX;


        const int nk = 1;
        const int ldim = 1;
        const int gdim = 1;

        /// 4 operator lists: uu, ud, du, dd
        std::vector<hamilt::Operator<T>*> ops;

        /// transition density matrix in AO representation
        /// Hxc only: size=1, calculate on the same address for each bands
        /// Hxc+Exx: size=nbands, store the result of each bands for common use
        std::unique_ptr<elecstate::DensityMatrix<T, T>> DM_trans;

        std::function<void(const int&, const T* const)> cal_dm_trans;
        const bool tdm_sym = false;     ///< whether to symmetrize the transition density matrix
    };
}