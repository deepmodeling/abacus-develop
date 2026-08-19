#pragma once
#include "source_hamilt/hamilt.h"
#include "source_estate/module_dm/density_matrix.h"
#include "source_lcao/module_lr/Grad/xc/pot_grad_xc.h"
#include "source_lcao/module_lr/potentials/pot_hxc_lrtd.h"
#include "source_lcao/module_lr/operator_casida/operator_lr_hxc.h"
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
                { 0 }, -2.0, ATYPE::CXC);
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
            if (exx_kernel_list().count(PARAM.inp.dft_functional))
            {
                hamilt::Operator<T>* op_ht_exx = new OperatorLREXX<T>(nspin, naos, nocc[0], nvirt[0], ucell, psi_ks,
                    *this->DM_diff, exx_lri, kv, pX[0], pc, pmat,
                    -2.0 * exx_alpha, //alpha; H=2K when D is symmetrized
                    ATYPE_EXX::CC_vo, {}, hamilt::calculation_type::lr_dmdiff_exx);
                this->ops->add(op_ht_exx);
            }
#endif

            // 3. $2\sum_{jb,kc} g^{xc}_{ia, jb, kc}X_{jb}X_{kc}$
            // singlet only: $K^T$ has no Hxc part, so the triplet Z-vector equation carries no
            // $g^{xc}$ term (see LR-Grad-formulas/LR-Grad-Zvector-Singlet-Triplet.md, "Z-Vector方程与乘子").
            if (LR_Util::has_local_xc(xc_kernel) && spin_type != "triplet")
            {            // !!  op_gxc has some bug now
                this->pot_grad = std::make_shared<PotGradXCLR>(pot.lock()->xc_kernel_components, pot.lock()->get_rho_basis(), ucell, pot.lock()->nrxx);
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

            this->cal_dm_trans = [&, this](const int& is, const T* X)->void
                {
                    const auto psi_ks_is = LR_Util::get_psi_spin(psi_ks, is, this->nk);
#ifdef __MPI
                    std::vector<ct::Tensor> dm_trans_2d = cal_dm_trans_pblas(X, this->pX[is], psi_ks_is, pc, naos, nocc[is], nvirt[is], pmat);
                    for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos, pmat);
#else
                    std::vector<ct::Tensor> dm_trans_2d = cal_dm_trans_blas(X, psi_ks_is, nocc[is], nvirt[is]);
                    for (auto& t : dm_trans_2d) LR_Util::matsym(t.data<T>(), naos);
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
}