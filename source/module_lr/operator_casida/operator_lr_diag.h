#pragma once
#include "module_lr/utils/lr_util.h"
#include "module_hsolver/kernels/math_kernel_op.h"
#include "module_hamilt_general/operator.h"
#ifdef __MPI
#include "module_base/parallel_common.h"
#endif
namespace LR
{
    /// @brief  Diag part of A operator: [AX]_iak = (e_ak - e_ik) X_iak
    template<typename T = double, typename Device = base_device::DEVICE_CPU>
    class OperatorLRDiag : public hamilt::Operator<T, Device>
    {
    public:
        OperatorLRDiag(const ModuleBase::matrix& eig_ks_in, const Parallel_2D* pX_in, const int& nk_in, const int& nocc_in, const int& nvirt_in)
            : eig_ks(eig_ks_in), pX(pX_in), nk(nk_in), nocc(nocc_in), nvirt(nvirt_in)
        {   // calculate the difference of eigenvalues
            ModuleBase::TITLE("OperatorLRDiag", "OperatorLRDiag");
#ifdef __MPI
            Parallel_Common::bcast_double(eig_ks.c, eig_ks.nr * eig_ks.nc);
#endif
            this->act_type = 2;
            this->cal_type = hamilt::calculation_type::no;
            this->eig_ks_diff.create(nk, pX->get_local_size(), false);
            for (int ik = 0;ik < nk;++ik)
                for (int io = 0;io < pX->get_col_size();++io)    //nocc_local
                    for (int iv = 0;iv < pX->get_row_size();++iv)    //nvirt_local
                    {
                        int io_g = pX->local2global_col(io);
                        int iv_g = pX->local2global_row(iv);
                        this->eig_ks_diff(ik, io * pX->get_row_size() + iv) = eig_ks(ik, nocc + iv_g) - eig_ks(ik, io_g);
                    }
        };
        void init(const int ik_in) override {};

        /// caution: put this operator at the head of the operator list,
        /// because vector_mul_vector_op directly assign to (rather than add on) psi_out.
        virtual void act(const psi::Psi<T>& psi_in, psi::Psi<T>& psi_out, const int nbands) const override
        {
            ModuleBase::TITLE("OperatorLRDiag", "act");
            assert(nbands <= psi_in.get_nbands());

            psi::Psi<T> psi_in_bfirst = LR_Util::k1_to_bfirst_wrapper(psi_in, this->nk, this->pX->get_local_size());
            psi::Psi<T> psi_out_bfirst = LR_Util::k1_to_bfirst_wrapper(psi_out, this->nk, this->pX->get_local_size());
            for (int ib = 0;ib < nbands;++ib)
            {
                psi_in_bfirst.fix_b(ib);
                psi_out_bfirst.fix_b(ib);
                hsolver::vector_mul_vector_op<T, Device>()(this->ctx,
                    psi_in_bfirst.get_nk() * psi_in_bfirst.get_nbasis(),
                    psi_out_bfirst.get_pointer(),
                    psi_in_bfirst.get_pointer(),
                    this->eig_ks_diff.c);
            }
        }
    private:
        const ModuleBase::matrix& eig_ks;
        const Parallel_2D* pX;
        ModuleBase::matrix eig_ks_diff;
        const int& nk;
        const int& nocc;
        const int& nvirt;
        Device* ctx = {};
    };
}