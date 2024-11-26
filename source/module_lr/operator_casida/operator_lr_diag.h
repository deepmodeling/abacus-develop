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
        OperatorLRDiag(const double* eig_ks, const Parallel_2D& pX_in, const int& nk_in, const int& nocc_in, const int& nvirt_in)
            : pX(pX_in), nk(nk_in), nocc(nocc_in), nvirt(nvirt_in)
        {   // calculate the difference of eigenvalues
            ModuleBase::TITLE("OperatorLRDiag", "OperatorLRDiag");
            const int nbands = nocc + nvirt;
            this->cal_type = hamilt::calculation_type::no;
            this->eig_ks_diff.create(nk, pX.get_local_size(), false);
            for (int ik = 0;ik < nk;++ik)
            {
                const int& istart = ik * nbands;
                for (int io = 0;io < pX.get_col_size();++io)    //nocc_local
                {
                    for (int iv = 0;iv < pX.get_row_size();++iv)    //nvirt_local
                    {
                        int io_g = pX.local2global_col(io);
                        int iv_g = pX.local2global_row(iv);
                        this->eig_ks_diff(ik, io * pX.get_row_size() + iv) = eig_ks[istart + nocc + iv_g] - eig_ks[istart + io_g];
                    }
                }
            }
        };
        void init(const int ik_in) override {};

        /// caution: put this operator at the head of the operator list,
        /// because vector_mul_vector_op directly assign to (rather than add on) psi_out.
        virtual void  act(const int nbands,
            const int nbasis,
            const int npol,
            const T* psi_in,
            T* hpsi,
            const int ngk_ik = 0,
            const bool is_first_node = false)const override
        {
            ModuleBase::TITLE("OperatorLRDiag", "act");
            hsolver::vector_mul_vector_op<T, Device>()(this->ctx,
                nk * pX.get_local_size(),   // local size of particle-hole basis
                hpsi,
                psi_in,
                this->eig_ks_diff.c);
        }
    private:
        const Parallel_2D& pX;
        ModuleBase::matrix eig_ks_diff;
        const int& nk;
        const int& nocc;
        const int& nvirt;
        Device* ctx = {};
    };
}