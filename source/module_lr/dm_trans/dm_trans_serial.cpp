#include "dm_trans.h"
#include "module_base/blas_connector.h"
#include "module_base/tool_title.h"
#include "module_base/global_function.h"
#include "module_lr/utils/lr_util.h"
namespace LR
{
    template<> std::vector<container::Tensor> cal_dm_trans_forloop_serial(
        const double* const X_istate,
        const psi::Psi<double>& c,
        const int& nocc,
        const int& nvirt,
        const bool renorm_k,
        const int nspin)
    {
        // cxc_out_test(X_istate, c);
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_trans_forloop");
        const int nks = c.get_nk();
        const int naos = c.get_nbasis();
        std::vector<container::Tensor> dm_trans(nks, container::Tensor(DAT::DT_DOUBLE, DEV::CpuDevice, { naos, naos }));
        for (auto& dm : dm_trans)ModuleBase::GlobalFunc::ZEROS(dm.data<double>(), naos * naos);
        // loop for AOs
        for (size_t isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int x_start = isk * nocc * nvirt;
            for (size_t mu = 0;mu < naos;++mu)
            {
                for (size_t nu = 0;nu < naos;++nu)
                {
                    // loop for ks states
                    for (size_t j = 0;j < nocc;++j)
                    {
                        for (size_t b = 0; b < nvirt;++b)
                            dm_trans[isk].data<double>()[mu * naos + nu] += c(j, mu) * X_istate[x_start + j * nvirt + b] * c(nocc + b, nu);
                    }
                }
            }
        }
        return dm_trans;
    }

    template<> std::vector<container::Tensor> cal_dm_trans_forloop_serial(
        const std::complex<double>* const X_istate,
        const psi::Psi<std::complex<double>>& c,
        const int& nocc,
        const int& nvirt,
        const bool renorm_k,
        const int nspin)
    {
        // cxc_out_test(X_istate, c);
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_trans_forloop");
        const int nks = c.get_nk();
        int naos = c.get_nbasis();
        std::vector<container::Tensor> dm_trans(nks, container::Tensor(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { naos, naos }));
        for (auto& dm : dm_trans)ModuleBase::GlobalFunc::ZEROS(dm.data<std::complex<double>>(), naos * naos);
        // loop for AOs
        for (size_t isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int x_start = isk * nocc * nvirt;
            for (size_t mu = 0;mu < naos;++mu)
            {
                for (size_t nu = 0;nu < naos;++nu)
                {
                    // loop for ks states
                    for (size_t j = 0;j < nocc;++j)
                    {
                        for (size_t b = 0; b < nvirt;++b)
                            dm_trans[isk].data<std::complex<double>>()[nu * naos + mu] +=
                            std::conj(c(j, mu)) * X_istate[x_start + j * nvirt + b] * c(nocc + b, nu) / static_cast<double>(nks / nspin);
                    }
                }
            }
        }
        return dm_trans;
    }


    template<> std::vector<container::Tensor> cal_dm_trans_blas(
        const double* const X_istate,
        const psi::Psi<double>& c,
        const int& nocc,
        const int& nvirt,
        const bool renorm_k,
        const int nspin)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_trans_blas");
        const int nks = c.get_nk();
        int naos = c.get_nbasis();
        std::vector<container::Tensor> dm_trans(nks, container::Tensor(DAT::DT_DOUBLE, DEV::CpuDevice, { naos, naos }));
        for (size_t isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int x_start = isk * nocc * nvirt;
            // 1. [X*C_occ^T]^T=C_occ*X^T
            char transa = 'N';
            char transb = 'T';
            const double alpha = 1.0;
            const double beta = 0.0;
            container::Tensor Xc(DAT::DT_DOUBLE, DEV::CpuDevice, { nvirt, naos });
            dgemm_(&transa, &transb, &naos, &nvirt, &nocc, &alpha,
                c.get_pointer(), &naos, X_istate + x_start, &nvirt,
                &beta, Xc.data<double>(), &naos);
            // 2. C_virt*[X*C_occ^T]
            dgemm_(&transa, &transb, &naos, &naos, &nvirt, &alpha,
                c.get_pointer(nocc), &naos, Xc.data<double>(), &naos, &beta,
                dm_trans[isk].data<double>(), &naos);
        }
        return dm_trans;
    }


    template<> std::vector<container::Tensor> cal_dm_trans_blas(
        const std::complex<double>* const X_istate,
        const psi::Psi<std::complex<double>>& c,
        const int& nocc,
        const int& nvirt,
        const bool renorm_k,
        const int nspin)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_trans_blas");
        const int nks = c.get_nk();
        int naos = c.get_nbasis();
        std::vector<container::Tensor> dm_trans(nks, container::Tensor(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { naos, naos }));
        for (size_t isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int x_start = isk * nocc * nvirt;

            char transa = 'N';
            char transb = 'C';
            std::complex<double> alpha(1.0, 0.0);
            const std::complex<double> beta(0.0, 0.0);

            // ============== C_virt * X * C_occ^\dagger=============
            // 1. [X*C_occ^\dagger]^\dagger=C_occ*X^\dagger
            // container::Tensor Xc(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { nvirt, naos });
            // zgemm_(&transa, &transb, &naos, &nvirt, &nocc, &alpha,
            //     c.get_pointer(), &naos, X_istate.get_pointer(), &nvirt,
            //     &beta, Xc.data<std::complex<double>>(), &naos);
            // // 2. C_virt*[X*C_occ^\dagger]
            // alpha = 1.0 / static_cast<double>(nks);
            // zgemm_(&transa, &transb, &naos, &naos, &nvirt, &alpha,
            //     c.get_pointer(nocc), &naos, Xc.data<std::complex<double>>(), &naos, &beta,
            //  dm_trans[isk].data<std::complex<double>>(), & naos);

            // ============== [C_virt * X * C_occ^\dagger]^T=============
            // ============== = [C_occ^* * X^T * C_virt^T]^T=============
            // 1. X*C_occ^\dagger
            container::Tensor Xc(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { naos, nvirt });
            zgemm_(&transa, &transb, &nvirt, &naos, &nocc, &alpha,
                X_istate + x_start, &nvirt, c.get_pointer(), &naos,
                &beta, Xc.data<std::complex<double>>(), &nvirt);
            // 2. [X*C_occ^\dagger]^TC_virt^T
            transa = transb = 'T';
            alpha.real(renorm_k ? 1.0 / static_cast<double>(nks / nspin) : 1.0);
            zgemm_(&transa, &transb, &naos, &naos, &nvirt, &alpha,
                Xc.data<std::complex<double>>(), &nvirt, c.get_pointer(nocc), &naos, &beta,
                dm_trans[isk].data<std::complex<double>>(), &naos);
        }
        return dm_trans;
    }

}
