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
        const double factor,
        const MO_TYPE type)
    {
        // cxc_out_test(X_istate, c);
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_trans_forloop");
        const int nks = c.get_nk();
        const int imo1 = type == MO_TYPE::VV ? nocc : 0;
        const int imo2 = type == MO_TYPE::OO ? 0 : nocc;
        const int nmo1 = type == MO_TYPE::VV ? nvirt : nocc;
        const int nmo2 = type == MO_TYPE::OO ? nocc : nvirt;

        const int naos = c.get_nbasis();
        std::vector<container::Tensor> dm_trans(nks, container::Tensor(DAT::DT_DOUBLE, DEV::CpuDevice, { naos, naos }));
        for (auto& dm : dm_trans)ModuleBase::GlobalFunc::ZEROS(dm.data<double>(), naos * naos);
        // loop for AOs
        for (size_t isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int x_start = isk * nmo1 * nmo2;
            for (size_t mu = 0;mu < naos;++mu)
            {
                for (size_t nu = 0;nu < naos;++nu)
                {
                    // loop for ks states
                    for (size_t p = 0;p < nmo1;++p)
                    {
                        for (size_t q = 0; q < nmo2;++q)
                            dm_trans[isk].data<double>()[mu * naos + nu] += c(imo1 + p, mu) * X_istate[x_start + p * nmo2 + q] * c(imo2 + q, nu) * factor;
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
        const std::complex<double> factor,
        const MO_TYPE type)
    {
        // cxc_out_test(X_istate, c);
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_trans_forloop");
        const int nks = c.get_nk();
        const int naos = c.get_nbasis();
        const int imo1 = type == MO_TYPE::VV ? nocc : 0;
        const int imo2 = type == MO_TYPE::OO ? 0 : nocc;
        const int nmo1 = type == MO_TYPE::VV ? nvirt : nocc;
        const int nmo2 = type == MO_TYPE::OO ? nocc : nvirt;
        std::vector<container::Tensor> dm_trans(nks, container::Tensor(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { naos, naos }));
        for (auto& dm : dm_trans)ModuleBase::GlobalFunc::ZEROS(dm.data<std::complex<double>>(), naos * naos);
        // loop for AOs
        for (size_t isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int x_start = isk * nmo1 * nmo2;
            for (size_t mu = 0;mu < naos;++mu)
            {
                for (size_t nu = 0;nu < naos;++nu)
                {
                    // loop for ks states
                    for (size_t p = 0;p < nmo1;++p)
                    {
                        for (size_t q = 0; q < nmo2;++q)
                            dm_trans[isk].data<std::complex<double>>()[nu * naos + mu] +=
                            std::conj(c(imo1 + p, mu)) * X_istate[x_start + p * nmo2 + q] * c(imo2 + q, nu) * factor;
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
        const double factor,
        const MO_TYPE type)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_trans_blas");
        const int nks = c.get_nk();
        const int naos = c.get_nbasis();
        const int imo1 = type == MO_TYPE::VV ? nocc : 0;
        const int imo2 = type == MO_TYPE::OO ? 0 : nocc;
        const int nmo1 = type == MO_TYPE::VV ? nvirt : nocc;
        const int nmo2 = type == MO_TYPE::OO ? nocc : nvirt;
        std::vector<container::Tensor> dm_trans(nks, container::Tensor(DAT::DT_DOUBLE, DEV::CpuDevice, { naos, naos }));
        for (size_t isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int x_start = isk * nmo1 * nmo2;
            // 1. [X*C_occ^T]^T=C_occ*X^T
            char transa = 'N';
            char transb = 'T';
            const double alpha = 1.0;
            const double beta = 0.0;
            container::Tensor Xc(DAT::DT_DOUBLE, DEV::CpuDevice, { nmo2, naos });
            dgemm_(&transa, &transb, &naos, &nmo2, &nmo1, &alpha,
                c.get_pointer(imo1), &naos, X_istate + x_start, &nmo2,
                &beta, Xc.data<double>(), &naos);
            // 2. C_virt*[X*C_occ^T]
            dgemm_(&transa, &transb, &naos, &naos, &nmo2, &factor,
                c.get_pointer(imo2), &naos, Xc.data<double>(), &naos, &beta,
                dm_trans[isk].data<double>(), &naos);
        }
        return dm_trans;
    }


    template<> std::vector<container::Tensor> cal_dm_trans_blas(
        const std::complex<double>* const X_istate,
        const psi::Psi<std::complex<double>>& c,
        const int& nocc,
        const int& nvirt,
        const std::complex<double> factor,
        const MO_TYPE type)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_trans_blas");
        const int nks = c.get_nk();
        const int naos = c.get_nbasis();
        const int imo1 = type == MO_TYPE::VV ? nocc : 0;
        const int imo2 = type == MO_TYPE::OO ? 0 : nocc;
        const int nmo1 = type == MO_TYPE::VV ? nvirt : nocc;
        const int nmo2 = type == MO_TYPE::OO ? nocc : nvirt;
        std::vector<container::Tensor> dm_trans(nks, container::Tensor(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { naos, naos }));
        for (size_t isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int x_start = isk * nmo1 * nmo2;

            char transa = 'N';
            char transb = 'C';
            const std::complex<double> alpha(1.0, 0.0);
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
            container::Tensor Xc(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { naos, nmo2 });
            zgemm_(&transa, &transb, &nmo2, &naos, &nmo1, &alpha,
                X_istate + x_start, &nmo2, c.get_pointer(imo1), &naos,
                &beta, Xc.data<std::complex<double>>(), &nmo2);
            // 2. [X*C_occ^\dagger]^TC_virt^T
            transa = transb = 'T';
            zgemm_(&transa, &transb, &naos, &naos, &nmo2, &factor,
                Xc.data<std::complex<double>>(), &nmo2, c.get_pointer(imo2), &naos, &beta,
                dm_trans[isk].data<std::complex<double>>(), &naos);
        }
        return dm_trans;
    }

}
