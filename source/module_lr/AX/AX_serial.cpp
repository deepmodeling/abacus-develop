#include "AX.h"
#include "module_base/blas_connector.h"
#include "module_base/tool_title.h"
#include "module_lr/utils/lr_util.h"
namespace LR
{
    void cal_AX_forloop_serial(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<double>& coeff,
        const int& nocc,
        const int& nvirt,
        double* mat_mo)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_AX_forloop");
        const int nks = mat_ao.size();
        int naos = coeff.get_nbasis();
        ModuleBase::GlobalFunc::ZEROS(mat_mo, nks * nocc * nvirt);

        for (int isk = 0;isk < nks;++isk)
        {
            coeff.fix_k(isk);
            const int ax_start = isk * nocc * nvirt;
            for (int i = 0;i < nocc;++i)
            {
                for (int a = 0;a < nvirt;++a)
                {
                    for (int nu = 0;nu < naos;++nu)
                    {
                        for (int mu = 0;mu < naos;++mu)
                        {
                            mat_mo[ax_start + i * nvirt + a] += coeff(nocc + a, mu) * mat_ao[isk].data<double>()[nu * naos + mu] * coeff(i, nu);
                        }
                    }
                }
            }
        }
    }
    void cal_AX_forloop_serial(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<std::complex<double>>& coeff,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* const mat_mo)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_AX_forloop");
        const int nks = mat_ao.size();
        int naos = coeff.get_nbasis();
        ModuleBase::GlobalFunc::ZEROS(mat_mo, nks * nocc * nvirt);

        for (int isk = 0;isk < nks;++isk)
        {
            coeff.fix_k(isk);
            const int ax_start = isk * nocc * nvirt;
            for (int i = 0;i < nocc;++i)
            {
                for (int a = 0;a < nvirt;++a)
                {
                    for (int nu = 0;nu < naos;++nu)
                    {
                        for (int mu = 0;mu < naos;++mu)
                        {
                            mat_mo[ax_start + i * nvirt + a] += std::conj(coeff(nocc + a, mu)) * mat_ao[isk].data<std::complex<double>>()[nu * naos + mu] * coeff(i, nu);
                        }
                    }
                }
            }
        }
    }

    void cal_AX_blas(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<double>& coeff,
        const int& nocc,
        const int& nvirt,
        double* mat_mo,
        const bool add_on)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_AX_blas");
        const int nks = mat_ao.size();
        int naos = coeff.get_nbasis();

        for (int isk = 0;isk < nks;++isk)
        {
            coeff.fix_k(isk);
            const int ax_start = isk * nocc * nvirt;

            // Vc[naos*nocc]
            container::Tensor Vc(DAT::DT_DOUBLE, DEV::CpuDevice, { nocc, naos });// (Vc)^T
            Vc.zero();
            char transa = 'N';
            char transb = 'N';  //coeff is col major
            const double alpha = 1.0;
            const double beta = add_on ? 1.0 : 0.0;
            dgemm_(&transa, &transb, &naos, &nocc, &naos, &alpha,
                mat_ao[isk].data<double>(), &naos, coeff.get_pointer(), &naos, &beta,
                Vc.data<double>(), &naos);

            transa = 'T';
            //mat_mo=coeff^TVc (nvirt major)
            dgemm_(&transa, &transb, &nvirt, &nocc, &naos, &alpha,
                coeff.get_pointer(nocc), &naos, Vc.data<double>(), &naos, &beta,
                mat_mo + ax_start, &nvirt);
        }
    }
    void cal_AX_blas(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<std::complex<double>>& coeff,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* const mat_mo,
        const bool add_on)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_AX_blas");
        const int nks = mat_ao.size();
        int naos = coeff.get_nbasis();

        for (int isk = 0;isk < nks;++isk)
        {
            coeff.fix_k(isk);
            const int ax_start = isk * nocc * nvirt;

            // Vc[naos*nocc] (V is hermitian)
            container::Tensor Vc(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { nocc, naos });// (Vc)^T
            Vc.zero();
            char transa = 'N';
            char transb = 'N';  //coeff is col major
            const std::complex<double> alpha(1.0, 0.0);
            const std::complex<double> beta = add_on ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0);
            zgemm_(&transa, &transb, &naos, &nocc, &naos, &alpha,
                mat_ao[isk].data<std::complex<double>>(), &naos, coeff.get_pointer(), &naos, &beta,
                Vc.data<std::complex<double>>(), &naos);

            transa = 'C';
            //mat_mo=coeff^\dagger Vc (nvirt major)
            zgemm_(&transa, &transb, &nvirt, &nocc, &naos, &alpha,
                coeff.get_pointer(nocc), &naos, Vc.data<std::complex<double>>(), &naos, &beta,
                mat_mo + ax_start, &nvirt);
        }
    }
}