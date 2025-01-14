#include "AX.h"
#include "module_base/blas_connector.h"
#include "module_base/tool_title.h"
#include "module_lr/utils/lr_util.h"
namespace LR
{
    void cal_AX_forloop_serial(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<double>& c,
        const int& nocc,
        const int& nvirt,
        double* AX_istate)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_AX_forloop");
        const int nks = V_istate.size();
        int naos = c.get_nbasis();
        ModuleBase::GlobalFunc::ZEROS(AX_istate, nks * nocc * nvirt);

        for (int isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int ax_start = isk * nocc * nvirt;
            for (int i = 0;i < nocc;++i)
            {
                for (int a = 0;a < nvirt;++a)
                {
                    for (int nu = 0;nu < naos;++nu)
                    {
                        for (int mu = 0;mu < naos;++mu)
                        {
                            AX_istate[ax_start + i * nvirt + a] += c(nocc + a, mu) * V_istate[isk].data<double>()[nu * naos + mu] * c(i, nu);
                        }
                    }
                }
            }
        }
    }
    void cal_AX_forloop_serial(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<std::complex<double>>& c,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* const AX_istate)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_AX_forloop");
        const int nks = V_istate.size();
        int naos = c.get_nbasis();
        ModuleBase::GlobalFunc::ZEROS(AX_istate, nks * nocc * nvirt);

        for (int isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int ax_start = isk * nocc * nvirt;
            for (int i = 0;i < nocc;++i)
            {
                for (int a = 0;a < nvirt;++a)
                {
                    for (int nu = 0;nu < naos;++nu)
                    {
                        for (int mu = 0;mu < naos;++mu)
                        {
                            AX_istate[ax_start + i * nvirt + a] += std::conj(c(nocc + a, mu)) * V_istate[isk].data<std::complex<double>>()[nu * naos + mu] * c(i, nu);
                        }
                    }
                }
            }
        }
    }

    void cal_AX_blas(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<double>& c,
        const int& nocc,
        const int& nvirt,
        double* AX_istate,
        const bool add_on)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_AX_blas");
        const int nks = V_istate.size();
        int naos = c.get_nbasis();

        for (int isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int ax_start = isk * nocc * nvirt;

            // Vc[naos*nocc]
            container::Tensor Vc(DAT::DT_DOUBLE, DEV::CpuDevice, { nocc, naos });// (Vc)^T
            Vc.zero();
            char transa = 'N';
            char transb = 'N';  //c is col major
            const double alpha = 1.0;
            const double beta = add_on ? 1.0 : 0.0;
            dgemm_(&transa, &transb, &naos, &nocc, &naos, &alpha,
                V_istate[isk].data<double>(), &naos, c.get_pointer(), &naos, &beta,
                Vc.data<double>(), &naos);

            transa = 'T';
            //AX_istate=c^TVc (nvirt major)
            dgemm_(&transa, &transb, &nvirt, &nocc, &naos, &alpha,
                c.get_pointer(nocc), &naos, Vc.data<double>(), &naos, &beta,
                AX_istate + ax_start, &nvirt);
        }
    }
    void cal_AX_blas(
        const std::vector<container::Tensor>& V_istate,
        const psi::Psi<std::complex<double>>& c,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* const AX_istate,
        const bool add_on)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_AX_blas");
        const int nks = V_istate.size();
        int naos = c.get_nbasis();

        for (int isk = 0;isk < nks;++isk)
        {
            c.fix_k(isk);
            const int ax_start = isk * nocc * nvirt;

            // Vc[naos*nocc] (V is hermitian)
            container::Tensor Vc(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { nocc, naos });// (Vc)^T
            Vc.zero();
            char transa = 'N';
            char transb = 'N';  //c is col major
            const std::complex<double> alpha(1.0, 0.0);
            const std::complex<double> beta = add_on ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0);
            zgemm_(&transa, &transb, &naos, &nocc, &naos, &alpha,
                V_istate[isk].data<std::complex<double>>(), &naos, c.get_pointer(), &naos, &beta,
                Vc.data<std::complex<double>>(), &naos);

            transa = 'C';
            //AX_istate=c^\dagger Vc (nvirt major)
            zgemm_(&transa, &transb, &nvirt, &nocc, &naos, &alpha,
                c.get_pointer(nocc), &naos, Vc.data<std::complex<double>>(), &naos, &beta,
                AX_istate + ax_start, &nvirt);
        }
    }
}