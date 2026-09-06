#include "ao_to_mo.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/tool_title.h"
#include "source_lcao/module_lr/utils/lr_util.h"
namespace LR
{
    template<>
    void ao_to_mo_forloop_serial(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<double>& coeff,
        const int& nocc,
        const int& nvirt,
        double* mat_mo,
        const LR_Util::MO_TYPE type)
    {
        ModuleBase::TITLE("LR", "ao_to_mo_forloop_serial");
        const int nks = mat_ao.size();
        const int naos = coeff.get_nbasis();
        int nmo1, nmo2, imo1, imo2;
        LR_Util::set_dim(type, nocc, nvirt, nmo1, nmo2, imo1, imo2);

        ModuleBase::GlobalFunc::ZEROS(mat_mo, nks * nmo1 * nmo2);

        for (int isk = 0;isk < nks;++isk)
        {
            coeff.fix_k(isk);
            const int start = isk * nmo1 * nmo2;
            for (int p = 0;p < nmo1;++p)
            {
                for (int q = 0;q < nmo2;++q)
                {
                    for (int nu = 0;nu < naos;++nu)
                    {
                        for (int mu = 0;mu < naos;++mu)
                        {
                            mat_mo[start + p * nmo2 + q] += coeff(imo2 + q, mu) * mat_ao[isk].data<double>()[nu * naos + mu] * coeff(imo1 + p, nu);
                        }
                    }
                }
            }
        }
    }
    template<>
    void ao_to_mo_forloop_serial(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<std::complex<double>>& coeff,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* const mat_mo,
        const LR_Util::MO_TYPE type)
    {
        ModuleBase::TITLE("LR", "ao_to_mo_forloop_serial");
        const int nks = mat_ao.size();
        const int naos = coeff.get_nbasis();
        int nmo1, nmo2, imo1, imo2;
        LR_Util::set_dim(type, nocc, nvirt, nmo1, nmo2, imo1, imo2);

        ModuleBase::GlobalFunc::ZEROS(mat_mo, nks * nmo1 * nmo2);

        for (int isk = 0;isk < nks;++isk)
        {
            coeff.fix_k(isk);
            const int start = isk * nmo1 * nmo2;
            for (int p = 0;p < nmo1;++p)
            {
                for (int q = 0;q < nmo2;++q)
                {
                    for (int nu = 0;nu < naos;++nu)
                    {
                        for (int mu = 0;mu < naos;++mu)
                        {
                            mat_mo[start + p * nmo2 + q] += std::conj(coeff(imo2 + q, mu)) * mat_ao[isk].data<std::complex<double>>()[nu * naos + mu] * coeff(imo1 + p, nu);
                        }
                    }
                }
            }
        }
    }
    template<>
    void ao_to_mo_blas(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<double>& coeff,
        const int& nocc,
        const int& nvirt,
        double* mat_mo,
        const bool add_on,
        const LR_Util::MO_TYPE type)
    {
        ModuleBase::TITLE("LR", "ao_to_mo_blas");
        const int nks = mat_ao.size();
        const int naos = coeff.get_nbasis();        
        int nmo1, nmo2, imo1, imo2;
        LR_Util::set_dim(type, nocc, nvirt, nmo1, nmo2, imo1, imo2);

        for (int isk = 0;isk < nks;++isk)
        {
            coeff.fix_k(isk);
            const int start = isk * nmo1 * nmo2;

            // Vc[naos*nocc]
            container::Tensor Vc(DAT::DT_DOUBLE, DEV::CpuDevice, { nmo1, naos });// (Vc)^T
            Vc.zero();
            char transa = 'N';
            char transb = 'N';  //coeff is col major
            const double alpha = 1.0;
            const double beta = add_on ? 1.0 : 0.0;
            BlasConnector::gemm(transb, transa, nmo1, naos, naos, alpha,
                coeff.get_pointer(imo1), naos, mat_ao[isk].data<double>(), naos, beta,
                Vc.data<double>(), naos);

            transa = 'T';
            //mat_mo=coeff^TVc (nvirt major)
            BlasConnector::gemm(transb, transa, nmo1, nmo2, naos, alpha,
                Vc.data<double>(), naos, coeff.get_pointer(imo2), naos, beta,
                mat_mo + start, nmo2);
        }
    }
    template<>
    void ao_to_mo_blas(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<std::complex<double>>& coeff,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* const mat_mo,
        const bool add_on,
        const LR_Util::MO_TYPE type)
    {
        ModuleBase::TITLE("LR", "ao_to_mo_blas");
        const int nks = mat_ao.size();
        const int naos = coeff.get_nbasis();
        int nmo1, nmo2, imo1, imo2;
        LR_Util::set_dim(type, nocc, nvirt, nmo1, nmo2, imo1, imo2);

        for (int isk = 0;isk < nks;++isk)
        {
            coeff.fix_k(isk);
            const int start = isk * nmo1 * nmo2;

            // Vc[naos*nocc] (V is hermitian)
            container::Tensor Vc(DAT::DT_COMPLEX_DOUBLE, DEV::CpuDevice, { nmo1, naos });// (Vc)^T
            Vc.zero();
            char transa = 'N';
            char transb = 'N';  //coeff is col major
            const std::complex<double> alpha(1.0, 0.0);
            const std::complex<double> beta = add_on ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0);
            BlasConnector::gemm(transb, transa, nmo1, naos, naos, alpha,
                coeff.get_pointer(imo1), naos, mat_ao[isk].data<std::complex<double>>(), naos, beta,
                Vc.data<std::complex<double>>(), naos);

            transa = 'C';
            //mat_mo=coeff^\dagger Vc (nvirt major)
            BlasConnector::gemm(transb, transa, nmo1, nmo2, naos, alpha,
                Vc.data<std::complex<double>>(), naos, coeff.get_pointer(imo2), naos, beta,
                mat_mo + start, nmo2);
        }
    }
}