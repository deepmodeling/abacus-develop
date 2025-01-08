#include "ao_to_mo.h"
#include "module_base/blas_connector.h"
#include "module_base/tool_title.h"
#include "module_lr/utils/lr_util.h"
namespace LR
{
    template<>
    void ao_to_mo_forloop_serial(
        const std::vector<container::Tensor>& mat_ao,
        const psi::Psi<double>& coeff,
        const int& nocc,
        const int& nvirt,
        double* mat_mo,
        MO_TYPE type)
    {
        ModuleBase::TITLE("hamilt_lrtd", "ao_to_mo_forloop_serial");
        const int nks = mat_ao.size();
        const int naos = coeff.get_nbasis();
        const int nmo1 = type == MO_TYPE::VV ? nvirt : nocc;
        const int nmo2 = type == MO_TYPE::OO ? nocc : nvirt;
        const int imo1 = type == MO_TYPE::VV ? nocc : 0;
        const int imo2 = type == MO_TYPE::OO ? 0 : nocc;

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
        MO_TYPE type)
    {
        ModuleBase::TITLE("hamilt_lrtd", "ao_to_mo_forloop_serial");
        const int nks = mat_ao.size();
        const int naos = coeff.get_nbasis();
        const int nmo1 = type == MO_TYPE::VV ? nvirt : nocc;
        const int nmo2 = type == MO_TYPE::OO ? nocc : nvirt;
        const int imo1 = type == MO_TYPE::VV ? nocc : 0;
        const int imo2 = type == MO_TYPE::OO ? 0 : nocc;

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
        MO_TYPE type)
    {
        ModuleBase::TITLE("hamilt_lrtd", "ao_to_mo_blas");
        const int nks = mat_ao.size();
        const int naos = coeff.get_nbasis();
        const int nmo1 = type == MO_TYPE::VV ? nvirt : nocc;
        const int nmo2 = type == MO_TYPE::OO ? nocc : nvirt;
        const int imo1 = type == MO_TYPE::VV ? nocc : 0;
        const int imo2 = type == MO_TYPE::OO ? 0 : nocc;

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
            dgemm_(&transa, &transb, &naos, &nmo1, &naos, &alpha,
                mat_ao[isk].data<double>(), &naos, coeff.get_pointer(imo1), &naos, &beta,
                Vc.data<double>(), &naos);

            transa = 'T';
            //mat_mo=coeff^TVc (nvirt major)
            dgemm_(&transa, &transb, &nmo2, &nmo1, &naos, &alpha,
                coeff.get_pointer(imo2), &naos, Vc.data<double>(), &naos, &beta,
                mat_mo + start, &nmo2);
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
        MO_TYPE type)
    {
        ModuleBase::TITLE("hamilt_lrtd", "ao_to_mo_blas");
        const int nks = mat_ao.size();
        const int naos = coeff.get_nbasis();
        const int nmo1 = type == MO_TYPE::VV ? nvirt : nocc;
        const int nmo2 = type == MO_TYPE::OO ? nocc : nvirt;
        const int imo1 = type == MO_TYPE::VV ? nocc : 0;
        const int imo2 = type == MO_TYPE::OO ? 0 : nocc;

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
            zgemm_(&transa, &transb, &naos, &nmo1, &naos, &alpha,
                mat_ao[isk].data<std::complex<double>>(), &naos, coeff.get_pointer(imo1), &naos, &beta,
                Vc.data<std::complex<double>>(), &naos);

            transa = 'C';
            //mat_mo=coeff^\dagger Vc (nvirt major)
            zgemm_(&transa, &transb, &nmo2, &nmo1, &naos, &alpha,
                coeff.get_pointer(imo2), &naos, Vc.data<std::complex<double>>(), &naos, &beta,
                mat_mo + start, &nmo2);
        }
    }
}