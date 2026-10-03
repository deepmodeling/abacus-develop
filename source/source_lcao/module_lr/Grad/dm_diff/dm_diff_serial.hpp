#pragma once
#include "source_base/module_container/ATen/core/tensor_types.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/tool_title.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include "dm_diff.h"
namespace LR
{
    template<typename T>
    inline void print_colfirst(const T* ptr, const std::string& name, const int& nrow, const int& ncol)
    {
        std::cout << name << std::endl;
        for (int i = 0;i < nrow;++i)
        {
            for (int j = 0;j < ncol;++j)
                std::cout << ptr[j * nrow + i] << " ";
            std::cout << std::endl;
        }
    }
    inline void CvX(
        const double* C,
        const double* X,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        double* CvX)
    {
        char transa = 'N';
        char transb = 'N';
        const double alpha = 1.0;
        const double beta = 0;
        dgemm_(&transa, &transb, &naos, &nocc, &nvirt,
            &alpha, C + nocc * naos, &naos,
            X, &nvirt,
            &beta, CvX, &naos);
    }
    inline void CvX(
        const std::complex<double>* C,
        const std::complex<double>* X,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* CvX)
    {
        char transa = 'N';
        char transb = 'N';
        const std::complex<double> alpha(1.0, 0.0);
        const std::complex<double> beta(0.0, 0.0);
        zgemm_(&transa, &transb, &naos, &nocc, &nvirt,
            &alpha, C + nocc * naos, &naos,
            X, &nvirt,
            &beta, CvX, &naos);
    }

    inline void CoXT(
        const double* C,
        const double* X,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        double* CoXT)
    {
        char transa = 'N';
        char transb = 'T';
        const double alpha = 1.0;
        const double beta = 0;
        dgemm_(&transa, &transb, &naos, &nvirt, &nocc,
            &alpha, C, &naos,
            X, &nvirt,
            &beta, CoXT, &naos);
    }
    inline void CoXT(
        const std::complex<double>* C,
        const std::complex<double>* X,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        std::complex<double>* CoXT)
    {
        char transa = 'N';
        char transb = 'T';
        const std::complex<double> alpha(1.0, 0.0);
        const std::complex<double> beta(0.0, 0.0);
        zgemm_(&transa, &transb, &naos, &nvirt, &nocc,
            &alpha, C, &naos,
            X, &nvirt,
            &beta, CoXT, &naos);
    }

    inline void AAT(
        const double* A,
        const int& nrow,
        const int& ncol,
        double* C,
        const bool add = false,
        const double& alpha = 1.0)
    {
        char transa = 'N';
        char transb = 'T';
        const double beta = add ? 1.0 : 0.0;
        dgemm_(&transa, &transb, &nrow, &nrow, &ncol,
            &alpha, A, &nrow,
            A, &nrow,
            &beta, C, &nrow);
    }
    inline void AAT(
        const std::complex<double>* A,
        const int& nrow,
        const int& ncol,
        std::complex<double>* C,
        const bool add = false,
        const std::complex<double>& alpha = std::complex<double>(1.0, 0.0))
    {
        char transa = 'N';
        char transb = 'T';
        const std::complex<double> beta(add ? 1.0 : 0.0, 0.0);
        // take conjugate of A
        std::vector<std::complex<double>> A_conj(nrow * ncol);
        for (int i = 0;i < A_conj.size();++i) A_conj[i] = std::conj(A[i]);
        zgemm_(&transa, &transb, &nrow, &nrow, &ncol,
            &alpha, A_conj.data(), &nrow,
            A, &nrow,
            &beta, C, &nrow);
    }

    //output: col first, consistent with blas
    // c: nao*nbands in para2d, nbands*nao in psi  (row-para and constructed: nao)
    // X: nvirt*nocc in para2d, nocc*nvirt in psi (row-para and constructed: nvirt)
    template<typename T>
    std::vector<ct::Tensor> cal_dm_diff_blas(
        const T* const X_istate,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const bool renorm_k,
        const int nspin)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_diff_blas");
        const int nks = c.get_nk();
        const int nk = nks / nspin;

        ct::Tensor cvx(ct::DataTypeToEnum<T>::value, DEV::CpuDevice, { nocc, naos });
        ct::Tensor coxt(ct::DataTypeToEnum<T>::value, DEV::CpuDevice, { nvirt, naos });
        std::vector<ct::Tensor> dm_diff(nks, ct::Tensor(ct::DataTypeToEnum<T>::value, DEV::CpuDevice, { naos, naos }));
        for (int iks = 0;iks < nks;++iks)
        {
            c.fix_k(iks);
            const int start = iks * nocc * nvirt;
            // 1. C_virt * X
            CvX(c.get_pointer(), X_istate + start, naos, nocc, nvirt, cvx.data<T>());
            // 2. C_occ * X^T
            CoXT(c.get_pointer(), X_istate + start, naos, nocc, nvirt, coxt.data<T>());
            // 3. cvx*cvx^T + coxt*coxt^T
            AAT(cvx.data<T>(), naos, nocc, dm_diff[iks].data<T>(), false, renorm_k ? (T)(1.0 / (double)nk) : (T)1.0);
            AAT(coxt.data<T>(), naos, nvirt, dm_diff[iks].data<T>(), true, renorm_k ? (T)(-1.0 / (double)nk) : (T)(-1.0));
        }
        return dm_diff;
    }

    inline double get_conj(const double& x) { return x; }
    inline std::complex<double> get_conj(const std::complex<double>& x) { return std::conj(x); }

    template<typename T>
    std::vector<ct::Tensor> cal_dm_diff_forloop(
        const T* const X_istate,
        const psi::Psi<T, base_device::DEVICE_CPU>& c,
        const int& naos,
        const int& nocc,
        const int& nvirt,
        const bool renorm_k,
        const int nspin)
    {
        ModuleBase::TITLE("hamilt_lrtd", "cal_dm_diff_forloop");
        const int nks = c.get_nk();
        const int nk = nks / nspin;

        std::vector<ct::Tensor> dm_diff(nks, ct::Tensor(ct::DataTypeToEnum<T>::value, DEV::CpuDevice, { naos, naos }));
        for (int iks = 0;iks < nks;++iks)
        {
            dm_diff[iks].zero();
            c.fix_k(iks);
            const int start = iks * nocc * nvirt;
            for (int nu = 0;nu < naos;++nu)//col
                for (int mu = 0;mu < naos;++mu)//row
                {
                    for (int i = 0;i < nocc;++i)
                        for (int a = 0;a < nvirt;++a)
                        {
                            for (int b = 0;b < nvirt;++b)
                                dm_diff[iks].data<T>()[nu * naos + mu]
                                += get_conj(c.get_pointer()[(nocc + a) * naos + mu] * X_istate[start + i * nvirt + a])
                                * c.get_pointer()[(nocc + b) * naos + nu] * X_istate[start + i * nvirt + b];
                            for (int j = 0;j < nocc;++j)
                                dm_diff[iks].data<T>()[nu * naos + mu]
                                -= get_conj(c.get_pointer()[i * naos + mu] * X_istate[start + i * nvirt + a])
                                * c.get_pointer()[j * naos + nu] * X_istate[start + j * nvirt + a];
                        }
                    if (renorm_k)
                        dm_diff[iks].data<T>()[nu * naos + mu] /= (double)nk;
                }
        }
        return dm_diff;
    }
}