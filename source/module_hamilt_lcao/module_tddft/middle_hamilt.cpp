#include "middle_hamilt.h"

#include "module_base/lapack_connector.h"
#include "module_base/module_container/ATen/kernels/blas.h"
#include "module_base/scalapack_connector.h"

#include <complex>
#include <iostream>

namespace module_tddft
{
#ifdef __MPI

void half_Hmatrix(const Parallel_Orbitals* pv,
                  const int nband,
                  const int nlocal,
                  std::complex<double>* Htmp,
                  std::complex<double>* Stmp,
                  const std::complex<double>* H_laststep,
                  const std::complex<double>* S_laststep,
                  const int print_matrix)
{
    if (print_matrix)
    {
        GlobalV::ofs_running << std::setprecision(10);
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H(t+dt) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << Htmp[i * pv->ncol + j].real() << "+" << Htmp[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H(t):" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << H_laststep[i * pv->ncol + j].real() << "+"
                                     << H_laststep[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    std::complex<double> alpha = {0.5, 0.0};
    std::complex<double> beta = {0.5, 0.0};
    ScalapackConnector::geadd('N', nlocal, nlocal, alpha, H_laststep, 1, 1, pv->desc, beta, Htmp, 1, 1, pv->desc);
    ScalapackConnector::geadd('N', nlocal, nlocal, alpha, S_laststep, 1, 1, pv->desc, beta, Stmp, 1, 1, pv->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H (t+dt/2) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << Htmp[i * pv->ncol + j].real() << "+" << Htmp[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }
}

void half_Hmatrix_tensor(const Parallel_Orbitals* pv,
                         const int nband,
                         const int nlocal,
                         container::Tensor& Htmp,
                         container::Tensor& Stmp,
                         const container::Tensor& H_laststep,
                         const container::Tensor& S_laststep,
                         const int print_matrix)
{
    if (print_matrix)
    {
        GlobalV::ofs_running << std::setprecision(10);
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H(t+dt) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << Htmp.data<std::complex<double>>()[i * pv->ncol + j].real() << "+"
                                     << Htmp.data<std::complex<double>>()[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H(t):" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << H_laststep.data<std::complex<double>>()[i * pv->ncol + j].real() << "+"
                                     << H_laststep.data<std::complex<double>>()[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    std::complex<double> alpha = {0.5, 0.0};
    std::complex<double> beta = {0.5, 0.0};

    // Perform the operation Htmp = alpha * H_laststep + beta * Htmp
    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              H_laststep.data<std::complex<double>>(),
                              1,
                              1,
                              pv->desc,
                              beta,
                              Htmp.data<std::complex<double>>(),
                              1,
                              1,
                              pv->desc);

    // Perform the operation Stmp = alpha * S_laststep + beta * Stmp
    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              S_laststep.data<std::complex<double>>(),
                              1,
                              1,
                              pv->desc,
                              beta,
                              Stmp.data<std::complex<double>>(),
                              1,
                              1,
                              pv->desc);

    if (print_matrix)
    {
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H (t+dt/2) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << Htmp.data<std::complex<double>>()[i * pv->ncol + j].real() << "+"
                                     << Htmp.data<std::complex<double>>()[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }
}

void half_Hmatrix_tensor_lapack(const Parallel_Orbitals* pv,
                                const int nband,
                                const int nlocal,
                                container::Tensor& Htmp,
                                container::Tensor& Stmp,
                                const container::Tensor& H_laststep,
                                const container::Tensor& S_laststep,
                                const int print_matrix)
{
    if (print_matrix)
    {
        GlobalV::ofs_running << std::setprecision(10);
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H(t+dt) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << Htmp.data<std::complex<double>>()[i * pv->ncol + j].real() << "+"
                                     << Htmp.data<std::complex<double>>()[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H(t):" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << H_laststep.data<std::complex<double>>()[i * pv->ncol + j].real() << "+"
                                     << H_laststep.data<std::complex<double>>()[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }

    std::complex<double> alpha = {0.5, 0.0};
    std::complex<double> beta = {0.5, 0.0};

    // Perform the operation Htmp = alpha * H_laststep + beta * Htmp
    // Scale Htmp by beta
    container::kernels::blas_scal<std::complex<double>, container::DEVICE_CPU>()(nlocal * nlocal,
                                                                                 &beta,
                                                                                 Htmp.data<std::complex<double>>(),
                                                                                 1);
    // Htmp = alpha * H_laststep + Htmp
    container::kernels::blas_axpy<std::complex<double>, container::DEVICE_CPU>()(
        nlocal * nlocal,
        &alpha,
        H_laststep.data<std::complex<double>>(),
        1,
        Htmp.data<std::complex<double>>(),
        1);

    // Perform the operation Stmp = alpha * S_laststep + beta * Stmp
    // Scale Stmp by beta
    container::kernels::blas_scal<std::complex<double>, container::DEVICE_CPU>()(nlocal * nlocal,
                                                                                 &beta,
                                                                                 Stmp.data<std::complex<double>>(),
                                                                                 1);
    // Stmp = alpha * S_laststep + Stmp
    container::kernels::blas_axpy<std::complex<double>, container::DEVICE_CPU>()(
        nlocal * nlocal,
        &alpha,
        S_laststep.data<std::complex<double>>(),
        1,
        Stmp.data<std::complex<double>>(),
        1);

    if (print_matrix)
    {
        GlobalV::ofs_running << std::endl;
        GlobalV::ofs_running << " H (t+dt/2) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            for (int j = 0; j < pv->ncol; j++)
            {
                GlobalV::ofs_running << Htmp.data<std::complex<double>>()[i * pv->ncol + j].real() << "+"
                                     << Htmp.data<std::complex<double>>()[i * pv->ncol + j].imag() << "i ";
            }
            GlobalV::ofs_running << std::endl;
        }
        GlobalV::ofs_running << std::endl;
    }
}

#endif
} // namespace module_tddft