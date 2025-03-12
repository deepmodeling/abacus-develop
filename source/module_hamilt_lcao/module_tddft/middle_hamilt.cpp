#include "middle_hamilt.h"

#include "module_base/lapack_connector.h"
#include "module_base/module_container/ATen/kernels/blas.h"
#include "module_base/module_device/memory_op.h" // memory operations
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
                  std::ofstream& ofs_running,
                  const int print_matrix)
{
    if (print_matrix)
    {
        ofs_running << std::setprecision(10);
        ofs_running << std::endl;
        ofs_running << " H(t+dt) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                ofs_running << Htmp[in + j].real() << "+" << Htmp[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " H(t):" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                ofs_running << H_laststep[in + j].real() << "+" << H_laststep[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    std::complex<double> alpha = {0.5, 0.0};
    std::complex<double> beta = {0.5, 0.0};
    ScalapackConnector::geadd('N', nlocal, nlocal, alpha, H_laststep, 1, 1, pv->desc, beta, Htmp, 1, 1, pv->desc);
    ScalapackConnector::geadd('N', nlocal, nlocal, alpha, S_laststep, 1, 1, pv->desc, beta, Stmp, 1, 1, pv->desc);

    if (print_matrix)
    {
        ofs_running << std::endl;
        ofs_running << " H (t+dt/2) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                ofs_running << Htmp[in + j].real() << "+" << Htmp[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }
}

void half_Hmatrix_tensor(const Parallel_Orbitals* pv,
                         const int nband,
                         const int nlocal,
                         ct::Tensor& Htmp,
                         ct::Tensor& Stmp,
                         const ct::Tensor& H_laststep,
                         const ct::Tensor& S_laststep,
                         std::ofstream& ofs_running,
                         const int print_matrix)
{
    if (print_matrix)
    {
        ofs_running << std::setprecision(10);
        ofs_running << std::endl;
        ofs_running << " H(t+dt) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                ofs_running << Htmp.data<std::complex<double>>()[in + j].real() << "+"
                            << Htmp.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " H(t):" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                ofs_running << H_laststep.data<std::complex<double>>()[in + j].real() << "+"
                            << H_laststep.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
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
        ofs_running << std::endl;
        ofs_running << " H (t+dt/2) :" << std::endl;
        for (int i = 0; i < pv->nrow; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                ofs_running << Htmp.data<std::complex<double>>()[in + j].real() << "+"
                            << Htmp.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }
}

template <typename Device>
void half_Hmatrix_tensor_lapack(const Parallel_Orbitals* pv,
                                const int nband,
                                const int nlocal,
                                ct::Tensor& Htmp,
                                ct::Tensor& Stmp,
                                const ct::Tensor& H_laststep,
                                const ct::Tensor& S_laststep,
                                std::ofstream& ofs_running,
                                const int print_matrix)
{
    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    if (print_matrix)
    {
        ct::Tensor Htmp_cpu = Htmp.to_device<ct::DEVICE_CPU>();
        ct::Tensor H_laststep_cpu = H_laststep.to_device<ct::DEVICE_CPU>();

        ofs_running << std::setprecision(10);
        ofs_running << std::endl;
        ofs_running << " H(t+dt) :" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                ofs_running << Htmp_cpu.data<std::complex<double>>()[in + j].real() << "+"
                            << Htmp_cpu.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " H(t):" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                ofs_running << H_laststep_cpu.data<std::complex<double>>()[in + j].real() << "+"
                            << H_laststep_cpu.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    std::complex<double> one_half = {0.5, 0.0};

    // Perform the operation Htmp = one_half * H_laststep + one_half * Htmp
    // Scale Htmp by one_half
    ct::kernels::blas_scal<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                              &one_half,
                                                              Htmp.data<std::complex<double>>(),
                                                              1);
    // Htmp = one_half * H_laststep + Htmp
    ct::kernels::blas_axpy<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                              &one_half,
                                                              H_laststep.data<std::complex<double>>(),
                                                              1,
                                                              Htmp.data<std::complex<double>>(),
                                                              1);

    // Perform the operation Stmp = one_half * S_laststep + one_half * Stmp
    // Scale Stmp by one_half
    ct::kernels::blas_scal<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                              &one_half,
                                                              Stmp.data<std::complex<double>>(),
                                                              1);
    // Stmp = one_half * S_laststep + Stmp
    ct::kernels::blas_axpy<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                              &one_half,
                                                              S_laststep.data<std::complex<double>>(),
                                                              1,
                                                              Stmp.data<std::complex<double>>(),
                                                              1);

    if (print_matrix)
    {
        ct::Tensor Htmp_cpu = Htmp.to_device<ct::DEVICE_CPU>();

        ofs_running << std::endl;
        ofs_running << " H (t+dt/2) :" << std::endl;
        for (int i = 0; i < nlocal; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                ofs_running << Htmp_cpu.data<std::complex<double>>()[in + j].real() << "+"
                            << Htmp_cpu.data<std::complex<double>>()[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }
}

// Explicit instantiation of template functions
template void half_Hmatrix_tensor_lapack<base_device::DEVICE_CPU>(const Parallel_Orbitals* pv,
                                                                  const int nband,
                                                                  const int nlocal,
                                                                  ct::Tensor& Htmp,
                                                                  ct::Tensor& Stmp,
                                                                  const ct::Tensor& H_laststep,
                                                                  const ct::Tensor& S_laststep,
                                                                  std::ofstream& ofs_running,
                                                                  const int print_matrix);
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void half_Hmatrix_tensor_lapack<base_device::DEVICE_GPU>(const Parallel_Orbitals* pv,
                                                                  const int nband,
                                                                  const int nlocal,
                                                                  ct::Tensor& Htmp,
                                                                  ct::Tensor& Stmp,
                                                                  const ct::Tensor& H_laststep,
                                                                  const ct::Tensor& S_laststep,
                                                                  std::ofstream& ofs_running,
                                                                  const int print_matrix);
#endif // __CUDA
#endif // __MPI
} // namespace module_tddft