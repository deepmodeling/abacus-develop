#include "upsi.h"

#include "module_base/lapack_connector.h"
#include "module_base/module_container/ATen/kernels/blas.h"
#include "module_base/scalapack_connector.h"

#include <complex>
#include <iostream>

namespace module_tddft
{
#ifdef __MPI
void upsi(const Parallel_Orbitals* pv,
          const int nband,
          const int nlocal,
          const std::complex<double>* U_operator,
          const std::complex<double>* psi_k_laststep,
          std::complex<double>* psi_k,
          std::ofstream& ofs_running,
          const int print_matrix)
{
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nlocal,
                             1.0,
                             U_operator,
                             1,
                             1,
                             pv->desc,
                             psi_k_laststep,
                             1,
                             1,
                             pv->desc_wfc,
                             0.0,
                             psi_k,
                             1,
                             1,
                             pv->desc_wfc);

    if (print_matrix)
    {
        ofs_running << std::endl;
        ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa = psi_k[in + j].real();
                double bb = psi_k[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " psi_k_laststep:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa = psi_k_laststep[in + j].real();
                double bb = psi_k_laststep[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }
}

void upsi_tensor(const Parallel_Orbitals* pv,
                 const int nband,
                 const int nlocal,
                 const ct::Tensor& U_operator,
                 const ct::Tensor& psi_k_laststep,
                 ct::Tensor& psi_k,
                 std::ofstream& ofs_running,
                 const int print_matrix)
{
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nlocal,
                             1.0,
                             U_operator.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc,
                             psi_k_laststep.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_wfc,
                             0.0,
                             psi_k.data<std::complex<double>>(),
                             1,
                             1,
                             pv->desc_wfc);

    if (print_matrix)
    {
        ofs_running << std::endl;
        ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa = psi_k.data<std::complex<double>>()[in + j].real();
                double bb = psi_k.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " psi_k_laststep:" << std::endl;
        for (int i = 0; i < pv->ncol_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol; j++)
            {
                double aa = psi_k_laststep.data<std::complex<double>>()[in + j].real();
                double bb = psi_k_laststep.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }
}

template <typename Device>
void upsi_tensor_lapack(const Parallel_Orbitals* pv,
                        const int nband,
                        const int nlocal,
                        const ct::Tensor& U_operator,
                        const ct::Tensor& psi_k_laststep,
                        ct::Tensor& psi_k,
                        std::ofstream& ofs_running,
                        const int print_matrix)
{
    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    // Perform the matrix multiplication: psi_k = U_operator * psi_k_laststep
    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta = {0.0, 0.0};

    ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                              'N',
                                                              nlocal,
                                                              nband,
                                                              nlocal,
                                                              &alpha,
                                                              U_operator.data<std::complex<double>>(),
                                                              nlocal,
                                                              psi_k_laststep.data<std::complex<double>>(),
                                                              nlocal,
                                                              &beta,
                                                              psi_k.data<std::complex<double>>(),
                                                              nlocal);

    if (print_matrix)
    {
        ct::Tensor psi_k_cpu = psi_k.to_device<ct::DEVICE_CPU>();
        ct::Tensor psi_k_laststep_cpu = psi_k_laststep.to_device<ct::DEVICE_CPU>();

        ofs_running << std::endl;
        ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < nband; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                double aa = psi_k_cpu.data<std::complex<double>>()[in + j].real();
                double bb = psi_k_cpu.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " psi_k_laststep:" << std::endl;
        for (int i = 0; i < nband; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nlocal; j++)
            {
                double aa = psi_k_laststep_cpu.data<std::complex<double>>()[in + j].real();
                double bb = psi_k_laststep_cpu.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < 1e-8)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < 1e-8)
                {
                    bb = 0.0;
                }
                ofs_running << aa << "+" << bb << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }
}

// Explicit instantiation of template functions
template void upsi_tensor_lapack<base_device::DEVICE_CPU>(const Parallel_Orbitals* pv,
                                                          const int nband,
                                                          const int nlocal,
                                                          const ct::Tensor& U_operator,
                                                          const ct::Tensor& psi_k_laststep,
                                                          ct::Tensor& psi_k,
                                                          std::ofstream& ofs_running,
                                                          const int print_matrix);
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void upsi_tensor_lapack<base_device::DEVICE_GPU>(const Parallel_Orbitals* pv,
                                                          const int nband,
                                                          const int nlocal,
                                                          const ct::Tensor& U_operator,
                                                          const ct::Tensor& psi_k_laststep,
                                                          ct::Tensor& psi_k,
                                                          std::ofstream& ofs_running,
                                                          const int print_matrix);
#endif // __CUDA
#endif // __MPI
} // namespace module_tddft
