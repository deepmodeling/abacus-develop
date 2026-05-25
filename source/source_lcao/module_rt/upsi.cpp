#include "upsi.h"

#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"

#include <cassert>
#include <complex>
#include <iostream>

namespace module_rt
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
                 const int print_matrix,
                 CublasMpResources& cublas_res)
{
#ifdef __CUBLASMP
    // 1. Resource validation
    if (!cublas_res.is_initialized || cublas_res.cublasmp_grid == nullptr)
    {
        return;
    }

    assert(U_operator.device_type() == ct::DeviceType::GpuDevice);
    assert(psi_k_laststep.device_type() == ct::DeviceType::GpuDevice);
    assert(psi_k.device_type() == ct::DeviceType::GpuDevice);

    // 2. Extract device pointers
    void* d_U = static_cast<void*>(const_cast<std::complex<double>*>(U_operator.data<std::complex<double>>()));
    void* d_Psi_old
        = static_cast<void*>(const_cast<std::complex<double>*>(psi_k_laststep.data<std::complex<double>>()));
    void* d_Psi_k = static_cast<void*>(psi_k.data<std::complex<double>>());

    // 3. Create matrix descriptor for U operator (N x N)
    int64_t m_u = pv->desc[2];
    int64_t n_u = pv->desc[3];
    int64_t mb_u = pv->desc[4];
    int64_t nb_u = pv->desc[5];
    int64_t lld_u = pv->desc[8];

    cublasMpMatrixDescriptor_t desc_U;
    cublasMpMatrixDescriptorCreate(m_u, n_u, mb_u, nb_u, 0, 0, lld_u, CUDA_C_64F, cublas_res.cublasmp_grid, &desc_U);

    // 4. Create matrix descriptor for Psi (N x nband)
    int64_t m_psi = pv->desc_wfc[2];
    int64_t n_psi = pv->desc_wfc[3];
    int64_t mb_psi = pv->desc_wfc[4];
    int64_t nb_psi = pv->desc_wfc[5];
    int64_t lld_psi = pv->desc_wfc[8];

    cublasMpMatrixDescriptor_t desc_Psi;
    cublasMpMatrixDescriptorCreate(m_psi,
                                   n_psi,
                                   mb_psi,
                                   nb_psi,
                                   0,
                                   0,
                                   lld_psi,
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_Psi);

    // 5. Query workspace size for GEMM: Psi_k = 1.0 * U * Psi_old + 0.0 * Psi_k
    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta = {0.0, 0.0};
    size_t ws_gemm_dev = 0;
    size_t ws_gemm_host = 0;

    cublasMpGemm_bufferSize(cublas_res.cublasmp_handle,
                            CUBLAS_OP_N,
                            CUBLAS_OP_N,
                            m_u,
                            n_psi,
                            n_u,
                            &alpha,
                            d_U,
                            1,
                            1,
                            desc_U,
                            d_Psi_old,
                            1,
                            1,
                            desc_Psi,
                            &beta,
                            d_Psi_k,
                            1,
                            1,
                            desc_Psi,
                            CUBLAS_COMPUTE_64F,
                            &ws_gemm_dev,
                            &ws_gemm_host);

    void* d_work = nullptr;
    void* h_work = nullptr;

    cudaMallocAsync(&d_work, ws_gemm_dev, cublas_res.stream);
    h_work = malloc(ws_gemm_host);

    // 6. Execute GEMM
    cublasMpGemm(cublas_res.cublasmp_handle,
                 CUBLAS_OP_N,
                 CUBLAS_OP_N,
                 m_u,
                 n_psi,
                 n_u,
                 &alpha,
                 d_U,
                 1,
                 1,
                 desc_U,
                 d_Psi_old,
                 1,
                 1,
                 desc_Psi,
                 &beta,
                 d_Psi_k,
                 1,
                 1,
                 desc_Psi,
                 CUBLAS_COMPUTE_64F,
                 d_work,
                 ws_gemm_dev,
                 h_work,
                 ws_gemm_host);

    // 7. Synchronize and clean up resources
    cudaStreamSynchronize(cublas_res.stream);

    cublasMpMatrixDescriptorDestroy(desc_U);
    cublasMpMatrixDescriptorDestroy(desc_Psi);
    cudaFreeAsync(d_work, cublas_res.stream);
    free(h_work);
#endif // __CUBLASMP
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
} // namespace module_rt
