#include "norm_psi.h"

#include "source_base/block_cyclic_distribution_information.h"
#include "source_base/global_function.h" // ModuleBase::GlobalFunc::ZEROS
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"

#ifdef __CUBLASMP
#include "kernels/cuda/norm_psi_kernel.cuh"
#endif

#include <cassert>
#include <complex>
#include <cstddef>
#include <iostream>

namespace module_rt
{
#ifdef __MPI

inline int globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock = 0, gIndex = 0;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

void norm_psi(const BlockCyclicDistributionInformation& matrix_layout,
              const BlockCyclicDistributionInformation& wfc_layout,
              const BlockCyclicDistributionInformation& mo_repr_layout,
              const int nband,
              const int nlocal,
              const std::complex<double>* Stmp,
              std::complex<double>* psi_k,
              std::ofstream& ofs_running,
              const int print_matrix)
{
    assert(wfc_layout.size() > 0 && matrix_layout.size() > 0 && mo_repr_layout.size() > 0);

    std::complex<double>* tmp1 = new std::complex<double>[wfc_layout.size()];
    ModuleBase::GlobalFunc::ZEROS(tmp1, wfc_layout.size());

    std::complex<double>* mo_ovlp = new std::complex<double>[mo_repr_layout.size()];
    ModuleBase::GlobalFunc::ZEROS(mo_ovlp, mo_repr_layout.size());

    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nlocal,
                             1.0,
                             Stmp,
                             1,
                             1,
                             matrix_layout.desc(),
                             psi_k,
                             1,
                             1,
                             wfc_layout.desc(),
                             0.0,
                             tmp1,
                             1,
                             1,
                             wfc_layout.desc());

    ScalapackConnector::gemm('C',
                             'N',
                             nband,
                             nband,
                             nlocal,
                             1.0,
                             psi_k,
                             1,
                             1,
                             wfc_layout.desc(),
                             tmp1,
                             1,
                             1,
                             wfc_layout.desc(),
                             0.0,
                             mo_ovlp,
                             1,
                             1,
                             mo_repr_layout.desc());

    if (print_matrix)
    {
        ofs_running << "original Cij :" << std::endl;
        for (int i = 0; i < mo_repr_layout.ncloc(); i++)
        {
            const int in = i * mo_repr_layout.desc()[8];
            for (int j = 0; j < mo_repr_layout.nrloc(); j++)
            {
                double aa = mo_ovlp[in + j].real();
                double bb = mo_ovlp[in + j].imag();
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

    int naroc[2] = {0, 0}; // maximum number of row or column

    for (int iprow = 0; iprow < mo_repr_layout.nrproc(); ++iprow)
    {
        for (int ipcol = 0; ipcol < mo_repr_layout.ncproc(); ++ipcol)
        {
            if (mo_repr_layout.blacs_in_this_processor(iprow, ipcol))
            {
                naroc[0] = mo_repr_layout.nrloc();
                naroc[1] = mo_repr_layout.ncloc();
                const int lld = mo_repr_layout.desc()[8];
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, mo_repr_layout.nb(), mo_repr_layout.ncproc(), ipcol);
                    if (igcol >= nband)
                    {
                        continue;
                    }
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, mo_repr_layout.nb(), mo_repr_layout.nrproc(), iprow);
                        if (igrow >= nband)
                        {
                            continue;
                        }
                        if (igcol == igrow)
                        {
                            mo_ovlp[j * lld + i] = {1.0 / sqrt(mo_ovlp[j * lld + i].real()), 0.0};
                        }
                        else
                        {
                            mo_ovlp[j * lld + i] = {0.0, 0.0};
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow

    BlasConnector::copy(wfc_layout.size(), psi_k, 1, tmp1, 1);

    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nband,
                             1.0,
                             tmp1,
                             1,
                             1,
                             wfc_layout.desc(),
                             mo_ovlp,
                             1,
                             1,
                             mo_repr_layout.desc(),
                             0.0,
                             psi_k,
                             1,
                             1,
                             wfc_layout.desc());

    if (print_matrix)
    {
        ofs_running << " Cij:" << std::endl;
        for (int i = 0; i < mo_repr_layout.ncloc(); i++)
        {
            const int in = i * mo_repr_layout.desc()[8];
            for (int j = 0; j < mo_repr_layout.nrloc(); j++)
            {
                ofs_running << mo_ovlp[in + j].real() << "+" << mo_ovlp[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " psi_k:" << std::endl;
        for (int i = 0; i < wfc_layout.ncloc(); i++)
        {
            const int in = i * wfc_layout.desc()[8];
            for (int j = 0; j < wfc_layout.nrloc(); j++)
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
        ofs_running << " psi_k before normalization:" << std::endl;
        for (int i = 0; i < wfc_layout.ncloc(); i++)
        {
            const int in = i * wfc_layout.desc()[8];
            for (int j = 0; j < wfc_layout.nrloc(); j++)
            {
                double aa = tmp1[in + j].real();
                double bb = tmp1[in + j].imag();
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
        ofs_running << std::endl;
    }

    delete[] tmp1;
    delete[] mo_ovlp;
}

void norm_psi_tensor(const BlockCyclicDistributionInformation& matrix_layout,
                     const BlockCyclicDistributionInformation& wfc_layout,
                     const BlockCyclicDistributionInformation& mo_repr_layout,
                     const int nband,
                     const int nlocal,
                     const ct::Tensor& Stmp,
                     ct::Tensor& psi_k,
                     std::ofstream& ofs_running,
                     const int print_matrix,
                     CublasMpResources& cublas_res)
{
#ifdef __CUBLASMP
    if (!cublas_res.is_initialized || cublas_res.cublasmp_grid == nullptr)
    {
        return;
    }

    assert(wfc_layout.size() > 0 && matrix_layout.size() > 0 && mo_repr_layout.size() > 0);

    void* d_S = static_cast<void*>(const_cast<std::complex<double>*>(Stmp.data<std::complex<double>>()));
    void* d_Psi = static_cast<void*>(psi_k.data<std::complex<double>>());
    int64_t psi_elems = psi_k.NumElements();

    ct::Tensor Tmp1_gpu(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, ct::TensorShape({psi_elems}));
    void* d_Tmp1 = static_cast<void*>(Tmp1_gpu.data<std::complex<double>>());

    const std::size_t cij_elems = mo_repr_layout.size();
    ct::Tensor mo_ovlp_gpu(ct::DataType::DT_COMPLEX_DOUBLE,
                                     ct::DeviceType::GpuDevice,
                                     ct::TensorShape({static_cast<int64_t>(cij_elems)}));
    void* d_mo_ovlp = static_cast<void*>(mo_ovlp_gpu.data<std::complex<double>>());

    cudaMemsetAsync(d_mo_ovlp, 0, cij_elems * sizeof(std::complex<double>), cublas_res.stream);

    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta = {0.0, 0.0};

    cublasMpMatrixDescriptor_t desc_S, desc_Psi, desc_mo_ovlp;

    cublasMpMatrixDescriptorCreate(nlocal,
                                   nlocal,
                                   matrix_layout.desc()[4],
                                   matrix_layout.desc()[5],
                                   0,
                                   0,
                                   matrix_layout.desc()[8],
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_S);

    cublasMpMatrixDescriptorCreate(nlocal,
                                   nband,
                                   wfc_layout.desc()[4],
                                   wfc_layout.desc()[5],
                                   0,
                                   0,
                                   wfc_layout.desc()[8],
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_Psi);

    cublasMpMatrixDescriptorCreate(nband,
                                   nband,
                                   mo_repr_layout.desc()[4],
                                   mo_repr_layout.desc()[5],
                                   0,
                                   0,
                                   mo_repr_layout.desc()[8],
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_mo_ovlp);

    size_t ws_dev = 0, ws_host = 0;
    void *d_work = nullptr, *h_work = nullptr;

    // GEMM 1: S * Psi -> Tmp1
    cublasMpGemm_bufferSize(cublas_res.cublasmp_handle,
                            CUBLAS_OP_N,
                            CUBLAS_OP_N,
                            nlocal,
                            nband,
                            nlocal,
                            &alpha,
                            d_S,
                            1,
                            1,
                            desc_S,
                            d_Psi,
                            1,
                            1,
                            desc_Psi,
                            &beta,
                            d_Tmp1,
                            1,
                            1,
                            desc_Psi,
                            CUBLAS_COMPUTE_64F,
                            &ws_dev,
                            &ws_host);

    cudaMallocAsync(&d_work, ws_dev, cublas_res.stream);
    h_work = malloc(ws_host);

    cublasMpGemm(cublas_res.cublasmp_handle,
                 CUBLAS_OP_N,
                 CUBLAS_OP_N,
                 nlocal,
                 nband,
                 nlocal,
                 &alpha,
                 d_S,
                 1,
                 1,
                 desc_S,
                 d_Psi,
                 1,
                 1,
                 desc_Psi,
                 &beta,
                 d_Tmp1,
                 1,
                 1,
                 desc_Psi,
                 CUBLAS_COMPUTE_64F,
                 d_work,
                 ws_dev,
                 h_work,
                 ws_host);

    cudaFreeAsync(d_work, cublas_res.stream);
    free(h_work);

    // GEMM 2: Psi^H * Tmp1 -> mo_ovlp
    cublasMpGemm_bufferSize(cublas_res.cublasmp_handle,
                            CUBLAS_OP_C,
                            CUBLAS_OP_N,
                            nband,
                            nband,
                            nlocal,
                            &alpha,
                            d_Psi,
                            1,
                            1,
                            desc_Psi,
                            d_Tmp1,
                            1,
                            1,
                            desc_Psi,
                            &beta,
                            d_mo_ovlp,
                            1,
                            1,
                            desc_mo_ovlp,
                            CUBLAS_COMPUTE_64F,
                            &ws_dev,
                            &ws_host);

    cudaMallocAsync(&d_work, ws_dev, cublas_res.stream);
    h_work = malloc(ws_host);

    cublasMpGemm(cublas_res.cublasmp_handle,
                 CUBLAS_OP_C,
                 CUBLAS_OP_N,
                 nband,
                 nband,
                 nlocal,
                 &alpha,
                 d_Psi,
                 1,
                 1,
                 desc_Psi,
                 d_Tmp1,
                 1,
                 1,
                 desc_Psi,
                 &beta,
                 d_mo_ovlp,
                 1,
                 1,
                 desc_mo_ovlp,
                 CUBLAS_COMPUTE_64F,
                 d_work,
                 ws_dev,
                 h_work,
                 ws_host);

    cudaFreeAsync(d_work, cublas_res.stream);
    free(h_work);

    // Launch GPU In-place Normalization using the C++ wrapper
    const int mo_ovlp_elems = static_cast<int>(mo_repr_layout.size());
    module_rt::gpu::launch_normalize_cij_kernel(reinterpret_cast<cuDoubleComplex*>(d_mo_ovlp),
                                                mo_repr_layout.desc()[8],
                                                mo_ovlp_elems,
                                                mo_repr_layout.desc()[4],
                                                mo_repr_layout.nrproc(),
                                                mo_repr_layout.ncproc(),
                                                mo_repr_layout.irproc(),
                                                mo_repr_layout.icproc(),
                                                nband,
                                                cublas_res.stream);

    // GEMM 3: Tmp1 * mo_ovlp -> Psi
    cudaMemcpyAsync(d_Tmp1,
                    d_Psi,
                    psi_elems * sizeof(std::complex<double>),
                    cudaMemcpyDeviceToDevice,
                    cublas_res.stream);

    cublasMpGemm_bufferSize(cublas_res.cublasmp_handle,
                            CUBLAS_OP_N,
                            CUBLAS_OP_N,
                            nlocal,
                            nband,
                            nband,
                            &alpha,
                            d_Tmp1,
                            1,
                            1,
                            desc_Psi,
                            d_mo_ovlp,
                            1,
                            1,
                            desc_mo_ovlp,
                            &beta,
                            d_Psi,
                            1,
                            1,
                            desc_Psi,
                            CUBLAS_COMPUTE_64F,
                            &ws_dev,
                            &ws_host);

    cudaMallocAsync(&d_work, ws_dev, cublas_res.stream);
    h_work = malloc(ws_host);

    cublasMpGemm(cublas_res.cublasmp_handle,
                 CUBLAS_OP_N,
                 CUBLAS_OP_N,
                 nlocal,
                 nband,
                 nband,
                 &alpha,
                 d_Tmp1,
                 1,
                 1,
                 desc_Psi,
                 d_mo_ovlp,
                 1,
                 1,
                 desc_mo_ovlp,
                 &beta,
                 d_Psi,
                 1,
                 1,
                 desc_Psi,
                 CUBLAS_COMPUTE_64F,
                 d_work,
                 ws_dev,
                 h_work,
                 ws_host);

    cudaStreamSynchronize(cublas_res.stream);

    cublasMpMatrixDescriptorDestroy(desc_S);
    cublasMpMatrixDescriptorDestroy(desc_Psi);
    cublasMpMatrixDescriptorDestroy(desc_mo_ovlp);

    cudaFreeAsync(d_work, cublas_res.stream);
    free(h_work);
#endif // __CUBLASMP
}

template <typename Device>
void norm_psi_tensor_lapack(const int nband,
                            const int nlocal,
                            const ct::Tensor& Stmp,
                            ct::Tensor& psi_k,
                            std::ofstream& ofs_running,
                            const int print_matrix)
{
    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    // Create Tensor objects for temporary data
    ct::Tensor tmp1(
        ct::DataType::DT_COMPLEX_DOUBLE,
        ct_device_type,
        ct::TensorShape({nlocal * nband})); // tmp1 shape: nlocal * nband
    tmp1.zero();

    ct::Tensor mo_ovlp(ct::DataType::DT_COMPLEX_DOUBLE,
                            ct_device_type,
                            ct::TensorShape({nlocal * nlocal})); // mo_ovlp shape: nlocal * nlocal
    mo_ovlp.zero();

    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta = {0.0, 0.0};

    // Perform matrix multiplication: tmp1 = Stmp * psi_k
    ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                              'N',
                                                              nlocal,
                                                              nband,
                                                              nlocal,
                                                              &alpha,
                                                              Stmp.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of Stmp
                                                              psi_k.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of psi_k
                                                              &beta,
                                                              tmp1.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of tmp1

    // Perform matrix multiplication: mo_ovlp = psi_k^dagger * tmp1
    ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('C',
                                                              'N',
                                                              nband,
                                                              nband,
                                                              nlocal,
                                                              &alpha,
                                                              psi_k.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of psi_k
                                                              tmp1.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of tmp1
                                                              &beta,
                                                              mo_ovlp.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of mo_ovlp

    // Normalize mo_ovlp: set diagonal elements to 1/sqrt(mo_ovlp[i][i]), off-diagonal elements to 0
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        // Step 1: Copy mo_ovlp from GPU to CPU
        ct::Tensor mo_ovlp_cpu = mo_ovlp.to_device<ct::DEVICE_CPU>();

        // Step 2: Perform normalization on CPU
        for (int i = 0; i < nband; ++i)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nband; ++j)
            {
                if (i == j)
                {
                    mo_ovlp_cpu.data<std::complex<double>>()[in + j]
                        = {1.0 / sqrt(mo_ovlp_cpu.data<std::complex<double>>()[in + j].real()), 0.0};
                }
                else
                {
                    mo_ovlp_cpu.data<std::complex<double>>()[in + j] = {0.0, 0.0};
                }
            }
        }

        // Step 3: Copy normalized mo_ovlp back to GPU
        mo_ovlp = mo_ovlp_cpu.to_device<ct_Device>();
    }
    else
    {
        // CPU implementation
        for (int i = 0; i < nband; ++i)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nband; ++j)
            {
                if (i == j)
                {
                    mo_ovlp.data<std::complex<double>>()[in + j]
                        = {1.0 / sqrt(mo_ovlp.data<std::complex<double>>()[in + j].real()), 0.0};
                }
                else
                {
                    mo_ovlp.data<std::complex<double>>()[in + j] = {0.0, 0.0};
                }
            }
        }
    }

    // Copy psi_k to tmp1 (using deep copy)
    // tmp1.CopyFrom(psi_k); // Does not work because this will cause tmp1 and psi_k to share the same data
    tmp1 = psi_k; // operator= overload for Tensor class

    // Perform matrix multiplication: psi_k = tmp1 * mo_ovlp
    ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                              'N',
                                                              nlocal,
                                                              nband,
                                                              nband,
                                                              &alpha,
                                                              tmp1.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of tmp1
                                                              mo_ovlp.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of mo_ovlp
                                                              &beta,
                                                              psi_k.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of psi_k
}

// Explicit instantiation of template functions
template void norm_psi_tensor_lapack<base_device::DEVICE_CPU>(const int nband,
                                                              const int nlocal,
                                                              const ct::Tensor& Stmp,
                                                              ct::Tensor& psi_k,
                                                              std::ofstream& ofs_running,
                                                              const int print_matrix);
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void norm_psi_tensor_lapack<base_device::DEVICE_GPU>(const int nband,
                                                              const int nlocal,
                                                              const ct::Tensor& Stmp,
                                                              ct::Tensor& psi_k,
                                                              std::ofstream& ofs_running,
                                                              const int print_matrix);
#endif // __CUDA
#endif // __MPI
} // namespace module_rt
