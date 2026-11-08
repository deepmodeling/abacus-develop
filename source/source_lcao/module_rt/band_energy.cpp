#include "band_energy.h"

#include "source_base/parallel_blkcyc.h"
#include "source_base/global_function.h"
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_device/memory_op.h"
#include "source_base/module_external/scalapack_connector.h"

#ifdef __CUBLASMP
#include "kernels/cuda/band_energy_kernel.cuh"
#endif

#include <algorithm>
#include <complex>
#include <cstddef>
#include <fstream>
#include <vector>

namespace module_rt
{
#ifdef __MPI
namespace
{
void gather_mo_hamilt(const BlockCyclicDistribution& mo_repr_layout,
                                  const int nband,
                                  std::complex<double>* mo_hamilt,
                                  const EijOutputInfo& output_info,
                                  std::ofstream& ofs)
{
    const bool is_root = mo_repr_layout.blacs_in_this_processor(0, 0);
    const size_t buffer_size = is_root ? static_cast<size_t>(nband) * nband : 1;
    std::vector<std::complex<double>> gathered(buffer_size);
    int source_desc[9];
    std::copy(mo_repr_layout.desc(), mo_repr_layout.desc() + 9, source_desc);
    int target_desc[9] = {1, source_desc[1], nband, nband, nband, nband, 0, 0, nband};
    // All grid ranks participate, including ranks whose running log is closed.
    Cpxgemr2d(nband, nband, mo_hamilt, 1, 1, source_desc, gathered.data(), 1, 1, target_desc, source_desc[1]);
    if (is_root)
    {
        ModuleIO::print_eij(gathered.data(), nband, nband, output_info, ofs);
    }
}
} // namespace

inline int globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

void compute_ekb(const BlockCyclicDistribution& matrix_layout,
                 const BlockCyclicDistribution& wfc_layout,
                 const BlockCyclicDistribution& mo_repr_layout,
                 const int nband,
                 const int nlocal,
                 const std::complex<double>* Htmp,
                 const std::complex<double>* psi_k,
                 double* ekb,
                 std::ofstream& ofs_running,
                 const EijOutputInfo& output_info)
{
    assert(wfc_layout.size() > 0 && matrix_layout.size() > 0 && mo_repr_layout.size() > 0);

    std::complex<double>* tmp1 = new std::complex<double>[wfc_layout.size()];
    ModuleBase::GlobalFunc::ZEROS(tmp1, wfc_layout.size());

    std::complex<double>* mo_hamilt = new std::complex<double>[mo_repr_layout.size()];
    ModuleBase::GlobalFunc::ZEROS(mo_hamilt, mo_repr_layout.size());

    ScalapackConnector::
        gemm('N', 'N', nlocal, nband, nlocal, 1.0, Htmp, 1, 1, matrix_layout.desc(), psi_k, 1, 1, wfc_layout.desc(), 0.0, tmp1, 1, 1, wfc_layout.desc());

    ScalapackConnector::
        gemm('C', 'N', nband, nband, nlocal, 1.0, psi_k, 1, 1, wfc_layout.desc(), tmp1, 1, 1, wfc_layout.desc(), 0.0, mo_hamilt, 1, 1, mo_repr_layout.desc());

    if (output_info.threshold >= 0.0)
    {
        gather_mo_hamilt(mo_repr_layout, nband, mo_hamilt, output_info, ofs_running);
    }

    int info = 0;
    int naroc[2] = {0, 0};

    assert(nband > 0);
    double* eii = new double[nband];
    ModuleBase::GlobalFunc::ZEROS(eii, nband);

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
                            eii[igcol] = mo_hamilt[j * lld + i].real();
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow
    info = MPI_Allreduce(eii, ekb, nband, MPI_DOUBLE, MPI_SUM, mo_repr_layout.comm());

    delete[] tmp1;
    delete[] mo_hamilt;
    delete[] eii;
}

void compute_ekb_tensor(const BlockCyclicDistribution& matrix_layout,
                        const BlockCyclicDistribution& wfc_layout,
                        const BlockCyclicDistribution& mo_repr_layout,
                        const int nband,
                        const int nlocal,
                        const ct::Tensor& Htmp,
                        const ct::Tensor& psi_k,
                        ct::Tensor& ekb,
                        std::ofstream& ofs_running,
                        CublasMpResources& cublas_res)
{
#ifdef __CUBLASMP
    // 1. Resource validation
    if (!cublas_res.is_initialized || cublas_res.cublasmp_grid == nullptr)
    {
        return;
    }

    assert(wfc_layout.size() > 0 && matrix_layout.size() > 0 && mo_repr_layout.size() > 0);
    assert(Htmp.device_type() == ct::DeviceType::GpuDevice);
    assert(psi_k.device_type() == ct::DeviceType::GpuDevice);
    assert(ekb.device_type() == ct::DeviceType::GpuDevice);

    // 2. Data Pointers
    void* d_H = static_cast<void*>(const_cast<std::complex<double>*>(Htmp.data<std::complex<double>>()));
    void* d_Psi = static_cast<void*>(const_cast<std::complex<double>*>(psi_k.data<std::complex<double>>()));

    int64_t psi_elems = psi_k.NumElements();
    ct::Tensor Tmp1_gpu(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, ct::TensorShape({psi_elems}));
    void* d_Tmp1 = static_cast<void*>(Tmp1_gpu.data<std::complex<double>>());

    const std::size_t mo_hamilt_elems = mo_repr_layout.size();
    ct::Tensor mo_hamilt_gpu(ct::DataType::DT_COMPLEX_DOUBLE,
                                         ct::DeviceType::GpuDevice,
                                         ct::TensorShape({static_cast<int64_t>(mo_hamilt_elems)}));
    void* d_mo_hamilt = static_cast<void*>(mo_hamilt_gpu.data<std::complex<double>>());

    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta = {0.0, 0.0};

    // 3. Matrix Descriptors Creation
    cublasMpMatrixDescriptor_t desc_H, desc_Psi, desc_mo_hamilt;

    // H descriptor: nlocal x nlocal
    cublasMpMatrixDescriptorCreate(matrix_layout.desc()[2],
                                   matrix_layout.desc()[3],
                                   matrix_layout.desc()[4],
                                   matrix_layout.desc()[5],
                                   0,
                                   0,
                                   matrix_layout.desc()[8],
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_H);

    // Psi descriptor: nlocal x nband
    cublasMpMatrixDescriptorCreate(wfc_layout.desc()[2],
                                   wfc_layout.desc()[3],
                                   wfc_layout.desc()[4],
                                   wfc_layout.desc()[5],
                                   0,
                                   0,
                                   wfc_layout.desc()[8],
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_Psi);

    // MO representation descriptor
    cublasMpMatrixDescriptorCreate(nband,
                                   nband,
                                   mo_repr_layout.desc()[4],
                                   mo_repr_layout.desc()[5],
                                   0,
                                   0,
                                   mo_repr_layout.desc()[8],
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_mo_hamilt);

    size_t ws_dev = 0, ws_host = 0;
    void *d_work = nullptr, *h_work = nullptr;

    // 4. GEMM 1: Tmp1 = H * Psi
    cublasMpGemm_bufferSize(cublas_res.cublasmp_handle,
                            CUBLAS_OP_N,
                            CUBLAS_OP_N,
                            matrix_layout.desc()[2],
                            wfc_layout.desc()[3],
                            matrix_layout.desc()[3],
                            &alpha,
                            d_H,
                            1,
                            1,
                            desc_H,
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
                 matrix_layout.desc()[2],
                 wfc_layout.desc()[3],
                 matrix_layout.desc()[3],
                 &alpha,
                 d_H,
                 1,
                 1,
                 desc_H,
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

    // 5. GEMM 2: Eij = Psi^H * Tmp1
    cublasMpGemm_bufferSize(cublas_res.cublasmp_handle,
                            CUBLAS_OP_C,
                            CUBLAS_OP_N,
                            wfc_layout.desc()[3],
                            wfc_layout.desc()[3],
                            wfc_layout.desc()[2],
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
                            d_mo_hamilt,
                            1,
                            1,
                            desc_mo_hamilt,
                            CUBLAS_COMPUTE_64F,
                            &ws_dev,
                            &ws_host);

    cudaMallocAsync(&d_work, ws_dev, cublas_res.stream);
    h_work = malloc(ws_host);

    cublasMpGemm(cublas_res.cublasmp_handle,
                 CUBLAS_OP_C,
                 CUBLAS_OP_N,
                 wfc_layout.desc()[3],
                 wfc_layout.desc()[3],
                 wfc_layout.desc()[2],
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
                 d_mo_hamilt,
                 1,
                 1,
                 desc_mo_hamilt,
                 CUBLAS_COMPUTE_64F,
                 d_work,
                 ws_dev,
                 h_work,
                 ws_host);

    cudaFreeAsync(d_work, cublas_res.stream);
    free(h_work);

    // 6. Extract Diagonal directly on GPU
    // Prepare a zero-initialized buffer on GPU to store the local parts of the diagonal
    ct::Tensor eii_gpu(ct::DataType::DT_DOUBLE, ct::DeviceType::GpuDevice, ct::TensorShape({nband}));
    double* d_eii = static_cast<double*>(eii_gpu.data<double>());
    cudaMemsetAsync(d_eii, 0, nband * sizeof(double), cublas_res.stream);

    // Launch the extraction kernel
    const int mo_repr_elems = static_cast<int>(mo_repr_layout.size());
    module_rt::gpu::launch_extract_ekb_kernel(reinterpret_cast<cuDoubleComplex*>(d_mo_hamilt),
                                              d_eii,
                                              mo_repr_layout.desc()[8],
                                              mo_repr_elems,
                                              mo_repr_layout.desc()[4],
                                              mo_repr_layout.nrproc(),
                                              mo_repr_layout.ncproc(),
                                              mo_repr_layout.irproc(),
                                              mo_repr_layout.icproc(),
                                              nband,
                                              cublas_res.stream);

    // 7. CUDA-aware MPI Reduction
    // VERY IMPORTANT: We must synchronize the stream before passing the GPU pointer
    // to MPI, because MPI operations are generally synchronous to the CPU thread.
    cudaStreamSynchronize(cublas_res.stream);

    double* d_ekb = static_cast<double*>(ekb.data<double>());

    // Direct GPU-to-GPU reduction using CUDA-aware MPI
    MPI_Allreduce(d_eii, d_ekb, nband, MPI_DOUBLE, MPI_SUM, mo_repr_layout.comm());

    // 8. Cleanup
    cublasMpMatrixDescriptorDestroy(desc_H);
    cublasMpMatrixDescriptorDestroy(desc_Psi);
    cublasMpMatrixDescriptorDestroy(desc_mo_hamilt);
#endif // __CUBLASMP
}

template <typename Device>
void compute_ekb_tensor_lapack(const int nband,
                               const int nlocal,
                               const ct::Tensor& Htmp,
                               const ct::Tensor& psi_k,
                               ct::Tensor& ekb,
                               std::ofstream& ofs_running,
                               const EijOutputInfo& output_info)
{
    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    // Create Tensor objects for temporary data
    ct::Tensor tmp1(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal * nband})); // tmp1 shape: nlocal * nband
    tmp1.zero();

    ct::Tensor mo_hamilt(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal * nlocal})); // mo_hamilt shape: nlocal * nlocal
    // Why not use nband * nband ?????
    mo_hamilt.zero();

    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta = {0.0, 0.0};

    // Perform matrix multiplication: tmp1 = Htmp * psi_k
    ct::kernels::blas_gemm<std::complex<double>, ct_Device>()('N',
                                                              'N',
                                                              nlocal,
                                                              nband,
                                                              nlocal,
                                                              &alpha,
                                                              Htmp.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of Htmp
                                                              psi_k.data<std::complex<double>>(),
                                                              nlocal, // Leading dimension of psi_k
                                                              &beta,
                                                              tmp1.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of tmp1

    // Perform matrix multiplication: mo_hamilt = psi_k^dagger * tmp1
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
                                                              mo_hamilt.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of mo_hamilt

    if (output_info.threshold >= 0.0)
    {
        ct::Tensor mo_hamilt_cpu = mo_hamilt.to_device<ct::DEVICE_CPU>();
        ModuleIO::print_eij(mo_hamilt_cpu.data<std::complex<double>>(), nband, nlocal, output_info, ofs_running);
    }

    // Extract diagonal elements of mo_hamilt into ekb
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        // GPU implementation
        for (int i = 0; i < nband; ++i)
        {
            base_device::memory::synchronize_memory_op<double, Device, Device>()(
                ekb.data<double>() + i,
                reinterpret_cast<const double*>(mo_hamilt.data<std::complex<double>>() + i * nlocal + i),
                1);
        }
    }
    else
    {
        // CPU implementation
        for (int i = 0; i < nband; ++i)
        {
            ekb.data<double>()[i] = mo_hamilt.data<std::complex<double>>()[i * nlocal + i].real();
        }
    }
}

// Explicit instantiation of template functions
template void compute_ekb_tensor_lapack<base_device::DEVICE_CPU>(const int nband,
                                                                 const int nlocal,
                                                                 const ct::Tensor& Htmp,
                                                                 const ct::Tensor& psi_k,
                                                                 ct::Tensor& ekb,
                                                                 std::ofstream& ofs_running,
                                                                 const EijOutputInfo& output_info);

#if ((defined __CUDA) /* || (defined __ROCM) */)
template void compute_ekb_tensor_lapack<base_device::DEVICE_GPU>(const int nband,
                                                                 const int nlocal,
                                                                 const ct::Tensor& Htmp,
                                                                 const ct::Tensor& psi_k,
                                                                 ct::Tensor& ekb,
                                                                 std::ofstream& ofs_running,
                                                                 const EijOutputInfo& output_info);
#endif // __CUDA
#endif // __MPI

} // namespace module_rt
