#include "band_energy.h"

#include "evolve_elec.h"
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_external/scalapack_connector.h"

#ifdef __CUBLASMP
#include "kernels/cuda/band_energy_kernel.cuh"
#endif

#include <complex>
#include <iostream>

namespace module_rt
{
#ifdef __MPI

inline int globalIndex(int localindex, int nblk, int nprocs, int myproc)
{
    int iblock, gIndex;
    iblock = localindex / nblk;
    gIndex = (iblock * nprocs + myproc) * nblk + localindex % nblk;
    return gIndex;
}

void compute_ekb(const Parallel_Orbitals* pv,
                 const int nband,
                 const int nlocal,
                 const std::complex<double>* Htmp,
                 const std::complex<double>* psi_k,
                 double* ekb,
                 std::ofstream& ofs_running)
{
    assert(pv->nloc_wfc > 0 && pv->nloc > 0);

    std::complex<double>* tmp1 = new std::complex<double>[pv->nloc_wfc];
    ModuleBase::GlobalFunc::ZEROS(tmp1, pv->nloc_wfc);

    std::complex<double>* eij = new std::complex<double>[pv->nloc];
    ModuleBase::GlobalFunc::ZEROS(eij, pv->nloc);

    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nband,
                             nlocal,
                             1.0,
                             Htmp,
                             1,
                             1,
                             pv->desc,
                             psi_k,
                             1,
                             1,
                             pv->desc_wfc,
                             0.0,
                             tmp1,
                             1,
                             1,
                             pv->desc_wfc);

    ScalapackConnector::gemm('C',
                             'N',
                             nband,
                             nband,
                             nlocal,
                             1.0,
                             psi_k,
                             1,
                             1,
                             pv->desc_wfc,
                             tmp1,
                             1,
                             1,
                             pv->desc_wfc,
                             0.0,
                             eij,
                             1,
                             1,
                             pv->desc_Eij);

    if (PARAM.inp.td_print_eij > 0.0)
    {
        ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
        ofs_running << " Eij:" << std::endl;
        for (int i = 0; i < pv->nrow_bands; i++)
        {
            const int in = i * pv->ncol;
            for (int j = 0; j < pv->ncol_bands; j++)
            {
                double aa = eij[in + j].real();
                double bb = eij[in + j].imag();
                if (std::abs(aa) < PARAM.inp.td_print_eij)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < PARAM.inp.td_print_eij)
                {
                    bb = 0.0;
                }
                if (std::abs(aa) > 0.0 || std::abs(bb) > 0.0)
                {
                    std::streamsize original_precision = ofs_running.precision();
                    ofs_running << std::fixed << std::setprecision(8);
                    ofs_running << "i = " << std::setw(2) << i << ", j = " << std::setw(2) << j
                                << ", Eij = " << std::setw(12) << aa << " + " << std::setw(12) << bb << " i"
                                << std::endl;
                    ofs_running.unsetf(std::ios_base::fixed);
                    ofs_running.precision(original_precision);
                }
            }
        }
        ofs_running << std::endl;
        ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
    }

    int info = 0;
    int naroc[2] = {0, 0};

    assert(nband > 0);
    double* eii = new double[nband];
    ModuleBase::GlobalFunc::ZEROS(eii, nband);

    for (int iprow = 0; iprow < pv->dim0; ++iprow)
    {
        for (int ipcol = 0; ipcol < pv->dim1; ++ipcol)
        {
            if (iprow == pv->coord[0] && ipcol == pv->coord[1])
            {
                naroc[0] = pv->nrow;
                naroc[1] = pv->ncol;
                for (int j = 0; j < naroc[1]; ++j)
                {
                    int igcol = globalIndex(j, pv->nb, pv->dim1, ipcol);
                    if (igcol >= nband)
                    {
                        continue;
                    }
                    for (int i = 0; i < naroc[0]; ++i)
                    {
                        int igrow = globalIndex(i, pv->nb, pv->dim0, iprow);
                        if (igrow >= nband)
                        {
                            continue;
                        }
                        if (igcol == igrow)
                        {
                            eii[igcol] = eij[j * naroc[0] + i].real();
                        }
                    }
                }
            }
        } // loop ipcol
    } // loop iprow
    info = MPI_Allreduce(eii, ekb, nband, MPI_DOUBLE, MPI_SUM, pv->comm());

    delete[] tmp1;
    delete[] eij;
    delete[] eii;
}

void compute_ekb_tensor(const Parallel_Orbitals* pv,
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

    assert(pv->nloc_wfc > 0 && pv->nloc > 0);
    assert(Htmp.device_type() == ct::DeviceType::GpuDevice);
    assert(psi_k.device_type() == ct::DeviceType::GpuDevice);
    assert(ekb.device_type() == ct::DeviceType::GpuDevice);

    // 2. Data Pointers
    void* d_H = static_cast<void*>(const_cast<std::complex<double>*>(Htmp.data<std::complex<double>>()));
    void* d_Psi = static_cast<void*>(const_cast<std::complex<double>*>(psi_k.data<std::complex<double>>()));

    int64_t psi_elems = psi_k.NumElements();
    ct::Tensor Tmp1_gpu(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, ct::TensorShape({psi_elems}));
    void* d_Tmp1 = static_cast<void*>(Tmp1_gpu.data<std::complex<double>>());

    int64_t eij_elems = pv->nloc;
    ct::Tensor Eij_gpu(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, ct::TensorShape({eij_elems}));
    void* d_Eij = static_cast<void*>(Eij_gpu.data<std::complex<double>>());

    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta = {0.0, 0.0};

    // 3. Matrix Descriptors Creation
    cublasMpMatrixDescriptor_t desc_H, desc_Psi, desc_Eij;

    // H descriptor: nlocal x nlocal
    cublasMpMatrixDescriptorCreate(pv->desc[2],
                                   pv->desc[3],
                                   pv->desc[4],
                                   pv->desc[5],
                                   0,
                                   0,
                                   pv->desc[8],
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_H);

    // Psi descriptor: nlocal x nband
    cublasMpMatrixDescriptorCreate(pv->desc_wfc[2],
                                   pv->desc_wfc[3],
                                   pv->desc_wfc[4],
                                   pv->desc_wfc[5],
                                   0,
                                   0,
                                   pv->desc_wfc[8],
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_Psi);

    // Eij descriptor: MUST use nband x nband physically, to match pv->desc_Eij expectations
    cublasMpMatrixDescriptorCreate(nband,
                                   nband,
                                   pv->desc_Eij[4],
                                   pv->desc_Eij[5],
                                   0,
                                   0,
                                   pv->desc_Eij[8],
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_Eij);

    size_t ws_dev = 0, ws_host = 0;
    void *d_work = nullptr, *h_work = nullptr;

    // 4. GEMM 1: Tmp1 = H * Psi
    cublasMpGemm_bufferSize(cublas_res.cublasmp_handle,
                            CUBLAS_OP_N,
                            CUBLAS_OP_N,
                            pv->desc[2],
                            pv->desc_wfc[3],
                            pv->desc[3],
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
                 pv->desc[2],
                 pv->desc_wfc[3],
                 pv->desc[3],
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
                            pv->desc_wfc[3],
                            pv->desc_wfc[3],
                            pv->desc_wfc[2],
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
                            d_Eij,
                            1,
                            1,
                            desc_Eij,
                            CUBLAS_COMPUTE_64F,
                            &ws_dev,
                            &ws_host);

    cudaMallocAsync(&d_work, ws_dev, cublas_res.stream);
    h_work = malloc(ws_host);

    cublasMpGemm(cublas_res.cublasmp_handle,
                 CUBLAS_OP_C,
                 CUBLAS_OP_N,
                 pv->desc_wfc[3],
                 pv->desc_wfc[3],
                 pv->desc_wfc[2],
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
                 d_Eij,
                 1,
                 1,
                 desc_Eij,
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
    module_rt::gpu::launch_extract_ekb_kernel(reinterpret_cast<cuDoubleComplex*>(d_Eij),
                                              d_eii,
                                              pv->desc_Eij[8],
                                              pv->nloc,
                                              pv->desc_Eij[4],
                                              pv->dim0,
                                              pv->dim1,
                                              pv->coord[0],
                                              pv->coord[1],
                                              nband,
                                              cublas_res.stream);

    // 7. CUDA-aware MPI Reduction
    // VERY IMPORTANT: We must synchronize the stream before passing the GPU pointer
    // to MPI, because MPI operations are generally synchronous to the CPU thread.
    cudaStreamSynchronize(cublas_res.stream);

    double* d_ekb = static_cast<double*>(ekb.data<double>());

    // Direct GPU-to-GPU reduction using CUDA-aware MPI
    MPI_Allreduce(d_eii, d_ekb, nband, MPI_DOUBLE, MPI_SUM, pv->comm());

    // 8. Cleanup
    cublasMpMatrixDescriptorDestroy(desc_H);
    cublasMpMatrixDescriptorDestroy(desc_Psi);
    cublasMpMatrixDescriptorDestroy(desc_Eij);
#endif // __CUBLASMP
}

template <typename Device>
void compute_ekb_tensor_lapack(const Parallel_Orbitals* pv,
                               const int nband,
                               const int nlocal,
                               const ct::Tensor& Htmp,
                               const ct::Tensor& psi_k,
                               ct::Tensor& ekb,
                               std::ofstream& ofs_running)
{
    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    // Create Tensor objects for temporary data
    ct::Tensor tmp1(ct::DataType::DT_COMPLEX_DOUBLE,
                    ct_device_type,
                    ct::TensorShape({nlocal * nband})); // tmp1 shape: nlocal * nband
    tmp1.zero();

    ct::Tensor eij(ct::DataType::DT_COMPLEX_DOUBLE,
                   ct_device_type,
                   ct::TensorShape({nlocal * nlocal})); // eij shape: nlocal * nlocal
    // Why not use nband * nband ?????
    eij.zero();

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

    // Perform matrix multiplication: eij = psi_k^dagger * tmp1
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
                                                              eij.data<std::complex<double>>(),
                                                              nlocal); // Leading dimension of eij

    if (PARAM.inp.td_print_eij >= 0.0)
    {
        ct::Tensor eij_cpu = eij.to_device<ct::DEVICE_CPU>();

        ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
        ofs_running << " Eij:" << std::endl;
        for (int i = 0; i < nband; i++)
        {
            const int in = i * nlocal;
            for (int j = 0; j < nband; j++)
            {
                double aa = eij_cpu.data<std::complex<double>>()[in + j].real();
                double bb = eij_cpu.data<std::complex<double>>()[in + j].imag();
                if (std::abs(aa) < PARAM.inp.td_print_eij)
                {
                    aa = 0.0;
                }
                if (std::abs(bb) < PARAM.inp.td_print_eij)
                {
                    bb = 0.0;
                }
                if (std::abs(aa) > 0.0 || std::abs(bb) > 0.0)
                {
                    std::streamsize original_precision = ofs_running.precision();
                    ofs_running << std::fixed << std::setprecision(8);
                    ofs_running << "i = " << std::setw(2) << i << ", j = " << std::setw(2) << j
                                << ", Eij = " << std::setw(12) << aa << " + " << std::setw(12) << bb << " i"
                                << std::endl;
                    ofs_running.unsetf(std::ios_base::fixed);
                    ofs_running.precision(original_precision);
                }
            }
        }
        ofs_running << std::endl;
        ofs_running
            << "------------------------------------------------------------------------------------------------"
            << std::endl;
    }

    // Extract diagonal elements of eij into ekb
    if (ct_device_type == ct::DeviceType::GpuDevice)
    {
        // GPU implementation
        for (int i = 0; i < nband; ++i)
        {
            base_device::memory::synchronize_memory_op<double, Device, Device>()(
                ekb.data<double>() + i,
                reinterpret_cast<const double*>(eij.data<std::complex<double>>() + i * nlocal + i),
                1);
        }
    }
    else
    {
        // CPU implementation
        for (int i = 0; i < nband; ++i)
        {
            ekb.data<double>()[i] = eij.data<std::complex<double>>()[i * nlocal + i].real();
        }
    }
}

// Explicit instantiation of template functions
template void compute_ekb_tensor_lapack<base_device::DEVICE_CPU>(const Parallel_Orbitals* pv,
                                                                 const int nband,
                                                                 const int nlocal,
                                                                 const ct::Tensor& Htmp,
                                                                 const ct::Tensor& psi_k,
                                                                 ct::Tensor& ekb,
                                                                 std::ofstream& ofs_running);

#if ((defined __CUDA) /* || (defined __ROCM) */)
template void compute_ekb_tensor_lapack<base_device::DEVICE_GPU>(const Parallel_Orbitals* pv,
                                                                 const int nband,
                                                                 const int nlocal,
                                                                 const ct::Tensor& Htmp,
                                                                 const ct::Tensor& psi_k,
                                                                 ct::Tensor& ekb,
                                                                 std::ofstream& ofs_running);
#endif // __CUDA
#endif // __MPI

} // namespace module_rt
