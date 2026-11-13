#include "middle_hamilt.h"

#include "source_base/parallel_blkcyc.h"
#include "source_base/global_variable.h"
#include "source_base/matrix2d_linalg.h"
#include "source_base/matrix2d_view.h"
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_device/memory_op.h" // memory operations
#include "source_base/timer.h"

#include <cassert>
#include <complex>
#include <iomanip>
#include <iostream>

namespace module_rt
{
#ifdef __MPI
namespace
{
void print_local_matrix(
    std::ofstream& ofs_running,
    const Matrix2dView<const std::complex<double>>& matrix,
    const char* label)
{
    ofs_running << std::setprecision(10);
    ofs_running << std::endl;
    ofs_running << " " << label << " :" << std::endl;
    for (int i = 0; i < matrix.nrloc(); ++i)
    {
        for (int j = 0; j < matrix.ncloc(); ++j)
        {
            const std::complex<double>& value = matrix.at_local(i, j);
            ofs_running << value.real() << "+" << value.imag() << "i ";
        }
        ofs_running << std::endl;
    }
    ofs_running << std::endl;
}
} // namespace

void half_Hmatrix(
                  Matrix2dView<const std::complex<double>> previous_hamiltonian,
                  Matrix2dView<const std::complex<double>> previous_overlap,
                  Matrix2dView<std::complex<double>> hamiltonian,
                  Matrix2dView<std::complex<double>> overlap,
                  std::ofstream& ofs_running,
                  const int print_matrix)
{
    if (print_matrix)
    {
        print_local_matrix(ofs_running, hamiltonian, "H(t+dt)");
        print_local_matrix(ofs_running, previous_hamiltonian, "H(t)");
    }

    const std::complex<double> half(0.5, 0.0);
    matrix2d_linalg::linear_combination(
        hamiltonian, half, previous_hamiltonian, half, hamiltonian);
    matrix2d_linalg::linear_combination(
        overlap, half, previous_overlap, half, overlap);

    if (print_matrix)
    {
        print_local_matrix(ofs_running, hamiltonian, "H(t+dt/2)");
    }
}

void half_Hmatrix_tensor(const Parallel_Blkcyclic& layout,
                         const int nband,
                         const int nlocal,
                         ct::Tensor& Htmp,
                         ct::Tensor& Stmp,
                         const ct::Tensor& H_laststep,
                         const ct::Tensor& S_laststep,
                         std::ofstream& ofs_running,
                         const int print_matrix,
                         CublasMpResources& cublas_res)
{
#ifdef __CUBLASMP
    // 1. Validate resources and ensure the grid is properly initialized
    if (!cublas_res.is_initialized || cublas_res.cublasmp_grid == nullptr)
    {
        return;
    }

    assert(Htmp.device_type() == ct::DeviceType::GpuDevice);
    assert(Stmp.device_type() == ct::DeviceType::GpuDevice);
    assert(H_laststep.device_type() == ct::DeviceType::GpuDevice);
    assert(S_laststep.device_type() == ct::DeviceType::GpuDevice);

    // 2. Extract device pointers
    void* d_Htmp = static_cast<void*>(Htmp.data<std::complex<double>>());
    void* d_Stmp = static_cast<void*>(Stmp.data<std::complex<double>>());
    void* d_H_last = static_cast<void*>(const_cast<std::complex<double>*>(H_laststep.data<std::complex<double>>()));
    void* d_S_last = static_cast<void*>(const_cast<std::complex<double>*>(S_laststep.data<std::complex<double>>()));

    int64_t m_global = layout.desc()[2];
    int64_t n_global = layout.desc()[3];
    int64_t mb = layout.desc()[4];
    int64_t nb = layout.desc()[5];
    int64_t rsrc = layout.desc()[6];
    int64_t csrc = layout.desc()[7];
    int64_t lld = layout.desc()[8];

    // 3. Create matrix descriptor
    cublasMpMatrixDescriptor_t desc_mat;
    cublasMpMatrixDescriptorCreate(m_global,
                                   n_global,
                                   mb,
                                   nb,
                                   rsrc,
                                   csrc,
                                   lld,
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_mat);

    std::complex<double> alpha = {0.5, 0.0};
    std::complex<double> beta = {0.5, 0.0};

    size_t ws_size_dev = 0;
    size_t ws_size_host = 0;

    // 4. Query workspace size
    cublasMpGeadd_bufferSize(cublas_res.cublasmp_handle,
                             CUBLAS_OP_N,
                             m_global,
                             n_global,
                             &alpha,
                             d_H_last,
                             1,
                             1,
                             desc_mat,
                             &beta,
                             d_Htmp,
                             1,
                             1,
                             desc_mat,
                             &ws_size_dev,
                             &ws_size_host);

    void* d_work = nullptr;
    void* h_work = nullptr;

    cudaMallocAsync(&d_work, ws_size_dev, cublas_res.stream);
    h_work = malloc(ws_size_host);

    // 5. Compute Htmp = 0.5 * H_last + 0.5 * Htmp
    cublasMpGeadd(cublas_res.cublasmp_handle,
                  CUBLAS_OP_N,
                  m_global,
                  n_global,
                  &alpha,
                  d_H_last,
                  1,
                  1,
                  desc_mat,
                  &beta,
                  d_Htmp,
                  1,
                  1,
                  desc_mat,
                  d_work,
                  ws_size_dev,
                  h_work,
                  ws_size_host);

    // 6. Compute Stmp = 0.5 * S_last + 0.5 * Stmp
    cublasMpGeadd(cublas_res.cublasmp_handle,
                  CUBLAS_OP_N,
                  m_global,
                  n_global,
                  &alpha,
                  d_S_last,
                  1,
                  1,
                  desc_mat,
                  &beta,
                  d_Stmp,
                  1,
                  1,
                  desc_mat,
                  d_work,
                  ws_size_dev,
                  h_work,
                  ws_size_host);

    // 7. Synchronize stream and release resources
    cudaStreamSynchronize(cublas_res.stream);

    cublasMpMatrixDescriptorDestroy(desc_mat);
    cudaFreeAsync(d_work, cublas_res.stream);
    free(h_work);
#endif // __CUBLASMP
}

template <typename Device>
void half_Hmatrix_tensor_lapack(const int nband,
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
}

// Explicit instantiation of template functions
template void half_Hmatrix_tensor_lapack<base_device::DEVICE_CPU>(const int nband,
                                                                  const int nlocal,
                                                                  ct::Tensor& Htmp,
                                                                  ct::Tensor& Stmp,
                                                                  const ct::Tensor& H_laststep,
                                                                  const ct::Tensor& S_laststep,
                                                                  std::ofstream& ofs_running,
                                                                  const int print_matrix);
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void half_Hmatrix_tensor_lapack<base_device::DEVICE_GPU>(const int nband,
                                                                  const int nlocal,
                                                                  ct::Tensor& Htmp,
                                                                  ct::Tensor& Stmp,
                                                                  const ct::Tensor& H_laststep,
                                                                  const ct::Tensor& S_laststep,
                                                                  std::ofstream& ofs_running,
                                                                  const int print_matrix);
#endif // __CUDA
#endif // __MPI
} // namespace module_rt
