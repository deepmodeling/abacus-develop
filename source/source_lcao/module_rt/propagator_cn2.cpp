#include "propagator.h"
#include "source_base/parallel_blkcyc.h"
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_container/ATen/kernels/lapack.h"
#include "source_base/module_container/ATen/kernels/memory.h" // memory operations (Tensor)
#include "source_base/module_device/memory_op.h"              // memory operations
#include "source_base/matrix2d_linalg.h"

#include <cassert>
#include <complex>
#include <cstddef>
#include <iostream>
#include <vector>

namespace module_rt
{
#ifdef __MPI
void Propagator::compute_propagator_cn2(
                                        Matrix2dView<const std::complex<double>> Stmp,
                                        Matrix2dView<const std::complex<double>> Htmp,
                                        Matrix2dView<std::complex<double>> U_operator,
                                        std::ofstream& ofs_running,
                                        const int print_matrix) const
{
    assert(Htmp.size() > 0);

    const int nrloc = Htmp.nrloc();
    const int ncloc = Htmp.ncloc();
    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double> beta1(0.0, -0.25 * this->dt);
    const std::complex<double> beta2(0.0, 0.25 * this->dt);
    matrix2d_linalg::linear_combination(
        U_operator, alpha, Stmp, beta1, Htmp);

    std::vector<std::complex<double>> denominator(Htmp.size());
    Matrix2dView<std::complex<double>> denominator_view(
        denominator.data(), nrloc, ncloc, Htmp.desc());
    matrix2d_linalg::linear_combination(
        denominator_view, alpha, Stmp, beta2, Htmp);

    if (print_matrix)
    {
        ofs_running << std::endl;
        ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < nrloc; ++i)
        {
            for (int j = 0; j < ncloc; ++j)
            {
                const std::complex<double>& value = Stmp.at_local(i, j);
                ofs_running << value.real() << "+" << value.imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < nrloc; ++i)
        {
            for (int j = 0; j < ncloc; ++j)
            {
                const std::complex<double>& value = Htmp.at_local(i, j);
                ofs_running << value.real() << "+" << value.imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (print_matrix)
    {
        ofs_running << " beta=" << beta1 << std::endl;
        ofs_running << " fenmu:" << std::endl;
        for (int i = 0; i < nrloc; ++i)
        {
            for (int j = 0; j < ncloc; ++j)
            {
                const std::complex<double>& value = denominator_view.at_local(i, j);
                ofs_running << value.real() << "+" << value.imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " fenzi:" << std::endl;
        for (int i = 0; i < nrloc; ++i)
        {
            for (int j = 0; j < ncloc; ++j)
            {
                const std::complex<double>& value = U_operator.at_local(i, j);
                ofs_running << value.real() << "+" << value.imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    matrix2d_linalg::solve(denominator_view, U_operator);

    if (print_matrix)
    {
        ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < nrloc; ++i)
        {
            for (int j = 0; j < ncloc; ++j)
            {
                double aa = U_operator.at_local(i, j).real();
                double bb = U_operator.at_local(i, j).imag();
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
    }

}

void Propagator::compute_propagator_cn2_tensor(
                                               const Parallel_Blkcyclic& layout,
                                               const int nlocal,
                                               const ct::Tensor& Stmp,
                                               const ct::Tensor& Htmp,
                                               ct::Tensor& U_operator,
                                               std::ofstream& ofs_running,
                                               const int print_matrix,
                                               CublasMpResources& cublas_res) const
{
#ifdef __CUBLASMP
    // 1. Resource Validation
    if (!cublas_res.is_initialized || cublas_res.cublasmp_grid == nullptr || cublas_res.cusolvermp_grid == nullptr)
    {
        return;
    }

    assert(Stmp.device_type() == ct::DeviceType::GpuDevice);
    assert(Htmp.device_type() == ct::DeviceType::GpuDevice);
    assert(U_operator.device_type() == ct::DeviceType::GpuDevice);

    // 2. Extract Pointers
    void* d_S = static_cast<void*>(Stmp.data<std::complex<double>>());
    void* d_H = static_cast<void*>(Htmp.data<std::complex<double>>());
    void* d_Num = static_cast<void*>(U_operator.data<std::complex<double>>());

    const std::size_t len_loc = layout.size();

    // Allocate temporary tensor for denominator matrix
    ct::Tensor Denominator_gpu(ct::DataType::DT_COMPLEX_DOUBLE,
                               ct::DeviceType::GpuDevice,
                               ct::TensorShape({static_cast<int64_t>(len_loc)}));
    void* d_Den = static_cast<void*>(Denominator_gpu.data<std::complex<double>>());

    // 3. Matrix Descriptors Creation
    int64_t m_global = layout.desc()[2];
    int64_t n_global = layout.desc()[3];
    int64_t mb = layout.desc()[4];
    int64_t nb = layout.desc()[5];
    int64_t rsrc = layout.desc()[6];
    int64_t csrc = layout.desc()[7];
    int64_t lld = layout.desc()[8];

    // 3.1 cuBLASMp Descriptor
    cublasMpMatrixDescriptor_t desc_blas;
    cublasMpMatrixDescriptorCreate(m_global,
                                   n_global,
                                   mb,
                                   nb,
                                   rsrc,
                                   csrc,
                                   lld,
                                   CUDA_C_64F,
                                   cublas_res.cublasmp_grid,
                                   &desc_blas);

    // 3.2 cuSOLVERMp Descriptor
    cusolverMpMatrixDescriptor_t desc_solver;
    cusolverMpCreateMatrixDesc(&desc_solver,
                               cublas_res.cusolvermp_grid,
                               CUDA_C_64F,
                               m_global,
                               n_global,
                               mb,
                               nb,
                               rsrc,
                               csrc,
                               lld);

    // 4. Construct A (Denominator) and B (Numerator) using Geadd
    std::complex<double> one = {1.0, 0.0};
    std::complex<double> coef_neg_i = {0.0, -0.25 * this->dt};
    std::complex<double> coef_pos_i = {0.0, 0.25 * this->dt};

    cudaMemcpyAsync(d_Num, d_S, len_loc * sizeof(std::complex<double>), cudaMemcpyDeviceToDevice, cublas_res.stream);
    cudaMemcpyAsync(d_Den, d_S, len_loc * sizeof(std::complex<double>), cudaMemcpyDeviceToDevice, cublas_res.stream);

    size_t ws_geadd_dev = 0, ws_geadd_host = 0;
    cublasMpGeadd_bufferSize(cublas_res.cublasmp_handle,
                             CUBLAS_OP_N,
                             m_global,
                             n_global,
                             &coef_neg_i,
                             d_H,
                             1,
                             1,
                             desc_blas,
                             &one,
                             d_Num,
                             1,
                             1,
                             desc_blas,
                             &ws_geadd_dev,
                             &ws_geadd_host);

    void *d_work_geadd = nullptr, *h_work_geadd = nullptr;
    cudaMallocAsync(&d_work_geadd, ws_geadd_dev, cublas_res.stream);
    h_work_geadd = malloc(ws_geadd_host);

    // B = S - i * (dt/4) * H
    cublasMpGeadd(cublas_res.cublasmp_handle,
                  CUBLAS_OP_N,
                  m_global,
                  n_global,
                  &coef_neg_i,
                  d_H,
                  1,
                  1,
                  desc_blas,
                  &one,
                  d_Num,
                  1,
                  1,
                  desc_blas,
                  d_work_geadd,
                  ws_geadd_dev,
                  h_work_geadd,
                  ws_geadd_host);

    // A = S + i * (dt/4) * H
    cublasMpGeadd(cublas_res.cublasmp_handle,
                  CUBLAS_OP_N,
                  m_global,
                  n_global,
                  &coef_pos_i,
                  d_H,
                  1,
                  1,
                  desc_blas,
                  &one,
                  d_Den,
                  1,
                  1,
                  desc_blas,
                  d_work_geadd,
                  ws_geadd_dev,
                  h_work_geadd,
                  ws_geadd_host);

    cudaFreeAsync(d_work_geadd, cublas_res.stream);
    free(h_work_geadd);

    // 5. QR Factorization of A (Denominator)
    int64_t tau_size = m_global + nb;
    void* d_tau = nullptr;
    cudaMallocAsync(&d_tau, tau_size * sizeof(std::complex<double>), cublas_res.stream);

    int* d_info = nullptr;
    cudaMallocAsync(&d_info, sizeof(int), cublas_res.stream);
    cudaMemsetAsync(d_info, 0, sizeof(int), cublas_res.stream);

    size_t ws_geqrf_dev = 0, ws_geqrf_host = 0;
    cusolverMpGeqrf_bufferSize(cublas_res.cusolvermp_handle,
                               m_global,
                               n_global,
                               d_Den,
                               1,
                               1,
                               desc_solver,
                               CUDA_C_64F,
                               &ws_geqrf_dev,
                               &ws_geqrf_host);

    void *d_work_geqrf = nullptr, *h_work_geqrf = nullptr;
    cudaMallocAsync(&d_work_geqrf, ws_geqrf_dev, cublas_res.stream);
    h_work_geqrf = malloc(ws_geqrf_host);

    cusolverMpGeqrf(cublas_res.cusolvermp_handle,
                    m_global,
                    n_global,
                    d_Den,
                    1,
                    1,
                    desc_solver,
                    d_tau,
                    CUDA_C_64F,
                    d_work_geqrf,
                    ws_geqrf_dev,
                    h_work_geqrf,
                    ws_geqrf_host,
                    d_info);

    cudaFreeAsync(d_work_geqrf, cublas_res.stream);
    free(h_work_geqrf);

    // Check QR Info
    int h_info = 0;
    cudaMemcpyAsync(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost, cublas_res.stream);
    cudaStreamSynchronize(cublas_res.stream);
    if (h_info != 0)
    {
        std::cerr << "CRITICAL: cusolverMpGeqrf failed with Info: " << h_info << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // 6. Apply Q^H to B (Numerator)
    size_t ws_ormqr_dev = 0, ws_ormqr_host = 0;
    cusolverMpOrmqr_bufferSize(cublas_res.cusolvermp_handle,
                               CUBLAS_SIDE_LEFT,
                               CUBLAS_OP_C,
                               m_global,
                               n_global,
                               n_global,
                               d_Den,
                               1,
                               1,
                               desc_solver,
                               d_tau,
                               d_Num,
                               1,
                               1,
                               desc_solver,
                               CUDA_C_64F,
                               &ws_ormqr_dev,
                               &ws_ormqr_host);

    void *d_work_ormqr = nullptr, *h_work_ormqr = nullptr;
    cudaMallocAsync(&d_work_ormqr, ws_ormqr_dev, cublas_res.stream);
    h_work_ormqr = malloc(ws_ormqr_host);

    cusolverMpOrmqr(cublas_res.cusolvermp_handle,
                    CUBLAS_SIDE_LEFT,
                    CUBLAS_OP_C,
                    m_global,
                    n_global,
                    n_global,
                    d_Den,
                    1,
                    1,
                    desc_solver,
                    d_tau,
                    d_Num,
                    1,
                    1,
                    desc_solver,
                    CUDA_C_64F,
                    d_work_ormqr,
                    ws_ormqr_dev,
                    h_work_ormqr,
                    ws_ormqr_host,
                    d_info);

    cudaFreeAsync(d_work_ormqr, cublas_res.stream);
    free(h_work_ormqr);

    // 7. Solve Triangular System (TRSM)
    size_t ws_trsm_dev = 0, ws_trsm_host = 0;
    std::complex<double> alpha_trsm = {1.0, 0.0};

    cublasMpTrsm_bufferSize(cublas_res.cublasmp_handle,
                            CUBLAS_SIDE_LEFT,
                            CUBLAS_FILL_MODE_UPPER,
                            CUBLAS_OP_N,
                            CUBLAS_DIAG_NON_UNIT,
                            m_global,
                            n_global,
                            &alpha_trsm,
                            d_Den,
                            1,
                            1,
                            desc_blas,
                            d_Num,
                            1,
                            1,
                            desc_blas,
                            CUBLAS_COMPUTE_64F,
                            &ws_trsm_dev,
                            &ws_trsm_host);

    void *d_work_trsm = nullptr, *h_work_trsm = nullptr;
    cudaMallocAsync(&d_work_trsm, ws_trsm_dev, cublas_res.stream);
    h_work_trsm = malloc(ws_trsm_host);

    cublasMpTrsm(cublas_res.cublasmp_handle,
                 CUBLAS_SIDE_LEFT,
                 CUBLAS_FILL_MODE_UPPER,
                 CUBLAS_OP_N,
                 CUBLAS_DIAG_NON_UNIT,
                 m_global,
                 n_global,
                 &alpha_trsm,
                 d_Den,
                 1,
                 1,
                 desc_blas,
                 d_Num,
                 1,
                 1,
                 desc_blas,
                 CUBLAS_COMPUTE_64F,
                 d_work_trsm,
                 ws_trsm_dev,
                 h_work_trsm,
                 ws_trsm_host);

    cudaFreeAsync(d_work_trsm, cublas_res.stream);
    free(h_work_trsm);

    // 8. Cleanup and Final Synchronization
    cudaStreamSynchronize(cublas_res.stream);

    cublasMpMatrixDescriptorDestroy(desc_blas);
    cusolverMpDestroyMatrixDesc(desc_solver);

    cudaFreeAsync(d_tau, cublas_res.stream);
    cudaFreeAsync(d_info, cublas_res.stream);
#endif // __CUBLASMP
}

template <typename Device>
void Propagator::compute_propagator_cn2_tensor_lapack(const int nlocal,
                                                      const ct::Tensor& Stmp,
                                                      const ct::Tensor& Htmp,
                                                      ct::Tensor& U_operator,
                                                      std::ofstream& ofs_running,
                                                      const int print_matrix) const
{
    // ct_device_type = ct::DeviceType::CpuDevice or ct::DeviceType::GpuDevice
    ct::DeviceType ct_device_type = ct::DeviceTypeToEnum<Device>::value;
    // ct_Device = ct::DEVICE_CPU or ct::DEVICE_GPU
    using ct_Device = typename ct::PsiToContainer<Device>::type;

    // Define coefficients
    // beta1 = -i * dt/4 (for Numerator)
    // beta2 = +i * dt/4 (for Denominator)
    std::complex<double> beta1 = {0.0, -0.25 * this->dt};
    std::complex<double> beta2 = {0.0, 0.25 * this->dt};

    // ========================================================================
    // Numerator = Stmp + beta1 * Htmp
    // ========================================================================

    // 1. Copy Stmp to U_operator
    base_device::memory::synchronize_memory_op<std::complex<double>, Device, Device>()(
        U_operator.data<std::complex<double>>(),
        Stmp.data<std::complex<double>>(),
        nlocal * nlocal);

    // 2. U_operator = beta1 * Htmp + U_operator
    ct::kernels::blas_axpy<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                              &beta1,
                                                              Htmp.data<std::complex<double>>(),
                                                              1,
                                                              U_operator.data<std::complex<double>>(),
                                                              1);

    // ========================================================================
    // Denominator = Stmp + beta2 * Htmp
    // ========================================================================

    ct::Tensor Denominator(ct::DataType::DT_COMPLEX_DOUBLE, ct_device_type, ct::TensorShape({nlocal * nlocal}));

    // 1. Copy Stmp to Denominator
    base_device::memory::synchronize_memory_op<std::complex<double>, Device, Device>()(
        Denominator.data<std::complex<double>>(),
        Stmp.data<std::complex<double>>(),
        nlocal * nlocal);

    // 2. Denominator = beta2 * Htmp + Denominator
    ct::kernels::blas_axpy<std::complex<double>, ct_Device>()(nlocal * nlocal,
                                                              &beta2,
                                                              Htmp.data<std::complex<double>>(),
                                                              1,
                                                              Denominator.data<std::complex<double>>(),
                                                              1);

    // ========================================================================
    // Solve D * U = N, result overwrites N (which is U_operator)
    // ========================================================================

    ct::Tensor ipiv(ct::DataType::DT_INT, ct_device_type, ct::TensorShape({nlocal}));
    // No need to zero ipiv, it is output only.

    // 1. LU Factorization of Denominator (In-place)
    ct::kernels::lapack_getrf<std::complex<double>, ct_Device>()(nlocal,
                                                                 nlocal,
                                                                 Denominator.data<std::complex<double>>(),
                                                                 nlocal,
                                                                 ipiv.data<int>());

    // 2. Solve D * X = B
    ct::kernels::lapack_getrs<std::complex<double>, ct_Device>()('N',
                                                                 nlocal,
                                                                 nlocal,
                                                                 Denominator.data<std::complex<double>>(),
                                                                 nlocal,
                                                                 ipiv.data<int>(),
                                                                 U_operator.data<std::complex<double>>(),
                                                                 nlocal);
}

// Explicit instantiation of template functions
template void Propagator::compute_propagator_cn2_tensor_lapack<base_device::DEVICE_CPU>(const int nlocal,
                                                                                        const ct::Tensor& Stmp,
                                                                                        const ct::Tensor& Htmp,
                                                                                        ct::Tensor& U_operator,
                                                                                        std::ofstream& ofs_running,
                                                                                        const int print_matrix) const;
#if ((defined __CUDA) /* || (defined __ROCM) */)
template void Propagator::compute_propagator_cn2_tensor_lapack<base_device::DEVICE_GPU>(const int nlocal,
                                                                                        const ct::Tensor& Stmp,
                                                                                        const ct::Tensor& Htmp,
                                                                                        ct::Tensor& U_operator,
                                                                                        std::ofstream& ofs_running,
                                                                                        const int print_matrix) const;
#endif // __CUDA
#endif // __MPI
} // namespace module_rt
