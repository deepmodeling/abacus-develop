#include "propagator.h"
#include "source_base/global_function.h"
#include "source_base/module_container/ATen/kernels/blas.h"
#include "source_base/module_container/ATen/kernels/lapack.h"
#include "source_base/module_container/ATen/kernels/memory.h" // memory operations (Tensor)
#include "source_base/module_device/memory_op.h"              // memory operations
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/timer.h"
#include "source_io/module_parameter/parameter.h"

#include <cassert>
#include <complex>
#include <iostream>

namespace module_rt
{
#ifdef __MPI
void Propagator::compute_propagator_cn2(const int nlocal,
                                        const std::complex<double>* Stmp,
                                        const std::complex<double>* Htmp,
                                        std::complex<double>* U_operator,
                                        std::ofstream& ofs_running,
                                        const int print_matrix) const
{
    assert(this->ParaV->nloc > 0);

    // (1) copy Htmp to Numerator & Denominator
    std::complex<double>* Numerator = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Numerator, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Htmp, 1, Numerator, 1);

    std::complex<double>* Denominator = new std::complex<double>[this->ParaV->nloc];
    ModuleBase::GlobalFunc::ZEROS(Denominator, this->ParaV->nloc);
    BlasConnector::copy(this->ParaV->nloc, Htmp, 1, Denominator, 1);

    if (print_matrix)
    {
        ofs_running << std::endl;
        ofs_running << " S matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Stmp[in + j].real() << "+" << Stmp[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << std::endl;
        ofs_running << " H matrix :" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Numerator[in + j].real() << "+" << Numerator[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    // ->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (2) compute Numerator & Denominator by GEADD
    // Numerator = Stmp - i*para * Htmp;     beta1 = - para = -0.25 * this->dt
    // Denominator = Stmp + i*para * Htmp;   beta2 = para = 0.25 * this->dt
    std::complex<double> alpha = {1.0, 0.0};
    std::complex<double> beta1 = {0.0, -0.25 * this->dt};
    std::complex<double> beta2 = {0.0, 0.25 * this->dt};

    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              Stmp,
                              1,
                              1,
                              this->ParaV->desc,
                              beta1,
                              Numerator,
                              1,
                              1,
                              this->ParaV->desc);
    ScalapackConnector::geadd('N',
                              nlocal,
                              nlocal,
                              alpha,
                              Stmp,
                              1,
                              1,
                              this->ParaV->desc,
                              beta2,
                              Denominator,
                              1,
                              1,
                              this->ParaV->desc);

    if (print_matrix)
    {
        ofs_running << " beta=" << beta1 << std::endl;
        ofs_running << " fenmu:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Denominator[in + j].real() << "+" << Denominator[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
    }

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // (3) Next, invert Denominator
    // What is the size of ipiv exactly? Need to check ScaLAPACK documentation!
    // But anyway, not this->ParaV->nloc
    int* ipiv = new int[this->ParaV->nrow + this->ParaV->nb];
    ModuleBase::GlobalFunc::ZEROS(ipiv, this->ParaV->nrow + this->ParaV->nb);
    int info = 0;
    // (3.1) compute ipiv
    ScalapackConnector::getrf(nlocal, nlocal, Denominator, 1, 1, this->ParaV->desc, ipiv, &info);

    // Print ipiv
    if (print_matrix)
    {
        ofs_running << " this->ParaV->nloc = " << this->ParaV->nloc << std::endl;
        ofs_running << " this->ParaV->nrow = " << this->ParaV->nrow << std::endl;
        ofs_running << " this->ParaV->ncol = " << this->ParaV->ncol << std::endl;
        ofs_running << " this->ParaV->nb = " << this->ParaV->nb << std::endl;
        ofs_running << " this->ParaV->get_block_size() = " << this->ParaV->get_block_size() << std::endl;
        ofs_running << " nlocal = " << nlocal << std::endl;
        ofs_running << " ipiv:" << std::endl;
        for (int i = 0; i < this->ParaV->nloc; i++)
        {
            ofs_running << ipiv[i] << " ";
        }
        ofs_running << std::endl;
    }

    int lwork = -1;
    int liwotk = -1;
    std::vector<std::complex<double>> work(1, 0);
    std::vector<int> iwork(1, 0);
    // (3.2) compute work
    ScalapackConnector::getri(nlocal,
                              Denominator,
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv,
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwotk,
                              &info);
    lwork = work[0].real();
    work.resize(lwork, 0);
    liwotk = iwork[0];
    iwork.resize(liwotk, 0);
    // (3.3) compute inverse matrix of Denominator
    ScalapackConnector::getri(nlocal,
                              Denominator,
                              1,
                              1,
                              this->ParaV->desc,
                              ipiv,
                              work.data(),
                              &lwork,
                              iwork.data(),
                              &liwotk,
                              &info);
    assert(0 == info);

    //->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    // (4) U_operator = Denominator * Numerator;
    ScalapackConnector::gemm('N',
                             'N',
                             nlocal,
                             nlocal,
                             nlocal,
                             1.0,
                             Denominator,
                             1,
                             1,
                             this->ParaV->desc,
                             Numerator,
                             1,
                             1,
                             this->ParaV->desc,
                             0.0,
                             U_operator,
                             1,
                             1,
                             this->ParaV->desc);

    if (print_matrix)
    {
        ofs_running << " fenmu^-1:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Denominator[in + j].real() << "+" << Denominator[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " fenzi:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                ofs_running << Numerator[in + j].real() << "+" << Numerator[in + j].imag() << "i ";
            }
            ofs_running << std::endl;
        }
        ofs_running << std::endl;
        ofs_running << " U operator:" << std::endl;
        for (int i = 0; i < this->ParaV->nrow; i++)
        {
            const int in = i * this->ParaV->ncol;
            for (int j = 0; j < this->ParaV->ncol; j++)
            {
                double aa = U_operator[in + j].real();
                double bb = U_operator[in + j].imag();
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

    delete[] Numerator;
    delete[] Denominator;
    delete[] ipiv;
}

void Propagator::compute_propagator_cn2_tensor(const int nlocal,
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

    int64_t len_loc = this->ParaV->nloc;

    // Allocate temporary tensor for denominator matrix
    ct::Tensor Denominator_gpu(ct::DataType::DT_COMPLEX_DOUBLE, ct::DeviceType::GpuDevice, ct::TensorShape({len_loc}));
    void* d_Den = static_cast<void*>(Denominator_gpu.data<std::complex<double>>());

    // 3. Matrix Descriptors Creation
    int64_t m_global = this->ParaV->desc[2];
    int64_t n_global = this->ParaV->desc[3];
    int64_t mb = this->ParaV->desc[4];
    int64_t nb = this->ParaV->desc[5];
    int64_t rsrc = this->ParaV->desc[6];
    int64_t csrc = this->ParaV->desc[7];
    int64_t lld = this->ParaV->desc[8];

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
