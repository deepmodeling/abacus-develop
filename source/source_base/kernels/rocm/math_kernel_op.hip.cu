#include "source_base/module_device/memory_op.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_psi/psi.h"
#include "source_base/tool_quit.h"

#include <base/macros/macros.h>
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <thrust/complex.h>
template <>
struct GetTypeReal<thrust::complex<float>> {
    using type = float; /**< The return type specialization for std::complex<double>. */
};
template <>
struct GetTypeReal<thrust::complex<double>> {
    using type = double; /**< The return type specialization for std::complex<double>. */
};

namespace ModuleBase {

template <typename T>
struct GetTypeThrust {
    using type = T;
};

template <>
struct GetTypeThrust<std::complex<float>> {
    using type = thrust::complex<float>; /**< The return type specialization for std::complex<float>. */
};

template <>
struct GetTypeThrust<std::complex<double>> {
    using type = thrust::complex<double>; /**< The return type specialization for std::complex<float>. */
};

static hipblasHandle_t cublas_handle = nullptr;

void xdot_wrapper(const int &n, const float * x, const int &incx, const float * y, const int &incy, float &result) {
    hipblasErrcheck(hipblasSdot(cublas_handle, n, x, incx, y, incy, &result));
}

void xdot_wrapper(const int &n, const double * x, const int &incx, const double * y, const int &incy, double &result) {
    hipblasErrcheck(hipblasDdot(cublas_handle, n, x, incx, y, incy, &result));
}

void createGpuBlasHandle(){
    if (cublas_handle == nullptr) {
        hipblasErrcheck(hipblasCreate(&cublas_handle));
    }
}

void destoryBLAShandle(){
    if (cublas_handle != nullptr) {
        hipblasErrcheck(hipblasDestroy(cublas_handle));
        cublas_handle = nullptr;
    }
}

template <>
void scal_op<float, base_device::DEVICE_GPU>::operator()(const int& N,
                                                         const std::complex<float>* alpha,
                                                         std::complex<float>* X,
                                                         const int& incx)
{
    hipblasErrcheck(hipblasCscal(cublas_handle, N, (hipblasComplex*)alpha, (hipblasComplex*)X, incx));
}

template <>
void scal_op<double, base_device::DEVICE_GPU>::operator()(const int& N,
                                                          const std::complex<double>* alpha,
                                                          std::complex<double>* X,
                                                          const int& incx)
{
    hipblasErrcheck(hipblasZscal(cublas_handle, N, (hipblasDoubleComplex*)alpha, (hipblasDoubleComplex*)X, incx));
}

template <>
void axpy_op<double, base_device::DEVICE_GPU>::operator()(const int& N,
                                                          const double* alpha,
                                                          const double* X,
                                                          const int& incX,
                                                          double* Y,
                                                          const int& incY)
{
    hipblasErrcheck(hipblasDaxpy(cublas_handle, N, alpha, X, incX, Y, incY));
}

template <>
void axpy_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const int& N,
                                                                       const std::complex<float>* alpha,
                                                                       const std::complex<float>* X,
                                                                       const int& incX,
                                                                       std::complex<float>* Y,
                                                                       const int& incY)
{
    hipblasErrcheck(
        hipblasCaxpy(cublas_handle, N, (hipblasComplex*)alpha, (hipblasComplex*)X, incX, (hipblasComplex*)Y, incY));
}

template <>
void axpy_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const int& N,
                                                                        const std::complex<double>* alpha,
                                                                        const std::complex<double>* X,
                                                                        const int& incX,
                                                                        std::complex<double>* Y,
                                                                        const int& incY)
{
    hipblasErrcheck(hipblasZaxpy(cublas_handle,
                                 N,
                                 (hipblasDoubleComplex*)alpha,
                                 (hipblasDoubleComplex*)X,
                                 incX,
                                 (hipblasDoubleComplex*)Y,
                                 incY));
}

template <typename T>
__launch_bounds__(1024)
__global__ void matrix_transpose_kernel(
        const int row,
        const int col,
    const T* in,
    T* out)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < row)
    {
        for (int j = 0; j < col; j++)
        {
            out[j * row + i] = in[i * col + j];
        }
    }
}

template <typename T>
__launch_bounds__(1024) __global__
    void matrix_copy_kernel(const int n1, const int n2, const T* A, const int LDA, T* B, const int LDB)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < n1 && j < n2)
    {
        B[i * LDB + j] = A[i * LDA + j];
    }
}

hipblasOperation_t judge_trans_op(bool is_complex, const char& trans, const char* name)
{
    if (trans == 'N')
    {
        return HIPBLAS_OP_N;
    }
    else if(trans == 'T')
    {
        return HIPBLAS_OP_T;
    }
    else if(is_complex && trans == 'C')
    {
        return HIPBLAS_OP_C;
    }
    else 
    {
        ModuleBase::WARNING_QUIT(name, std::string("Unknown trans type ") + trans + std::string(" !"));
    }
}

template <>
void gemv_op<double, base_device::DEVICE_GPU>::operator()(const char& trans,
                                                          const int& m,
                                                          const int& n,
                                                          const double* alpha,
                                                          const double* A,
                                                          const int& lda,
                                                          const double* X,
                                                          const int& incx,
                                                          const double* beta,
                                                          double* Y,
                                                          const int& incy)
{
    hipblasOperation_t cutrans = judge_trans_op(false, trans, "gemv_op");
    hipblasErrcheck(hipblasDgemv(cublas_handle, cutrans, m, n, alpha, A, lda, X, incx, beta, Y, incx));
}

template <>
void gemv_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const char& trans,
                                                                       const int& m,
                                                                       const int& n,
                                                                       const std::complex<float>* alpha,
                                                                       const std::complex<float>* A,
                                                                       const int& lda,
                                                                       const std::complex<float>* X,
                                                                       const int& incx,
                                                                       const std::complex<float>* beta,
                                                                       std::complex<float>* Y,
                                                                       const int& incy)
{
    hipblasOperation_t cutrans = judge_trans_op(true, trans, "gemv_op");
    hipblasErrcheck(hipblasCgemv(cublas_handle, cutrans, m, n, (hipblasComplex*)alpha, (hipblasComplex*)A, lda, (hipblasComplex*)X, incx, (hipblasComplex*)beta, (hipblasComplex*)Y, incx));
}

template <>
void gemv_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const char& trans,
                                                                        const int& m,
                                                                        const int& n,
                                                                        const std::complex<double>* alpha,
                                                                        const std::complex<double>* A,
                                                                        const int& lda,
                                                                        const std::complex<double>* X,
                                                                        const int& incx,
                                                                        const std::complex<double>* beta,
                                                                        std::complex<double>* Y,
                                                                        const int& incy)
{
    hipblasOperation_t cutrans = judge_trans_op(true, trans, "gemv_op");
    hipblasErrcheck(hipblasZgemv(cublas_handle, cutrans, m, n, (hipblasDoubleComplex*)alpha, (hipblasDoubleComplex*)A, lda, (hipblasDoubleComplex*)X, incx, (hipblasDoubleComplex*)beta, (hipblasDoubleComplex*)Y, incx));
}

template <>
void gemm_op<float, base_device::DEVICE_GPU>::operator()(const char& transa,
                                                         const char& transb,
                                                         const int& m,
                                                         const int& n,
                                                         const int& k,
                                                         const float* alpha,
                                                         const float* a,
                                                         const int& lda,
                                                         const float* b,
                                                         const int& ldb,
                                                         const float* beta,
                                                         float* c,
                                                         const int& ldc)
{
    hipblasOperation_t cutransA = judge_trans_op(false, transa, "gemm_op");
    hipblasOperation_t cutransB = judge_trans_op(false, transb, "gemm_op");
    hipblasErrcheck(hipblasSgemm(cublas_handle, cutransA, cutransB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc));
}

template <>
void gemm_op<double, base_device::DEVICE_GPU>::operator()(const char& transa,
                                                          const char& transb,
                                                          const int& m,
                                                          const int& n,
                                                          const int& k,
                                                          const double* alpha,
                                                          const double* a,
                                                          const int& lda,
                                                          const double* b,
                                                          const int& ldb,
                                                          const double* beta,
                                                          double* c,
                                                          const int& ldc)
{
    hipblasOperation_t cutransA = judge_trans_op(false, transa, "gemm_op");
    hipblasOperation_t cutransB = judge_trans_op(false, transb, "gemm_op");
    hipblasErrcheck(hipblasDgemm(cublas_handle, cutransA, cutransB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc));
}

template <>
void gemm_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const char& transa,
                                                                       const char& transb,
                                                                       const int& m,
                                                                       const int& n,
                                                                       const int& k,
                                                                       const std::complex<float>* alpha,
                                                                       const std::complex<float>* a,
                                                                       const int& lda,
                                                                       const std::complex<float>* b,
                                                                       const int& ldb,
                                                                       const std::complex<float>* beta,
                                                                       std::complex<float>* c,
                                                                       const int& ldc)
{
    hipblasOperation_t cutransA = judge_trans_op(true, transa, "gemm_op");
    hipblasOperation_t cutransB = judge_trans_op(true, transb, "gemm_op");
    hipblasErrcheck(hipblasCgemm(cublas_handle, cutransA, cutransB, m, n ,k, (hipblasComplex*)alpha, (hipblasComplex*)a , lda, (hipblasComplex*)b, ldb, (hipblasComplex*)beta, (hipblasComplex*)c, ldc));
}

template <>
void gemm_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const char& transa,
                                                                        const char& transb,
                                                                        const int& m,
                                                                        const int& n,
                                                                        const int& k,
                                                                        const std::complex<double>* alpha,
                                                                        const std::complex<double>* a,
                                                                        const int& lda,
                                                                        const std::complex<double>* b,
                                                                        const int& ldb,
                                                                        const std::complex<double>* beta,
                                                                        std::complex<double>* c,
                                                                        const int& ldc)
{
    hipblasOperation_t cutransA = judge_trans_op(true, transa, "gemm_op");
    hipblasOperation_t cutransB = judge_trans_op(true, transb, "gemm_op");
    hipblasErrcheck(hipblasZgemm(cublas_handle, cutransA, cutransB, m, n ,k, (hipblasDoubleComplex*)alpha, (hipblasDoubleComplex*)a , lda, (hipblasDoubleComplex*)b, ldb, (hipblasDoubleComplex*)beta, (hipblasDoubleComplex*)c, ldc));
}

template <>
void matrixTranspose_op<double, base_device::DEVICE_GPU>::operator()(const int& row,
                                                                     const int& col,
                                                                     const double* input_matrix,
                                                                     double* output_matrix)
{
    double* device_temp = nullptr;
    base_device::memory::resize_memory_op<double, base_device::DEVICE_GPU>()(device_temp, row * col);

    if (row == col)
    {
        double ONE = 1.0, ZERO = 0.0;
        // use 'geam' API todo transpose.
        hipblasErrcheck(hipblasDgeam(cublas_handle, HIPBLAS_OP_T, HIPBLAS_OP_N, col, row, &ONE, input_matrix, col, &ZERO, input_matrix, col, device_temp, col));
    }
    else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_transpose_kernel<double>), dim3(block), dim3(thread), 0, 0, row, col, input_matrix, device_temp);
        hipCheckOnDebug();
    }

    base_device::memory::synchronize_memory_op<double, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(
        output_matrix,
        device_temp,
        row * col);

    base_device::memory::delete_memory_op<double, base_device::DEVICE_GPU>()(device_temp);
}

template <>
void matrixTranspose_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(
    const int& row,
    const int& col,
    const std::complex<float>* input_matrix,
    std::complex<float>* output_matrix)
{
    std::complex<float>* device_temp = nullptr;
    base_device::memory::resize_memory_op<std::complex<float>, base_device::DEVICE_GPU>()(device_temp, row * col);

    if (row == col)
    {
        float2 ONE, ZERO;
        ONE.x = 1.0;
        ONE.y = 0.0;
        ZERO.x = ZERO.y = 0.0;

        // use 'geam' API todo transpose.
        hipblasErrcheck(hipblasCgeam(cublas_handle, HIPBLAS_OP_T, HIPBLAS_OP_N, col, row,
                                   reinterpret_cast<const hipblasComplex *>(&ONE), (hipblasComplex*)input_matrix, col,
                                   reinterpret_cast<const hipblasComplex *>(&ZERO), (hipblasComplex*)input_matrix, col, (hipblasComplex*)device_temp, col));
    } else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_transpose_kernel<thrust::complex<float>>), dim3(block), dim3(thread), 0, 0, row, col, (thrust::complex<float>*)input_matrix, (thrust::complex<float>*)device_temp);
        hipCheckOnDebug();
    }

    base_device::memory::synchronize_memory_op<std::complex<float>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(
        output_matrix,
        device_temp,
        row * col);

    base_device::memory::delete_memory_op<std::complex<float>, base_device::DEVICE_GPU>()(device_temp);
}

template <>
void matrixTranspose_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(
    const int& row,
    const int& col,
    const std::complex<double>* input_matrix,
    std::complex<double>* output_matrix)
{
    std::complex<double>* device_temp = nullptr;
    base_device::memory::resize_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(device_temp, row * col);

    if (row == col)
    {
        hipblasDoubleComplex ONE{1.0, 0.0}, ZERO{0.0, 0.0};
        // use 'geam' API todo transpose.
        hipblasErrcheck(hipblasZgeam(cublas_handle, HIPBLAS_OP_T, HIPBLAS_OP_N, col, row, &ONE, (hipblasDoubleComplex*)input_matrix, col, &ZERO, (hipblasDoubleComplex*)input_matrix, col, (hipblasDoubleComplex*)device_temp, col));
    } else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_transpose_kernel<thrust::complex<double>>), dim3(block), dim3(thread), 0, 0, row, col, (thrust::complex<double>*)input_matrix, (thrust::complex<double>*)device_temp);
        hipCheckOnDebug();
    }

    base_device::memory::synchronize_memory_op<std::complex<double>,
                                               base_device::DEVICE_GPU,
                                               base_device::DEVICE_GPU>()(output_matrix, device_temp, row * col);

    base_device::memory::delete_memory_op<std::complex<double>, base_device::DEVICE_GPU>()(device_temp);
}

template <>
void matrixCopy<double, base_device::DEVICE_GPU>::operator()(const int& n1,
                                                             const int& n2,
                                                             const double* A,
                                                             const int& LDA,
                                                             double* B,
                                                             const int& LDB)
{
    const dim3 blockSize(16, 16);
    const dim3 gridSize((n1 + blockSize.x - 1) / blockSize.x, (n2 + blockSize.y - 1) / blockSize.y);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_copy_kernel<double>), gridSize, blockSize, 0, 0, n1, n2, A, LDA, B, LDB);
    hipCheckOnDebug();
}
template <>
void matrixCopy<std::complex<float>, base_device::DEVICE_GPU>::operator()(const int& n1,
                                                                          const int& n2,
                                                                          const std::complex<float>* A,
                                                                          const int& LDA,
                                                                          std::complex<float>* B,
                                                                          const int& LDB)
{
    const dim3 blockSize(16, 16);
    const dim3 gridSize((n1 + blockSize.x - 1) / blockSize.x, (n2 + blockSize.y - 1) / blockSize.y);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_copy_kernel<thrust::complex<float>>), gridSize, blockSize, 0, 0, n1, n2, reinterpret_cast<const thrust::complex<float>*>(A), LDA, reinterpret_cast<thrust::complex<float>*>(B), LDB);
    hipCheckOnDebug();
}
template <>
void matrixCopy<std::complex<double>, base_device::DEVICE_GPU>::operator()(const int& n1,
                                                                           const int& n2,
                                                                           const std::complex<double>* A,
                                                                           const int& LDA,
                                                                           std::complex<double>* B,
                                                                           const int& LDB)
{
    const dim3 blockSize(16, 16);
    const dim3 gridSize((n1 + blockSize.x - 1) / blockSize.x, (n2 + blockSize.y - 1) / blockSize.y);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(matrix_copy_kernel<thrust::complex<double>>), gridSize, blockSize, 0, 0, n1, n2, reinterpret_cast<const thrust::complex<double>*>(A), LDA, reinterpret_cast<thrust::complex<double>*>(B), LDB);
    hipCheckOnDebug();
}



// Explicitly instantiate functors for the types of functor registered.
template struct matrixCopy<double, base_device::DEVICE_GPU>;
template struct matrixCopy<std::complex<float>, base_device::DEVICE_GPU>;
template struct matrixCopy<std::complex<double>, base_device::DEVICE_GPU>;
}  // namespace ModuleBase
