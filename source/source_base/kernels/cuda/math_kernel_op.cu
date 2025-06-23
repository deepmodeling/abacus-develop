#include "source_base/module_device/memory_op.h"
#include "source_base/kernels/math_kernel_op.h"
#include "module_psi/psi.h"
#include "source_base/tool_quit.h"

#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <thrust/execution_policy.h>
#include <thrust/inner_product.h>
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

static cublasHandle_t cublas_handle = nullptr;

void xdot_wrapper(const int &n, const float * x, const int &incx, const float * y, const int &incy, float &result) {
    cublasErrcheck(cublasSdot(cublas_handle, n, x, incx, y, incy, &result));
}

void xdot_wrapper(const int &n, const double * x, const int &incx, const double * y, const int &incy, double &result) {
    cublasErrcheck(cublasDdot(cublas_handle, n, x, incx, y, incy, &result));
}

void createGpuBlasHandle(){
    if (cublas_handle == nullptr) {
        cublasErrcheck(cublasCreate(&cublas_handle));
    }
}

void destoryBLAShandle(){
    if (cublas_handle != nullptr) {
        cublasErrcheck(cublasDestroy(cublas_handle));
        cublas_handle = nullptr;
    }
}

// template <typename FPTYPE>
// __forceinline__ __device__ void warp_reduce(FPTYPE& val) {
//     for (int offset = 16; offset > 0; offset >>= 1)
//         val += __shfl_down_sync(full_mask, val, offset);
// }
template <>
void scal_op<float, base_device::DEVICE_GPU>::operator()(const int& N,
                                                         const std::complex<float>* alpha,
                                                         std::complex<float>* X,
                                                         const int& incx)
{
    cublasErrcheck(cublasCscal(cublas_handle, N, (float2*)alpha, (float2*)X, incx));
}

template <>
void scal_op<double, base_device::DEVICE_GPU>::operator()(const int& N,
                                                          const std::complex<double>* alpha,
                                                          std::complex<double>* X,
                                                          const int& incx)
{
    cublasErrcheck(cublasZscal(cublas_handle, N, (double2*)alpha, (double2*)X, incx));
}

template <>
void axpy_op<double, base_device::DEVICE_GPU>::operator()(const int& N,
                                                          const double* alpha,
                                                          const double* X,
                                                          const int& incX,
                                                          double* Y,
                                                          const int& incY)
{
    cublasErrcheck(cublasDaxpy(cublas_handle, N, alpha, X, incX, Y, incY));
}

template <>
void axpy_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const int& N,
                                                                       const std::complex<float>* alpha,
                                                                       const std::complex<float>* X,
                                                                       const int& incX,
                                                                       std::complex<float>* Y,
                                                                       const int& incY)
{
    cublasErrcheck(cublasCaxpy(cublas_handle, N, (float2*)alpha, (float2*)X, incX, (float2*)Y, incY));
}

template <>
void axpy_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const int& N,
                                                                        const std::complex<double>* alpha,
                                                                        const std::complex<double>* X,
                                                                        const int& incX,
                                                                        std::complex<double>* Y,
                                                                        const int& incY)
{
    cublasErrcheck(cublasZaxpy(cublas_handle, N, (double2*)alpha, (double2*)X, incX, (double2*)Y, incY));
}


template <typename T>
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
__global__ void matrix_copy_kernel(const int n1, const int n2, const T* A, const int LDA, T* B, const int LDB)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < n1 && j < n2)
    {
        B[i * LDB + j] = A[i * LDA + j];
    }
}

cublasOperation_t judge_trans_op(bool is_complex, const char& trans, const char* name)
{
    if (trans == 'N')
    {
        return CUBLAS_OP_N;
    }
    else if(trans == 'T')
    {
        return CUBLAS_OP_T;
    }
    else if(is_complex && trans == 'C')
    {
        return CUBLAS_OP_C;
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
    cublasOperation_t cutrans = judge_trans_op(false, trans, "gemv_op");
    cublasErrcheck(cublasDgemv(cublas_handle, cutrans, m, n, alpha, A, lda, X, incx, beta, Y, incx));
}

template <>
void gemv_op<std::complex<float>, base_device::DEVICE_GPU>::operator()(const char& trans,
                                                                       const int& m,
                                                                       const int& n,
                                                                       const std::complex<float>* alpha_in,
                                                                       const std::complex<float>* A,
                                                                       const int& lda,
                                                                       const std::complex<float>* X,
                                                                       const int& incx,
                                                                       const std::complex<float>* beta_in,
                                                                       std::complex<float>* Y,
                                                                       const int& incy)
{
    cublasOperation_t cutrans = judge_trans_op(true, trans, "gemv_op");
    cuFloatComplex alpha = make_cuFloatComplex(alpha_in->real(), alpha_in->imag());
    cuFloatComplex beta = make_cuFloatComplex(beta_in->real(), beta_in->imag());
    cublasErrcheck(cublasCgemv(cublas_handle, cutrans, m, n, &alpha, (cuFloatComplex*)A, lda, (cuFloatComplex*)X, incx, &beta, (cuFloatComplex*)Y, incx));
}

template <>
void gemv_op<std::complex<double>, base_device::DEVICE_GPU>::operator()(const char& trans,
                                                                        const int& m,
                                                                        const int& n,
                                                                        const std::complex<double>* alpha_in,
                                                                        const std::complex<double>* A,
                                                                        const int& lda,
                                                                        const std::complex<double>* X,
                                                                        const int& incx,
                                                                        const std::complex<double>* beta_in,
                                                                        std::complex<double>* Y,
                                                                        const int& incy)
{
    cublasOperation_t cutrans = judge_trans_op(true, trans, "gemv_op");
    cuDoubleComplex alpha = make_cuDoubleComplex(alpha_in->real(), alpha_in->imag());
    cuDoubleComplex beta = make_cuDoubleComplex(beta_in->real(), beta_in->imag());
    // icpc and nvcc have some compatible problems
    // We must use cuDoubleComplex instead of converting std::complex<double>* to cuDoubleComplex*
    cublasErrcheck(cublasZgemv(cublas_handle, cutrans, m, n, &alpha, (cuDoubleComplex*)A, lda, (cuDoubleComplex*)X, incx, &beta, (cuDoubleComplex*)Y, incx));
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
    cublasOperation_t cutransA = judge_trans_op(false, transa, "gemm_op");
    cublasOperation_t cutransB = judge_trans_op(false, transb, "gemm_op");
    cublasErrcheck(cublasSgemm(cublas_handle, cutransA, cutransB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc));
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
    cublasOperation_t cutransA = judge_trans_op(false, transa, "gemm_op");
    cublasOperation_t cutransB = judge_trans_op(false, transb, "gemm_op");
    cublasErrcheck(cublasDgemm(cublas_handle, cutransA, cutransB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc));
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
    cublasOperation_t cutransA = judge_trans_op(true, transa, "gemm_op");
    cublasOperation_t cutransB = judge_trans_op(true, transb, "gemm_op");
    cublasErrcheck(cublasCgemm(cublas_handle, cutransA, cutransB, m, n ,k, (float2*)alpha, (float2*)a , lda, (float2*)b, ldb, (float2*)beta, (float2*)c, ldc));
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
    cublasOperation_t cutransA = judge_trans_op(true, transa, "gemm_op");
    cublasOperation_t cutransB = judge_trans_op(true, transb, "gemm_op");
    cublasErrcheck(cublasZgemm(cublas_handle, cutransA, cutransB, m, n ,k, (double2*)alpha, (double2*)a , lda, (double2*)b, ldb, (double2*)beta, (double2*)c, ldc));
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
        cublasErrcheck(cublasDgeam(cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, col, row, &ONE, input_matrix, col, &ZERO, input_matrix, col, device_temp, col));
    }
    else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        matrix_transpose_kernel<double> <<<block, thread >>> (row, col, input_matrix, device_temp);

        cudaCheckOnDebug();
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
        double2 ONE, ZERO;
        ONE.x = 1.0;
        ONE.y = 0.0;
        ZERO.x = ZERO.y = 0.0;

        // use 'geam' API todo transpose.
        cublasErrcheck(cublasCgeam(cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, col, row,
                                   reinterpret_cast<const float2 *>(&ONE), (float2*)input_matrix, col,
                                   reinterpret_cast<const float2 *>(&ZERO), (float2*)input_matrix, col, (float2*)device_temp, col));
    } else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        matrix_transpose_kernel<thrust::complex<float>> <<<block, thread >>> (row, col, (thrust::complex<float>*)input_matrix, (thrust::complex<float>*)device_temp);

        cudaCheckOnDebug();
    }

    base_device::memory::synchronize_memory_op<std::complex<float>, base_device::DEVICE_GPU, base_device::DEVICE_GPU>()(
        output_matrix,
        device_temp,
        row * col);

    base_device::memory::delete_memory_op<std::complex<float>, base_device::DEVICE_GPU>()(device_temp);

    cudaCheckOnDebug();

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
        double2 ONE, ZERO;
        ONE.x = 1.0;
        ONE.y = 0.0;
        ZERO.x = ZERO.y = 0.0;

        // use 'geam' API todo transpose.
        cublasErrcheck(cublasZgeam(cublas_handle, CUBLAS_OP_T, CUBLAS_OP_N, col, row, &ONE, (double2*)input_matrix, col, &ZERO, (double2*)input_matrix, col, (double2*)device_temp, col));
    } else
    {
        int thread = 1024;
        int block = (row + col + thread - 1) / thread;
        matrix_transpose_kernel<thrust::complex<double>> <<<block, thread >>> (row, col, (thrust::complex<double>*)input_matrix, (thrust::complex<double>*)device_temp);
        cudaCheckOnDebug();
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
    matrix_copy_kernel<double> <<<gridSize, blockSize >>> (n1, n2, A, LDA, B, LDB);
    cudaCheckOnDebug();
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
    matrix_copy_kernel<thrust::complex<float>> <<<gridSize, blockSize >>> (n1, n2, reinterpret_cast<const thrust::complex<float>*>(A), LDA, reinterpret_cast<thrust::complex<float>*>(B), LDB);
    cudaCheckOnDebug();

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
    matrix_copy_kernel<thrust::complex<double>> <<<gridSize, blockSize >>> (n1, n2, reinterpret_cast<const thrust::complex<double>*>(A), LDA, reinterpret_cast<thrust::complex<double>*>(B), LDB);
    cudaCheckOnDebug();
}


// Explicitly instantiate functors for the types of functor registered.

template struct matrixCopy<std::complex<float>, base_device::DEVICE_GPU>;
template struct matrixCopy<double, base_device::DEVICE_GPU>;
template struct matrixCopy<std::complex<double>, base_device::DEVICE_GPU>;
}  // namespace ModuleBase
