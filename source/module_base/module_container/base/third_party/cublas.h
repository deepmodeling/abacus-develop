#ifndef BASE_THIRD_PARTY_CUBLAS_H_
#define BASE_THIRD_PARTY_CUBLAS_H_

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <base/macros/cuda.h>

namespace container {
namespace cuBlasConnector {

static inline
void dot(cublasHandle_t& handle, const int& n, const float *x, const int& incx, const float *y, const int& incy, float* result)
{
    cublasSdot(handle, n, x, incx, y, incy, result);
}
static inline
void dot(cublasHandle_t& handle, const int& n, const double *x, const int& incx, const double *y, const int& incy, double* result)
{
    cublasDdot(handle, n, x, incx, y, incy, result);
}
static inline
void dot(cublasHandle_t& handle, const int& n, const std::complex<float> *x, const int& incx, const std::complex<float> *y, const int& incy, std::complex<float>* result)
{
    cublasCdotc(handle, n, reinterpret_cast<const cuComplex*>(x), incx, reinterpret_cast<const cuComplex*>(y), incy, reinterpret_cast<cuComplex*>(result));
}
static inline
void dot(cublasHandle_t& handle, const int& n, const std::complex<double> *x, const int& incx, const std::complex<double> *y, const int& incy, std::complex<double>* result)
{
    cublasZdotc(handle, n, reinterpret_cast<const cuDoubleComplex*>(x), incx, reinterpret_cast<const cuDoubleComplex*>(y), incy, reinterpret_cast<cuDoubleComplex*>(result));
}

static inline
void axpy(cublasHandle_t& handle, const int& n, const float& alpha, const float *x, const int& incx, float *y, const int& incy)
{
    cublasSaxpy(handle, n, &alpha, x, incx, y, incy);
}
static inline
void axpy(cublasHandle_t& handle, const int& n, const double& alpha, const double *x, const int& incx, double *y, const int& incy)
{
    cublasDaxpy(handle, n, &alpha, x, incx, y, incy);
}
static inline
void axpy(cublasHandle_t& handle, const int& n, const std::complex<float>& alpha, const std::complex<float> *x, const int& incx, std::complex<float> *y, const int& incy)
{
    cublasCaxpy(handle, n, reinterpret_cast<const cuComplex*>(&alpha), reinterpret_cast<const cuComplex*>(x), incx, reinterpret_cast<cuComplex*>(y), incy);
}
static inline
void axpy(cublasHandle_t& handle, const int& n, const std::complex<double>& alpha, const std::complex<double> *x, const int& incx, std::complex<double> *y, const int& incy)
{
    cublasZaxpy(handle, n, reinterpret_cast<const cuDoubleComplex*>(&alpha), reinterpret_cast<const cuDoubleComplex*>(x), incx, reinterpret_cast<cuDoubleComplex*>(y), incy);
}

static inline
void scal(cublasHandle_t& handle, const int& n,  const float& alpha, float *x, const int& incx)
{
    cublasSscal(handle, n, &alpha, x, incx);
}
static inline
void scal(cublasHandle_t& handle, const int& n, const double& alpha, double *x, const int& incx)
{
    cublasDscal(handle, n, &alpha, x, incx);
}
static inline
void scal(cublasHandle_t& handle, const int& n, const std::complex<float>& alpha, std::complex<float> *x, const int& incx)
{
    cublasCscal(handle, n, reinterpret_cast<const cuComplex*>(&alpha), reinterpret_cast<cuComplex*>(x), incx);
}
static inline
void scal(cublasHandle_t& handle, const int& n, const std::complex<double>& alpha, std::complex<double> *x, const int& incx)
{
    cublasZscal(handle, n, reinterpret_cast<const cuDoubleComplex*>(&alpha), reinterpret_cast<cuDoubleComplex*>(x), incx);
}

static inline
void gemv(cublasHandle_t& handle, const char& trans, const int& m, const int& n,
          const float& alpha, const float *A, const int& lda, const float *x, const int& incx,
          const float& beta, float *y, const int& incy)
{
    cublasSgemv(handle, GetCublasOperation(trans), m, n, &alpha, A, lda, x, incx, &beta, y, incy);
}
static inline
void gemv(cublasHandle_t& handle, const char& trans, const int& m, const int& n,
          const double& alpha, const double *A, const int& lda, const double *x, const int& incx,
          const double& beta, double *y, const int& incy)
{
    cublasDgemv(handle, GetCublasOperation(trans), m, n, &alpha, A, lda, x, incx, &beta, y, incy);
}
static inline
void gemv(cublasHandle_t& handle, const char& trans, const int& m, const int& n,
          const std::complex<float>& alpha, const std::complex<float> *A, const int& lda, const std::complex<float> *x, const int& incx,
          const std::complex<float>& beta, std::complex<float> *y, const int& incy)
{
    cublasCgemv(handle, GetCublasOperation(trans), m, n, reinterpret_cast<const cuComplex*>(&alpha), 
    reinterpret_cast<const cuComplex*>(A), lda, reinterpret_cast<const cuComplex*>(x), incx, reinterpret_cast<const cuComplex*>(&beta), reinterpret_cast<cuComplex*>(y), incy);
}
static inline
void gemv(cublasHandle_t& handle, const char& trans, const int& m, const int& n,
          const std::complex<double>& alpha, const std::complex<double> *A, const int& lda, const std::complex<double> *x, const int& incx,
          const std::complex<double>& beta, std::complex<double> *y, const int& incy)
{
    cublasZgemv(handle, GetCublasOperation(trans), m, n, reinterpret_cast<const cuDoubleComplex*>(&alpha), 
    reinterpret_cast<const cuDoubleComplex*>(A), lda, reinterpret_cast<const cuDoubleComplex*>(x), incx, reinterpret_cast<const cuDoubleComplex*>(&beta), reinterpret_cast<cuDoubleComplex*>(y), incy);
}

template <typename T>
static inline
void gemv_batched(cublasHandle_t& handle, const char& trans, const int& m, const int& n,
          const T& alpha, T** A, const int& lda, T** x, const int& incx,
          const T& beta, T** y, const int& incy, const int& batch_size)
{
    for (int ii = 0; ii < batch_size; ++ii) {
        // Call the single GEMV for each pair of matrix A[ii] and vector x[ii]
        cuBlasConnector::gemv(handle, trans, m, n, alpha, A[ii], lda, x[ii], incy, beta, y[ii], incy);
    }
}

template <typename T>
static inline
void gemv_batched_strided(cublasHandle_t& handle, const char& transa, const int& m, const int& n,
          const T& alpha, const T* A, const int& lda, const int& stride_a, const T* x, const int& incx, const int& stride_x,
          const T& beta, T* y, const int& incy, const int& stride_y, const int& batch_size)
{
    for (int ii = 0; ii < batch_size; ii++) {
        // Call the single GEMV for each pair of matrix A[ii] and vector x[ii]
        cuBlasConnector::gemv(handle, transa, m, n, alpha, A + ii * stride_a, lda, x + ii * stride_x, incx, beta, y + ii * stride_y, incy);
    }
}

static inline
void gemm(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const float& alpha, const float* A, const int& lda, const float* B, const int& ldb,
          const float& beta, float* C, const int& ldc)
{
    cublasSgemm(handle, GetCublasOperation(transa), GetCublasOperation(transb),
                m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
static inline
void gemm(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const double& alpha, const double* A, const int& lda, const double* B, const int& ldb,
          const double& beta, double* C, const int& ldc)
{
    cublasDgemm(handle, GetCublasOperation(transa), GetCublasOperation(transb),
                m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc);
}
static inline
void gemm(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<float>& alpha, const std::complex<float>* A, const int& lda, const std::complex<float>* B, const int& ldb,
          const std::complex<float>& beta, std::complex<float>* C, const int& ldc)
{
    cublasCgemm(handle, GetCublasOperation(transa), GetCublasOperation(transb),
                m, n, k, 
                reinterpret_cast<const cuComplex*>(&alpha), 
                reinterpret_cast<const cuComplex*>(A), lda, 
                reinterpret_cast<const cuComplex*>(B), ldb, 
                reinterpret_cast<const cuComplex*>(&beta), 
                reinterpret_cast<cuComplex*>(C), ldc);
}
static inline
void gemm(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<double>& alpha, const std::complex<double>* A, const int& lda, const std::complex<double>* B, const int& ldb,
          const std::complex<double>& beta, std::complex<double>* C, const int& ldc)
{
    cublasZgemm(handle, GetCublasOperation(transa), GetCublasOperation(transb),
                m, n, k,
                reinterpret_cast<const cuDoubleComplex*>(&alpha),  
                reinterpret_cast<const cuDoubleComplex*>(A), lda, 
                reinterpret_cast<const cuDoubleComplex*>(B), ldb, 
                reinterpret_cast<const cuDoubleComplex*>(&beta), 
                reinterpret_cast<cuDoubleComplex*>(C), ldc);
}

template <typename T>
static inline 
T** allocate_(T** in, const int& batch_size)
{
    T** out = nullptr;
    cudaMalloc(reinterpret_cast<void **>(&out), sizeof(T*) * batch_size);
    cudaMemcpy(out, in, sizeof(T*) * batch_size, cudaMemcpyHostToDevice);
    return out;
}

static inline
void gemm_batched(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const float& alpha, float** A, const int& lda, float** B, const int& ldb,
          const float& beta, float** C, const int& ldc, const int& batch_size)
{
    float** d_A = allocate_(A, batch_size);
    float** d_B = allocate_(B, batch_size);
    float** d_C = allocate_(C, batch_size);
    cublasSgemmBatched(handle, GetCublasOperation(transa), GetCublasOperation(transb),
                       m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc, batch_size);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}
static inline
void gemm_batched(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const double& alpha, double** A, const int& lda, double** B, const int& ldb,
          const double& beta, double** C, const int& ldc, const int& batch_size)
{
    double** d_A = allocate_(A, batch_size);
    double** d_B = allocate_(B, batch_size);
    double** d_C = allocate_(C, batch_size);
    cublasDgemmBatched(handle, GetCublasOperation(transa), GetCublasOperation(transb),
                       m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc, batch_size);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}
static inline
void gemm_batched(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<float>& alpha, std::complex<float>** A, const int& lda, std::complex<float>** B, const int& ldb,
          const std::complex<float>& beta, std::complex<float>** C, const int& ldc, const int& batch_size)
{
    std::complex<float>** d_A = allocate_(A, batch_size);
    std::complex<float>** d_B = allocate_(B, batch_size);
    std::complex<float>** d_C = allocate_(C, batch_size);
    cublasCgemmBatched(handle, GetCublasOperation(transa), GetCublasOperation(transb),
                       m, n, k, 
                       reinterpret_cast<const cuComplex*>(&alpha), 
                       reinterpret_cast<cuComplex**>(d_A), lda, 
                       reinterpret_cast<cuComplex**>(d_B), ldb, 
                       reinterpret_cast<const cuComplex*>(&beta), 
                       reinterpret_cast<cuComplex**>(d_C), ldc, batch_size);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}
static inline
void gemm_batched(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<double>& alpha, std::complex<double>** A, const int& lda, std::complex<double>** B, const int& ldb,
          const std::complex<double>& beta, std::complex<double>** C, const int& ldc, const int& batch_size)
{
    std::complex<double>** d_A = allocate_(A, batch_size);
    std::complex<double>** d_B = allocate_(B, batch_size);
    std::complex<double>** d_C = allocate_(C, batch_size);
    cublasZgemmBatched(handle, GetCublasOperation(transa), GetCublasOperation(transb),
                       m, n, k, 
                       reinterpret_cast<const cuDoubleComplex*>(&alpha), 
                       reinterpret_cast<cuDoubleComplex**>(d_A), lda, 
                       reinterpret_cast<cuDoubleComplex**>(d_B), ldb, 
                       reinterpret_cast<const cuDoubleComplex*>(&beta), 
                       reinterpret_cast<cuDoubleComplex**>(d_C), ldc, batch_size);
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);
}

static inline
void gemm_batched_strided(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const float& alpha, const float* A, const int& lda, const int& stride_a, const float* B, const int& ldb, const int& stride_b,
          const float& beta, float* C, const int& ldc, const int& stride_c, const int& batch_size)
{
    cublasSgemmStridedBatched(
            handle, 
            GetCublasOperation(transa), 
            GetCublasOperation(transb),
            m, n, k, 
            &alpha, 
            A, lda, stride_a, 
            B, ldb, stride_b, 
            &beta, 
            C, ldc, stride_c,
            batch_size);
}
static inline
void gemm_batched_strided(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const double& alpha, const double* A, const int& lda, const int& stride_a, const double* B, const int& ldb, const int& stride_b,
          const double& beta, double* C, const int& ldc, const int& stride_c, const int& batch_size)
{
    cublasDgemmStridedBatched(
            handle, 
            GetCublasOperation(transa), 
            GetCublasOperation(transb),
            m, n, k, 
            &alpha, 
            A, lda, stride_a, 
            B, ldb, stride_b, 
            &beta, 
            C, ldc, stride_c,
            batch_size);
}
static inline
void gemm_batched_strided(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<float>& alpha, const std::complex<float>* A, const int& lda, const int& stride_a, const std::complex<float>* B, const int& ldb, const int& stride_b,
          const std::complex<float>& beta, std::complex<float>* C, const int& ldc, const int& stride_c, const int& batch_size)
{
    cublasCgemmStridedBatched(
            handle, 
            GetCublasOperation(transa), 
            GetCublasOperation(transb),
            m, n, k, 
            reinterpret_cast<const cuComplex*>(&alpha), 
            reinterpret_cast<const cuComplex*>(A), lda, stride_a, 
            reinterpret_cast<const cuComplex*>(B), ldb, stride_b, 
            reinterpret_cast<const cuComplex*>(&beta), 
            reinterpret_cast<cuComplex*>(C), ldc, stride_c,
            batch_size);
}
static inline
void gemm_batched_strided(cublasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<double>& alpha, const std::complex<double>* A, const int& lda, const int& stride_a, const std::complex<double>* B, const int& ldb, const int& stride_b,
          const std::complex<double>& beta, std::complex<double>* C, const int& ldc, const int& stride_c, const int& batch_size)
{
    cublasZgemmStridedBatched(
            handle, 
            GetCublasOperation(transa), 
            GetCublasOperation(transb),
            m, n, k, 
            reinterpret_cast<const cuDoubleComplex*>(&alpha), 
            reinterpret_cast<const cuDoubleComplex*>(A), lda, stride_a, 
            reinterpret_cast<const cuDoubleComplex*>(B), ldb, stride_b, 
            reinterpret_cast<const cuDoubleComplex*>(&beta), 
            reinterpret_cast<cuDoubleComplex*>(C), ldc, stride_c,
            batch_size);
}

} // namespace cuBlasConnector
} // namespace container

#endif // BASE_THIRD_PARTY_CUBLAS_H_