#ifndef BASE_THIRD_PARTY_HIPBLAS_H_
#define BASE_THIRD_PARTY_HIPBLAS_H_

#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <base/macros/rocm.h>

namespace container {
namespace hipBlasConnector {

static inline
void dot(hipblasHandle_t& handle, const int& n, const float *x, const int& incx, const float *y, const int& incy, float* result)
{
    hipblasErrcheck(hipblasSdot(handle, n, x, incx, y, incy, result));
}
static inline
void dot(hipblasHandle_t& handle, const int& n, const double *x, const int& incx, const double *y, const int& incy, double* result)
{
    hipblasErrcheck(hipblasDdot(handle, n, x, incx, y, incy, result));
}
static inline
void dot(hipblasHandle_t& handle, const int& n, const std::complex<float> *x, const int& incx, const std::complex<float> *y, const int& incy, std::complex<float>* result)
{
    hipblasErrcheck(hipblasCdotc(handle, n, reinterpret_cast<const hipblasComplex*>(x), incx, reinterpret_cast<const hipblasComplex*>(y), incy, reinterpret_cast<hipblasComplex*>(result)));
}
static inline
void dot(hipblasHandle_t& handle, const int& n, const std::complex<double> *x, const int& incx, const std::complex<double> *y, const int& incy, std::complex<double>* result)
{
    hipblasErrcheck(hipblasZdotc(handle, n, reinterpret_cast<const hipblasDoubleComplex*>(x), incx, reinterpret_cast<const hipblasDoubleComplex*>(y), incy, reinterpret_cast<hipblasDoubleComplex*>(result)));
}

static inline
void axpy(hipblasHandle_t& handle, const int& n, const float& alpha, const float *x, const int& incx, float *y, const int& incy)
{
    hipblasErrcheck(hipblasSaxpy(handle, n, &alpha, x, incx, y, incy));
}
static inline
void axpy(hipblasHandle_t& handle, const int& n, const double& alpha, const double *x, const int& incx, double *y, const int& incy)
{
    hipblasErrcheck(hipblasDaxpy(handle, n, &alpha, x, incx, y, incy));
}
static inline
void axpy(hipblasHandle_t& handle, const int& n, const std::complex<float>& alpha, const std::complex<float> *x, const int& incx, std::complex<float> *y, const int& incy)
{
    hipblasErrcheck(hipblasCaxpy(handle, n, reinterpret_cast<const hipblasComplex*>(&alpha), reinterpret_cast<const hipblasComplex*>(x), incx, reinterpret_cast<hipblasComplex*>(y), incy));
}
static inline
void axpy(hipblasHandle_t& handle, const int& n, const std::complex<double>& alpha, const std::complex<double> *x, const int& incx, std::complex<double> *y, const int& incy)
{
    hipblasErrcheck(hipblasZaxpy(handle, n, reinterpret_cast<const hipblasDoubleComplex*>(&alpha), reinterpret_cast<const hipblasDoubleComplex*>(x), incx, reinterpret_cast<hipblasDoubleComplex*>(y), incy));
}

static inline
void scal(hipblasHandle_t& handle, const int& n,  const float& alpha, float *x, const int& incx)
{
    hipblasErrcheck(hipblasSscal(handle, n, &alpha, x, incx));
}
static inline
void scal(hipblasHandle_t& handle, const int& n, const double& alpha, double *x, const int& incx)
{
    hipblasErrcheck(hipblasDscal(handle, n, &alpha, x, incx));
}
static inline
void scal(hipblasHandle_t& handle, const int& n, const std::complex<float>& alpha, std::complex<float> *x, const int& incx)
{
    hipblasErrcheck(hipblasCscal(handle, n, reinterpret_cast<const hipblasComplex*>(&alpha), reinterpret_cast<hipblasComplex*>(x), incx));
}
static inline
void scal(hipblasHandle_t& handle, const int& n, const std::complex<double>& alpha, std::complex<double> *x, const int& incx)
{
    hipblasErrcheck(hipblasZscal(handle, n, reinterpret_cast<const hipblasDoubleComplex*>(&alpha), reinterpret_cast<hipblasDoubleComplex*>(x), incx));
}

static inline
void gemv(hipblasHandle_t& handle, const char& trans, const int& m, const int& n,
          const float& alpha, const float *A, const int& lda, const float *x, const int& incx,
          const float& beta, float *y, const int& incy)
{
    hipblasErrcheck(hipblasSgemv(handle, GetHipblasOperation(trans), m, n, &alpha, A, lda, x, incx, &beta, y, incy));
}
static inline
void gemv(hipblasHandle_t& handle, const char& trans, const int& m, const int& n,
          const double& alpha, const double *A, const int& lda, const double *x, const int& incx,
          const double& beta, double *y, const int& incy)
{
    hipblasErrcheck(hipblasDgemv(handle, GetHipblasOperation(trans), m, n, &alpha, A, lda, x, incx, &beta, y, incy));
}
static inline
void gemv(hipblasHandle_t& handle, const char& trans, const int& m, const int& n,
          const std::complex<float>& alpha, const std::complex<float> *A, const int& lda, const std::complex<float> *x, const int& incx,
          const std::complex<float>& beta, std::complex<float> *y, const int& incy)
{
    hipblasErrcheck(hipblasCgemv(handle, GetHipblasOperation(trans), m, n, reinterpret_cast<const hipblasComplex*>(&alpha), 
    reinterpret_cast<const hipblasComplex*>(A), lda, reinterpret_cast<const hipblasComplex*>(x), incx, reinterpret_cast<const hipblasComplex*>(&beta), reinterpret_cast<hipblasComplex*>(y), incy));
}
static inline
void gemv(hipblasHandle_t& handle, const char& trans, const int& m, const int& n,
          const std::complex<double>& alpha, const std::complex<double> *A, const int& lda, const std::complex<double> *x, const int& incx,
          const std::complex<double>& beta, std::complex<double> *y, const int& incy)
{
    hipblasErrcheck(hipblasZgemv(handle, GetHipblasOperation(trans), m, n, reinterpret_cast<const hipblasDoubleComplex*>(&alpha), 
    reinterpret_cast<const hipblasDoubleComplex*>(A), lda, reinterpret_cast<const hipblasDoubleComplex*>(x), incx, reinterpret_cast<const hipblasDoubleComplex*>(&beta), reinterpret_cast<hipblasDoubleComplex*>(y), incy));
}

template <typename T>
static inline
void gemv_batched(hipblasHandle_t& handle, const char& trans, const int& m, const int& n,
          const T& alpha, T** A, const int& lda, T** x, const int& incx,
          const T& beta, T** y, const int& incy, const int& batch_size)
{
    for (int ii = 0; ii < batch_size; ++ii) {
        // Call the single GEMV for each pair of matrix A[ii] and vector x[ii]
        hipBlasConnector::gemv(handle, trans, m, n, alpha, A[ii], lda, x[ii], incy, beta, y[ii], incy);
    }
}

template <typename T>
static inline
void gemv_batched_strided(hipblasHandle_t& handle, const char& transa, const int& m, const int& n,
          const T& alpha, const T* A, const int& lda, const int& stride_a, const T* x, const int& incx, const int& stride_x,
          const T& beta, T* y, const int& incy, const int& stride_y, const int& batch_size)
{
    for (int ii = 0; ii < batch_size; ii++) {
        // Call the single GEMV for each pair of matrix A[ii] and vector x[ii]
        hipBlasConnector::gemv(handle, transa, m, n, alpha, A + ii * stride_a, lda, x + ii * stride_x, incx, beta, y + ii * stride_y, incy);
    }
}

static inline
void gemm(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const float& alpha, const float* A, const int& lda, const float* B, const int& ldb,
          const float& beta, float* C, const int& ldc)
{
    hipblasErrcheck(hipblasSgemm(handle, GetHipblasOperation(transa), GetHipblasOperation(transb),
                m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc));
}
static inline
void gemm(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const double& alpha, const double* A, const int& lda, const double* B, const int& ldb,
          const double& beta, double* C, const int& ldc)
{
    hipblasErrcheck(hipblasDgemm(handle, GetHipblasOperation(transa), GetHipblasOperation(transb),
                m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc));
}
static inline
void gemm(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<float>& alpha, const std::complex<float>* A, const int& lda, const std::complex<float>* B, const int& ldb,
          const std::complex<float>& beta, std::complex<float>* C, const int& ldc)
{
    hipblasErrcheck(hipblasCgemm(handle, GetHipblasOperation(transa), GetHipblasOperation(transb),
                m, n, k, 
                reinterpret_cast<const hipblasComplex*>(&alpha), 
                reinterpret_cast<const hipblasComplex*>(A), lda, 
                reinterpret_cast<const hipblasComplex*>(B), ldb, 
                reinterpret_cast<const hipblasComplex*>(&beta), 
                reinterpret_cast<hipblasComplex*>(C), ldc));
}
static inline
void gemm(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<double>& alpha, const std::complex<double>* A, const int& lda, const std::complex<double>* B, const int& ldb,
          const std::complex<double>& beta, std::complex<double>* C, const int& ldc)
{
    hipblasErrcheck(hipblasZgemm(handle, GetHipblasOperation(transa), GetHipblasOperation(transb),
                m, n, k,
                reinterpret_cast<const hipblasDoubleComplex*>(&alpha),  
                reinterpret_cast<const hipblasDoubleComplex*>(A), lda, 
                reinterpret_cast<const hipblasDoubleComplex*>(B), ldb, 
                reinterpret_cast<const hipblasDoubleComplex*>(&beta), 
                reinterpret_cast<hipblasDoubleComplex*>(C), ldc));
}

template <typename T>
static inline 
T** allocate_(T** in, const int& batch_size)
{
    T** out = nullptr;
    hipErrcheck(hipMalloc(reinterpret_cast<void **>(&out), sizeof(T*) * batch_size));
    hipErrcheck(hipMemcpy(out, in, sizeof(T*) * batch_size, hipMemcpyHostToDevice));
    return out;
}

static inline
void gemm_batched(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const float& alpha, float** A, const int& lda, float** B, const int& ldb,
          const float& beta, float** C, const int& ldc, const int& batch_size)
{
    float** d_A = allocate_(A, batch_size);
    float** d_B = allocate_(B, batch_size);
    float** d_C = allocate_(C, batch_size);
    hipblasErrcheck(hipblasSgemmBatched(handle, GetHipblasOperation(transa), GetHipblasOperation(transb),
                       m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc, batch_size));
    hipErrcheck(hipFree(d_A));
    hipErrcheck(hipFree(d_B));
    hipErrcheck(hipFree(d_C));
}
static inline
void gemm_batched(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const double& alpha, double** A, const int& lda, double** B, const int& ldb,
          const double& beta, double** C, const int& ldc, const int& batch_size)
{
    double** d_A = allocate_(A, batch_size);
    double** d_B = allocate_(B, batch_size);
    double** d_C = allocate_(C, batch_size);
    hipblasErrcheck(hipblasDgemmBatched(handle, GetHipblasOperation(transa), GetHipblasOperation(transb),
                       m, n, k, &alpha, d_A, lda, d_B, ldb, &beta, d_C, ldc, batch_size));
    hipErrcheck(hipFree(d_A));
    hipErrcheck(hipFree(d_B));
    hipErrcheck(hipFree(d_C));
}
static inline
void gemm_batched(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<float>& alpha, std::complex<float>** A, const int& lda, std::complex<float>** B, const int& ldb,
          const std::complex<float>& beta, std::complex<float>** C, const int& ldc, const int& batch_size)
{
    std::complex<float>** d_A = allocate_(A, batch_size);
    std::complex<float>** d_B = allocate_(B, batch_size);
    std::complex<float>** d_C = allocate_(C, batch_size);
    hipblasErrcheck(hipblasCgemmBatched(handle, GetHipblasOperation(transa), GetHipblasOperation(transb),
                       m, n, k, 
                       reinterpret_cast<const hipblasComplex*>(&alpha), 
                       reinterpret_cast<hipblasComplex**>(d_A), lda, 
                       reinterpret_cast<hipblasComplex**>(d_B), ldb, 
                       reinterpret_cast<const hipblasComplex*>(&beta), 
                       reinterpret_cast<hipblasComplex**>(d_C), ldc, batch_size));
    hipErrcheck(hipFree(d_A));
    hipErrcheck(hipFree(d_B));
    hipErrcheck(hipFree(d_C));
}
static inline
void gemm_batched(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<double>& alpha, std::complex<double>** A, const int& lda, std::complex<double>** B, const int& ldb,
          const std::complex<double>& beta, std::complex<double>** C, const int& ldc, const int& batch_size)
{
    std::complex<double>** d_A = allocate_(A, batch_size);
    std::complex<double>** d_B = allocate_(B, batch_size);
    std::complex<double>** d_C = allocate_(C, batch_size);
    hipblasErrcheck(hipblasZgemmBatched(handle, GetHipblasOperation(transa), GetHipblasOperation(transb),
                       m, n, k, 
                       reinterpret_cast<const hipblasDoubleComplex*>(&alpha), 
                       reinterpret_cast<hipblasDoubleComplex**>(d_A), lda, 
                       reinterpret_cast<hipblasDoubleComplex**>(d_B), ldb, 
                       reinterpret_cast<const hipblasDoubleComplex*>(&beta), 
                       reinterpret_cast<hipblasDoubleComplex**>(d_C), ldc, batch_size));
    hipErrcheck(hipFree(d_A));
    hipErrcheck(hipFree(d_B));
    hipErrcheck(hipFree(d_C));
}

static inline
void gemm_batched_strided(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const float& alpha, const float* A, const int& lda, const int& stride_a, const float* B, const int& ldb, const int& stride_b,
          const float& beta, float* C, const int& ldc, const int& stride_c, const int& batch_size)
{
    hipblasErrcheck(hipblasSgemmStridedBatched(
            handle, 
            GetHipblasOperation(transa), 
            GetHipblasOperation(transb),
            m, n, k, 
            &alpha, 
            A, lda, stride_a, 
            B, ldb, stride_b, 
            &beta, 
            C, ldc, stride_c,
            batch_size));
}
static inline
void gemm_batched_strided(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const double& alpha, const double* A, const int& lda, const int& stride_a, const double* B, const int& ldb, const int& stride_b,
          const double& beta, double* C, const int& ldc, const int& stride_c, const int& batch_size)
{
    hipblasErrcheck(hipblasDgemmStridedBatched(
            handle, 
            GetHipblasOperation(transa), 
            GetHipblasOperation(transb),
            m, n, k, 
            &alpha, 
            A, lda, stride_a, 
            B, ldb, stride_b, 
            &beta, 
            C, ldc, stride_c,
            batch_size));
}
static inline
void gemm_batched_strided(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<float>& alpha, const std::complex<float>* A, const int& lda, const int& stride_a, const std::complex<float>* B, const int& ldb, const int& stride_b,
          const std::complex<float>& beta, std::complex<float>* C, const int& ldc, const int& stride_c, const int& batch_size)
{
    hipblasErrcheck(hipblasCgemmStridedBatched(
            handle, 
            GetHipblasOperation(transa), 
            GetHipblasOperation(transb),
            m, n, k, 
            reinterpret_cast<const hipblasComplex*>(&alpha), 
            reinterpret_cast<const hipblasComplex*>(A), lda, stride_a, 
            reinterpret_cast<const hipblasComplex*>(B), ldb, stride_b, 
            reinterpret_cast<const hipblasComplex*>(&beta), 
            reinterpret_cast<hipblasComplex*>(C), ldc, stride_c,
            batch_size));
}
static inline
void gemm_batched_strided(hipblasHandle_t& handle, const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<double>& alpha, const std::complex<double>* A, const int& lda, const int& stride_a, const std::complex<double>* B, const int& ldb, const int& stride_b,
          const std::complex<double>& beta, std::complex<double>* C, const int& ldc, const int& stride_c, const int& batch_size)
{
    hipblasErrcheck(hipblasZgemmStridedBatched(
            handle, 
            GetHipblasOperation(transa), 
            GetHipblasOperation(transb),
            m, n, k, 
            reinterpret_cast<const hipblasDoubleComplex*>(&alpha), 
            reinterpret_cast<const hipblasDoubleComplex*>(A), lda, stride_a, 
            reinterpret_cast<const hipblasDoubleComplex*>(B), ldb, stride_b, 
            reinterpret_cast<const hipblasDoubleComplex*>(&beta), 
            reinterpret_cast<hipblasDoubleComplex*>(C), ldc, stride_c,
            batch_size));
}

} // namespace hipBlasConnector
} // namespace container

#endif // BASE_THIRD_PARTY_HIPBLAS_H_