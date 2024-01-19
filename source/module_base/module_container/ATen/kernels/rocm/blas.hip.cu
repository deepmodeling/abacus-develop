#include <ATen/kernels/blas.h>
#include <base/third_party/blas.h>

#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>

namespace container {
namespace kernels {

static hipblasHandle_t hipblas_handle = nullptr;

void createGpuBlasHandle() {
    if (hipblas_handle == nullptr) {
        hipblasErrcheckInternal(hipblasCreate(&hipblas_handle));
    }
}

void destroyGpuBlasHandle() {
    if (hipblas_handle != nullptr) {
        hipblasErrcheckInternal(hipblasDestroy(hipblas_handle));
        hipblas_handle = nullptr;
    }
}


template <typename T>
struct blas_dot<T, DEVICE_GPU> {
    void operator()(
        const int& n,
        const T* x,
        const int& incx,
        const T* y,
        const int& incy,
        T* result)
    {
        hipBlasConnector::dot(hipblas_handle, n, x, incx, y, incy, result);
    }
};

template <typename T>
struct blas_scal<T, DEVICE_GPU> {
    void operator()(
        const int& n,
        const T* alpha,
        T* x,
        const int& incx)
    {
        hipBlasConnector::scal(hipblas_handle, n, *alpha, x, incx);
    }
};

template <typename T>
struct blas_axpy<T, DEVICE_GPU> {
    void operator()(
        const int& n,
        const T* alpha,
        const T* x,
        const int& incx,
        T* y,
        const int& incy)
    {
        hipBlasConnector::axpy(hipblas_handle, n, *alpha, x, incx, y, incy);
    }
};

template <typename T>
struct blas_gemv<T, DEVICE_GPU> {
    void operator()(
        const char& trans,
        const int& m,
        const int& n,
        const T* alpha,
        const T* A,
        const int& lda,
        const T* x,
        const int& incx,
        const T* beta,
        T* y,
        const int& incy) 
    {
        hipBlasConnector::gemv(hipblas_handle, trans, m, n, *alpha, A, lda, x, incx, *beta, y, incy);
    }
};


template <typename T>
struct blas_gemv_batched<T, DEVICE_GPU> {
    void operator()(
        const char& trans,
        const int& m,
        const int& n,
        const T* alpha,
        T** A,
        const int& lda,
        T** x,
        const int& incx,
        const T* beta,
        T** y,
        const int& incy,
        const int& batch_size)
    {
        hipBlasConnector::gemv_batched(hipblas_handle, trans, m, n, *alpha, A, lda, x, incx, *beta, y, incy, batch_size);
    }
};


template <typename T>
struct blas_gemv_batched_strided<T, DEVICE_GPU> {
    void operator()(
        const char& trans,
        const int& m,
        const int& n,
        const T* alpha,
        const T* A,
        const int& lda,
        const int64_t& stride_a,
        const T* x,
        const int& incx,
        const int64_t& stride_x,
        const T* beta,
        T* y,
        const int& incy,
        const int64_t& stride_y,
        const int& batch_size)
    {
        hipBlasConnector::gemv_batched_strided(hipblas_handle, trans, m, n, *alpha, A, lda, stride_a, x, incx, stride_x, *beta, y, incy, stride_y, batch_size);
    }
};

template <typename T>
struct blas_gemm<T, DEVICE_GPU> {
    void operator()(
        const char& transa,
        const char& transb,
        const int& m,
        const int& n,
        const int& k,
        const T* alpha,
        const T* A,
        const int& lda,
        const T* B,
        const int& ldb,
        const T* beta,
        T* C,
        const int& ldc)
    {
        hipBlasConnector::gemm(hipblas_handle, transa, transb, m, n, k, *alpha, A, lda, B, ldb, *beta, C, ldc);
    }
};

template <typename T>
struct blas_gemm_batched<T, DEVICE_GPU> {
    void operator()(
        const char& transa,
        const char& transb,
        const int& m,
        const int& n,
        const int& k,
        const T* alpha,
        T** A,
        const int& lda,
        T** B,
        const int& ldb,
        const T* beta,
        T** C,
        const int& ldc,
        const int& batch_size)
    {
        hipBlasConnector::gemm_batched(hipblas_handle, transa, transb, m, n, k, *alpha, A, lda, B, ldb, *beta, C, ldc, batch_size);
    }
};

template <typename T>
struct blas_gemm_batched_strided<T, DEVICE_GPU> {
    void operator()(
        const char& transa,
        const char& transb,
        const int& m,
        const int& n,
        const int& k,
        const T* alpha,
        const T* A,
        const int& lda,
        const int& stride_a,
        const T* B,
        const int& ldb,
        const int& stride_b,
        const T* beta,
        T* C,
        const int& ldc,
        const int& stride_c,
        const int& batch_size)
    {
        hipBlasConnector::gemm_batched_strided(hipblas_handle, transa, transb, m, n, k, *alpha, A, lda, stride_a, B, ldb, stride_b, *beta, C, ldc, stride_c, batch_size);
    }
};

// Explicitly instantiate functors for the types of functor registered.
template struct blas_dot<float , DEVICE_GPU>;
template struct blas_dot<double, DEVICE_GPU>;
template struct blas_dot<std::complex<float> , DEVICE_GPU>;
template struct blas_dot<std::complex<double>, DEVICE_GPU>;

template struct blas_scal<float , DEVICE_GPU>;
template struct blas_scal<double, DEVICE_GPU>;
template struct blas_scal<std::complex<float> , DEVICE_GPU>;
template struct blas_scal<std::complex<double>, DEVICE_GPU>;

template struct blas_axpy<float , DEVICE_GPU>;
template struct blas_axpy<double, DEVICE_GPU>;
template struct blas_axpy<std::complex<float> , DEVICE_GPU>;
template struct blas_axpy<std::complex<double>, DEVICE_GPU>;

template struct blas_gemv<float , DEVICE_GPU>;
template struct blas_gemv<double, DEVICE_GPU>;
template struct blas_gemv<std::complex<float >, DEVICE_GPU>;
template struct blas_gemv<std::complex<double>, DEVICE_GPU>;

template struct blas_gemv_batched<float , DEVICE_GPU>;
template struct blas_gemv_batched<double, DEVICE_GPU>;
template struct blas_gemv_batched<std::complex<float >, DEVICE_GPU>;
template struct blas_gemv_batched<std::complex<double>, DEVICE_GPU>;

template struct blas_gemv_batched_strided<float , DEVICE_GPU>;
template struct blas_gemv_batched_strided<double, DEVICE_GPU>;
template struct blas_gemv_batched_strided<std::complex<float >, DEVICE_GPU>;
template struct blas_gemv_batched_strided<std::complex<double>, DEVICE_GPU>;

template struct blas_gemm<float , DEVICE_GPU>;
template struct blas_gemm<double, DEVICE_GPU>;
template struct blas_gemm<std::complex<float >, DEVICE_GPU>;
template struct blas_gemm<std::complex<double>, DEVICE_GPU>;

template struct blas_gemm_batched<float , DEVICE_GPU>;
template struct blas_gemm_batched<double, DEVICE_GPU>;
template struct blas_gemm_batched<std::complex<float >, DEVICE_GPU>;
template struct blas_gemm_batched<std::complex<double>, DEVICE_GPU>;

template struct blas_gemm_batched_strided<float , DEVICE_GPU>;
template struct blas_gemm_batched_strided<double, DEVICE_GPU>;
template struct blas_gemm_batched_strided<std::complex<float >, DEVICE_GPU>;
template struct blas_gemm_batched_strided<std::complex<double>, DEVICE_GPU>;

} // namespace kernels
} // namespace container