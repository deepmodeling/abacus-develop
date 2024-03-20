#ifndef ATEN_KERNELS_BLAS_H_
#define ATEN_KERNELS_BLAS_H_

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_types.h>

#include <base/third_party/blas.h>

namespace container {
namespace kernels {

template <typename T, typename Device>
struct blas_dot {
    void operator()(
        const int& n,
        const T* x,
        const int& incx,
        const T* y,
        const int& incy,
        T* result);
};


template <typename T, typename Device>
struct blas_scal {
    void operator()(
        const int& n,
        const T* alpha,
        T* x,
        const int& incx);
};


template <typename T, typename Device>
struct blas_axpy {
    void operator()(
        const int& n,
        const T* alpha,
        const T* x,
        const int& incx,
        T* y,
        const int& incy);
};


template <typename T, typename Device>
struct blas_gemv {
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
        const int& incy);
};


template <typename T, typename Device>
struct blas_gemv_batched {
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
        const int& batch_size);
};


template <typename T, typename Device>
struct blas_gemv_batched_strided {
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
        const int& batch_size);
};


template <typename T, typename Device>
struct blas_gemm {
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
        const int& ldc);
};


template <typename T, typename Device>
struct blas_gemm_batched {
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
        const int& batch_size);
};


template <typename T, typename Device>
struct blas_gemm_batched_strided {
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
        const int& batch_size);
};

#if __CUDA || __ROCM
void createGpuBlasHandle();  // create blas handle
void destroyGpuBlasHandle(); // destory blas handle
#endif // __CUDA || __UT_USE_CUDA

} // namespace kernels
} // namespace container

#endif // ATEN_KERNELS_BLAS_H_