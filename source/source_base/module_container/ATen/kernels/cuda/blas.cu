#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>
#include <base/third_party/blas.h>

#include <cuda_runtime.h>
#include <cublas_v2.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace container {
namespace kernels {

static cublasHandle_t cublas_handle = nullptr;

void createGpuBlasHandle() {
    if (cublas_handle == nullptr) {
        CHECK_CUBLAS(cublasCreate(&cublas_handle));
    }
}

void destroyGpuBlasHandle() {
    if (cublas_handle != nullptr) {
        CHECK_CUBLAS(cublasDestroy(cublas_handle));
        cublas_handle = nullptr;
    }
}

namespace {
template <typename T>
__global__ void set_identity(T* matrix, const int n)
{
    const size_t count = static_cast<size_t>(n) * n;
    const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < count; i += stride)
    {
        matrix[i] = i % n == i / n ? T(1) : T(0);
    }
}

template <typename T>
__global__ void copy_triangle(const T* inverse, T* matrix, const int n, const int lda, const bool upper)
{
    const size_t count = static_cast<size_t>(n) * n;
    const size_t stride = static_cast<size_t>(blockDim.x) * gridDim.x;
    for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x; i < count; i += stride)
    {
        const size_t row = i % n;
        const size_t column = i / n;
        if (upper ? row <= column : row >= column)
        {
            matrix[column * lda + row] = inverse[i];
        }
    }
}
} // namespace

template <typename T>
struct lapack_trtri<T, DEVICE_GPU> {
    void operator()(const char& uplo, const char& diag, const int& n, T* matrix, const int& lda)
    {
        if (n < 0 || lda < std::max(1, n))
        {
            throw std::invalid_argument("trtri: invalid matrix dimensions");
        }
        const bool upper = cublas_fill_mode(uplo) == CUBLAS_FILL_MODE_UPPER;
        const bool unit_diag = cublas_diag_type(diag) == CUBLAS_DIAG_UNIT;
        if (n == 0)
        {
            return;
        }
        if (!unit_diag)
        {
            // TRSM does not report singularity; preserve TRTRI's zero-diagonal rejection.
            std::vector<T> diagonal(n);
            const size_t pitch = (static_cast<size_t>(lda) + 1) * sizeof(T);
            CHECK_CUDA(cudaMemcpy2D(diagonal.data(), sizeof(T), matrix, pitch, sizeof(T), n, cudaMemcpyDeviceToHost));
            for (int i = 0; i < n; ++i)
            {
                if (diagonal[i] == T(0))
                {
                    throw std::runtime_error("trtri: zero diagonal at index " + std::to_string(i + 1));
                }
            }
        }

        // Solve A * inverse = I without the Xtrtri workspace overrun seen in CUDA 11.8/12.1.
        Tensor inverse(DataTypeToEnum<T>::value, DeviceType::GpuDevice, {n, n});
        T* data = inverse.data<T>();
        using CudaType = typename GetTypeThrust<T>::type;
        CudaType* cuda_data = reinterpret_cast<CudaType*>(data);
        CudaType* cuda_matrix = reinterpret_cast<CudaType*>(matrix);
        const size_t count = static_cast<size_t>(n) * n;
        const size_t blocks_needed = (count + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        const int blocks = std::min<size_t>(65535, blocks_needed);
        set_identity<<<blocks, THREADS_PER_BLOCK>>>(cuda_data, n);
        CHECK_CUDA(cudaGetLastError());
        cuBlasConnector::trsm(cublas_handle, uplo, diag, n, matrix, lda, data);
        copy_triangle<<<blocks, THREADS_PER_BLOCK>>>(cuda_data, cuda_matrix, n, lda, upper);
        CHECK_CUDA_SYNC();
    }
};

template struct lapack_trtri<float, DEVICE_GPU>;
template struct lapack_trtri<double, DEVICE_GPU>;
template struct lapack_trtri<std::complex<float>, DEVICE_GPU>;
template struct lapack_trtri<std::complex<double>, DEVICE_GPU>;

template <typename T>
struct blas_nrm2<T, DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    Real operator()(
        const int n,
        const T *x,
        const int incx)
    {
        Real result;
        cuBlasConnector::nrm2(cublas_handle, n, x, incx, &result);
        return result;
    }
};

template <typename T>
struct blas_copy<T, DEVICE_GPU> {
    void operator()(
        const int n,
        const T * x,
        const int incx,
        T *y,
        const int incy)
    {
        cuBlasConnector::copy(cublas_handle, n, x, incx, y, incy);
    }
};

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
        cuBlasConnector::dot(cublas_handle, n, x, incx, y, incy, result);
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
        cuBlasConnector::scal(cublas_handle, n, *alpha, x, incx);
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
        cuBlasConnector::axpy(cublas_handle, n, *alpha, x, incx, y, incy);
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
        cuBlasConnector::gemv(cublas_handle, trans, m, n, *alpha, A, lda, x, incx, *beta, y, incy);
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
        cuBlasConnector::gemv_batched(cublas_handle, trans, m, n, *alpha, A, lda, x, incx, *beta, y, incy, batch_size);
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
        cuBlasConnector::gemv_batched_strided(cublas_handle, trans, m, n, *alpha, A, lda, stride_a, x, incx, stride_x, *beta, y, incy, stride_y, batch_size);
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
        cuBlasConnector::gemm(cublas_handle, transa, transb, m, n, k, *alpha, A, lda, B, ldb, *beta, C, ldc);
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
        cuBlasConnector::gemm_batched(cublas_handle, transa, transb, m, n, k, *alpha, A, lda, B, ldb, *beta, C, ldc, batch_size);
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
        cuBlasConnector::gemm_batched_strided(cublas_handle, transa, transb, m, n, k, *alpha, A, lda, stride_a, B, ldb, stride_b, *beta, C, ldc, stride_c, batch_size);
    }
};

// Explicitly instantiate functors for the types of functor registered.



template struct blas_copy<float , DEVICE_GPU>;
template struct blas_copy<double, DEVICE_GPU>;
template struct blas_copy<std::complex<float> , DEVICE_GPU>;
template struct blas_copy<std::complex<double>, DEVICE_GPU>;

template struct blas_nrm2<float , DEVICE_GPU>;
template struct blas_nrm2<double, DEVICE_GPU>;
template struct blas_nrm2<std::complex<float> , DEVICE_GPU>;
template struct blas_nrm2<std::complex<double>, DEVICE_GPU>;

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
