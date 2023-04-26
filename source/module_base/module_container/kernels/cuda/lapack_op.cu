#include "../blas_op.h"
#include "../lapack_op.h"

#include <cassert>
#include <complex>

namespace container {
namespace op {

static cusolverDnHandle_t cusolver_H = nullptr;

void createCUSOLVERhandle() {
    if (cusolver_H == nullptr) {
        cusolverErrcheck(cusolverDnCreate(&cusolver_H));
    }
}

void destoryCUSOLVERhandle() {
    if (cusolver_H != nullptr) {
        cusolverErrcheck(cusolverDnDestroy(cusolver_H));
        cusolver_H = nullptr;
    }
}

static inline
void xhegvd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<float> * A, const int& lda,
        std::complex<float> * B, const int& ldb,
        float * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    // int * devInfo = nullptr;
    int lwork = 0;
    // float2 * work = nullptr;
    // cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));
    Tensor devInfo(DataType::DT_INT, DeviceType::GpuDevice, {1});

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnChegvd_bufferSize(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const float2 *>(A), lda,
                                                 reinterpret_cast<const float2 *>(B), ldb, W, &lwork));
    // allocate memery
    // cudaErrcheck(cudaMalloc((void**)&work, sizeof(float2) * lwork));
    Tensor work(DataType::DT_COMPLEX, DeviceType::GpuDevice, {lwork});

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnChegvd(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<float2 *>(A), lda, reinterpret_cast<float2 *>(B), ldb, W,
                                      reinterpret_cast<cuComplex *>(work.data<std::complex<float>>()), lwork, devInfo.data<int>()));

    // cudaErrcheck(cudaMemcpy(&info_gpu, devInfo.data<int>(), sizeof(int), cudaMemcpyDeviceToHost));
    // assert(0 == info_gpu);
    devInfo.to_device<DEVICE_CPU>();
    assert(0 == devInfo.data<int>()[0]);
    // free the buffer
}

static inline
void xhegvd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<double> * A, const int& lda,
        std::complex<double> * B, const int& ldb,
        double * W) {
    // prepare some values for cusolverDnZhegvd_bufferSize
    // int * devInfo = nullptr;
    int lwork = 0;
    // double2 * work = nullptr;
    // cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));
    Tensor devInfo(DataType::DT_INT, DeviceType::GpuDevice, {1});

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZhegvd_bufferSize(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const double2 *>(A), lda,
                                                 reinterpret_cast<const double2 *>(B), ldb, W, &lwork));
    // allocate memery
    // cudaErrcheck(cudaMalloc((void**)&work, sizeof(double2) * lwork));
    Tensor work(DataType::DT_COMPLEX_DOUBLE, DeviceType::GpuDevice, {lwork});

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnZhegvd(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<double2 *>(A), lda, reinterpret_cast<double2 *>(B), ldb, W,
                                      reinterpret_cast<cuDoubleComplex *>(work.data<std::complex<double>>()), lwork,
                                      devInfo.data<int>()));

    // cudaErrcheck(cudaMemcpy(&info_gpu, devInfo.data<int>(), sizeof(int), cudaMemcpyDeviceToHost));
    // assert(0 == info_gpu);
    // free the buffer
    devInfo.to_device<DEVICE_CPU>();
    assert(0 == devInfo.data<int>()[0]);
}

static inline
void xheevd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<float> * A, const int& lda,
        float * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    // int * devInfo = nullptr;
    int lwork = 0;
    // float2 * work = nullptr;
    // cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));
    Tensor devInfo(DataType::DT_INT, DeviceType::GpuDevice, {1});

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnCheevd_bufferSize(cusolver_H, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                     reinterpret_cast<const float2 *>(A), lda, W, &lwork));
    // allocate memery
    // cudaErrcheck(cudaMalloc((void**)&work, sizeof(float2) * lwork));
    // compute eigenvalues and eigenvectors.
    Tensor work(DataType::DT_COMPLEX, DeviceType::GpuDevice, {lwork});
    cusolverErrcheck(cusolverDnCheevd(cusolver_H, CUSOLVER_EIG_MODE_VECTOR, uplo, n, reinterpret_cast<float2 *>(A), lda, W,
                                      reinterpret_cast<cuComplex *>(work.data<std::complex<float>>()), lwork, devInfo.data<int>()));

    // cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    devInfo.to_device<DEVICE_CPU>();
    assert(0 == devInfo.data<int>()[0]);
}

static inline
void xheevd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<double> * A, const int& lda,
        double * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork = 0;
    // double2 * work = nullptr;
    // cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));
    Tensor devInfo(DataType::DT_INT, DeviceType::GpuDevice, {1});

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZheevd_bufferSize(cusolver_H, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                     reinterpret_cast<const double2 *>(A), lda, W, &lwork));
    // allocate memery
    // cudaErrcheck(cudaMalloc((void**)&work, sizeof(double2) * lwork));
    // compute eigenvalues and eigenvectors.
    Tensor work(DataType::DT_COMPLEX, DeviceType::GpuDevice, {lwork});

    cusolverErrcheck(cusolverDnZheevd(cusolver_H, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<double2 *>(A), lda, W,
                                      reinterpret_cast<cuDoubleComplex *>(work.data<std::complex<double>>()), lwork, devInfo.data<int>()));

    // cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    // assert(0 == info_gpu);
    devInfo.to_device<DEVICE_CPU>();
    assert(0 == devInfo.data<int>()[0]);
}

template <typename T>
struct dngvd_op<T, DEVICE_GPU> {
    void operator()(
            const int nstart,
            const int ldh,
            const std::complex<T> *A, // hcc
            const std::complex<T> *B, // scc
            T *W, // eigenvalue
            std::complex<T> *V)
    {
        assert(nstart == ldh);
        // A to V
        cudaErrcheck(cudaMemcpy(V, A, sizeof(std::complex<T>) * ldh * nstart, cudaMemcpyDeviceToDevice));
        xhegvd_wrapper(CUBLAS_FILL_MODE_UPPER, nstart, V, ldh,
                       (std::complex<T> *)B, ldh, W);
    }
};

template <typename T>
struct dnevx_op<T, DEVICE_GPU> {
    void operator()(
            const int nstart,
            const int ldh,
            const std::complex<T> *A, // hcc
            const int m,
            T *W, // eigenvalue
            std::complex<T> *V)
    {
        assert(nstart <= ldh);
        // A to V
        cudaErrcheck(cudaMemcpy(V, A, sizeof(std::complex<T>) * nstart * ldh, cudaMemcpyDeviceToDevice));
        xheevd_wrapper(CUBLAS_FILL_MODE_LOWER, nstart, V, ldh, W);
    }
};

template struct dngvd_op<float, DEVICE_GPU>;
template struct dnevx_op<float, DEVICE_GPU>;
template struct dngvd_op<double, DEVICE_GPU>;
template struct dnevx_op<double, DEVICE_GPU>;

} // namespace op
} // namespace container