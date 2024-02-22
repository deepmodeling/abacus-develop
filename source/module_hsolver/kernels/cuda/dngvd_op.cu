#include "module_hsolver/kernels/dngvd_op.h"
#include "helper_cuda.h"

#include <base/macros/macros.h>

#include <cusolverDn.h>

#include "module_base/memory_cuda.h"

namespace hsolver
{

static cusolverDnHandle_t cusolver_H = nullptr;

void createGpuSolverHandle()
{
    if (cusolver_H == nullptr)
    {
        cusolverErrcheck(cusolverDnCreate(&cusolver_H));
    }
}

void destroyGpuSolverHandle()
{
    if (cusolver_H != nullptr)
    {
        cusolverErrcheck(cusolverDnDestroy(cusolver_H));
        cusolver_H = nullptr;
    }
}

static inline
void xhegvd_wrapper(
    const cublasFillMode_t& uplo,
    const int& n,
    double* A, const int& lda,
    double* B, const int& ldb,
    double* W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int* devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    double* work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::xhegvdW",sizeof(int));
    

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnDsygvd_bufferSize(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
        A, lda, B, ldb, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(double) * lwork));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::xhegvdW",sizeof(double) * lwork);

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnDsygvd(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
        A, lda, B, ldb, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    // free the buffer
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
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
    int * devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    float2 * work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::xhegvdW",sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnChegvd_bufferSize(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const float2 *>(A), lda,
                                                 reinterpret_cast<const float2 *>(B), ldb, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(float2) * lwork));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::xhegvdW",sizeof(float2) * lwork);

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnChegvd(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<float2 *>(A), lda, reinterpret_cast<float2 *>(B), ldb, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    // free the buffer
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
}

static inline
void xhegvd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<double> * A, const int& lda,
        std::complex<double> * B, const int& ldb,
        double * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int * devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    double2 * work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::xhegvdW",sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZhegvd_bufferSize(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const double2 *>(A), lda,
                                                 reinterpret_cast<const double2 *>(B), ldb, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(double2) * lwork));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::xhegvdW",sizeof(double2) * lwork);

    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnZhegvd(cusolver_H, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<double2 *>(A), lda, reinterpret_cast<double2 *>(B), ldb, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    // free the buffer
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
}

static inline
void xheevd_wrapper(
    const cublasFillMode_t& uplo,
    const int& n,
    double* A, const int& lda,
    double* W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int* devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    double* work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::Dngvd::xheevdW",sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnDsyevd_bufferSize(cusolver_H, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
        A, lda, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(double) * lwork));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::Dngvd::xheevdW",sizeof(double) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnDsyevd(cusolver_H, CUSOLVER_EIG_MODE_VECTOR, uplo, n, A, lda, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
}

static inline
void xheevd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<float> * A, const int& lda,
        float * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int * devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    float2 * work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::xheevdW",sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnCheevd_bufferSize(cusolver_H, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const float2 *>(A), lda, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(float2) * lwork));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::xheevdW",sizeof(float2) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnCheevd(cusolver_H, CUSOLVER_EIG_MODE_VECTOR, uplo, n, reinterpret_cast<float2 *>(A), lda, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
}

static inline
void xheevd_wrapper (
        const cublasFillMode_t& uplo,
        const int& n,
        std::complex<double> * A, const int& lda,
        double * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int * devInfo = nullptr;
    int lwork = 0, info_gpu = 0;
    double2 * work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&devInfo, sizeof(int)));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::xheevdW",sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZheevd_bufferSize(cusolver_H, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                                 reinterpret_cast<const double2 *>(A), lda, W, &lwork));
    // allocate memery
    cudaErrcheck(cudaMalloc((void**)&work, sizeof(double2) * lwork));
    ModuleBase::Memory_CUDA::record("dngvd_op","Dngvd::xheevdW",sizeof(double2) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnZheevd(cusolver_H, CUSOLVER_EIG_MODE_VECTOR, uplo, n,
                                      reinterpret_cast<double2 *>(A), lda, W, work, lwork, devInfo));

    cudaErrcheck(cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost));
    assert(0 == info_gpu);
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(devInfo));
}

template <typename T>
struct dngvd_op<T, psi::DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const psi::DEVICE_GPU* d,
        const int nstart,
        const int ldh,
        const T* A, // hcc
        const T* B, // scc
        Real* W, // eigenvalue
        T* V)
    {
        assert(nstart == ldh);
        // A to V
        cudaErrcheck(cudaMemcpy(V, A, sizeof(T) * ldh * nstart, cudaMemcpyDeviceToDevice));
        //ModuleBase::Memory_CUDA::record("DngvdOp","DngvdOp",sizeof(T) * ldh * nstart);
        xhegvd_wrapper(CUBLAS_FILL_MODE_UPPER, nstart, V, ldh,
            (T*)B, ldh, W);
    }
};

template <typename T>
struct dnevx_op<T, psi::DEVICE_GPU> {
    using Real = typename GetTypeReal<T>::type;
    void operator()(
        const psi::DEVICE_GPU* d,
        const int nstart,
        const int ldh,
        const T* A, // hcc
        const int m,
        Real* W, // eigenvalue
        T* V)
    {
        assert(nstart <= ldh);
        // A to V
        cudaErrcheck(cudaMemcpy(V, A, sizeof(T) * nstart * ldh, cudaMemcpyDeviceToDevice));
        //ModuleBase::Memory_CUDA::record("DngvdOp","DnevxOp",sizeof(T) * nstart * ldh);
        xheevd_wrapper(CUBLAS_FILL_MODE_LOWER, nstart, V, ldh, W);
    }
};

template struct dngvd_op<std::complex<float>, psi::DEVICE_GPU>;
template struct dnevx_op<std::complex<float>, psi::DEVICE_GPU>;
template struct dngvd_op<std::complex<double>, psi::DEVICE_GPU>;
template struct dnevx_op<std::complex<double>, psi::DEVICE_GPU>;

#ifdef __LCAO
template struct dngvd_op<double, psi::DEVICE_GPU>;
template struct dnevx_op<double, psi::DEVICE_GPU>;
#endif

} // namespace hsolver