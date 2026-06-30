#ifndef BASE_THIRD_PARTY_CUSOLVER_H_
#define BASE_THIRD_PARTY_CUSOLVER_H_

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <iostream>
#include <memory>

// #include <base/third_party/cusolver_utils.h> // traits, needed if generic API is used.
// header provided by cusolver, including some data types and macros.
// see https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuSOLVER/utils/cusolver_utils.h
// The cuSolverDN library provides two different APIs; legacy and generic.
// https://docs.nvidia.com/cuda/cusolver/index.html#naming-conventions
// now only legacy APIs are used, while the general APIs have the potential to simplify code implementation.
// for example, cucusolverDnXpotrf/getrf/geqrf/sytrf
// More tests are needed to confirm that the generic APIs are operating normally, as they are not yet fully supported.

#include <base/macros/cuda.h>

namespace container {
namespace cuSolverConnector {

struct CudaDeleter {
    void operator()(void* ptr) const noexcept {
        if (ptr) cudaFree(ptr);
    }
};

template <typename T>
using unique_cuda_ptr = std::unique_ptr<T, CudaDeleter>;

struct HostDeleter {
    void operator()(void* ptr) const noexcept {
        if (ptr) free(ptr);
    }
};

template <typename T>
using unique_host_ptr = std::unique_ptr<T, HostDeleter>;

#if CUDA_VERSION >= 11000
// Generic API (CUDA 11.0+)
template <typename T>
static inline
void trtri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, T* A, const int& lda)
{
    size_t d_lwork = 0, h_lwork = 0;
    using Type = typename GetTypeThrust<T>::type;
    CHECK_CUSOLVER(cusolverDnXtrtri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), cublas_diag_type(diag), n, GetTypeCuda<T>::cuda_data_type, reinterpret_cast<Type*>(A), lda, &d_lwork, &h_lwork));
    void* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, d_lwork));
    unique_cuda_ptr<void> d_work(raw_d_work);
    unique_host_ptr<void> h_work;
    if (h_lwork) {
        void* raw_h_work = malloc(h_lwork);
        if (raw_h_work == nullptr) {
            throw std::bad_alloc();
        }
        h_work.reset(raw_h_work);
    }
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);
    CHECK_CUSOLVER(cusolverDnXtrtri(cusolver_handle, cublas_fill_mode(uplo), cublas_diag_type(diag), n, GetTypeCuda<T>::cuda_data_type, reinterpret_cast<Type*>(A), n, d_work.get(), d_lwork, h_work.get(), h_lwork, d_info.get()));
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("trtri: failed to invert matrix");
    }
}
#else
#error "CUDA version < 11.0 is not supported. cusolverDnXtrtri (CUDA 11.0+) is required."
#endif

static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, float * A, const int& lda)
{
    int lwork;
    CHECK_CUSOLVER(cusolverDnSpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    float* raw_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_work, lwork * sizeof(float)));
    unique_cuda_ptr<float> work(raw_work);
    CHECK_CUSOLVER(cusolverDnSpotri(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work.get(), lwork, nullptr));
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, double * A, const int& lda)
{
    int lwork;
    CHECK_CUSOLVER(cusolverDnDpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    double* raw_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_work, lwork * sizeof(double)));
    unique_cuda_ptr<double> work(raw_work);
    CHECK_CUSOLVER(cusolverDnDpotri(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work.get(), lwork, nullptr));
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    CHECK_CUSOLVER(cusolverDnCpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex *>(A), n, &lwork));
    cuComplex* raw_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_work, lwork * sizeof(cuComplex)));
    unique_cuda_ptr<cuComplex> work(raw_work);
    CHECK_CUSOLVER(cusolverDnCpotri(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex *>(A), n, work.get(), lwork, nullptr));
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    CHECK_CUSOLVER(cusolverDnZpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex *>(A), n, &lwork));
    cuDoubleComplex* raw_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_work, lwork * sizeof(cuDoubleComplex)));
    unique_cuda_ptr<cuDoubleComplex> work(raw_work);
    CHECK_CUSOLVER(cusolverDnZpotri(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex *>(A), n, work.get(), lwork, nullptr));
}


static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, float * A, const int& lda)
{
    int lwork;
    int* raw_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_info, 1 * sizeof(int)));
    unique_cuda_ptr<int> info(raw_info);
    CHECK_CUSOLVER(cusolverDnSpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    float* raw_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_work, lwork * sizeof(float)));
    unique_cuda_ptr<float> work(raw_work);
    CHECK_CUSOLVER(cusolverDnSpotrf(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work.get(), lwork, info.get()));
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, double * A, const int& lda)
{
    int lwork;
    int* raw_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_info, 1 * sizeof(int)));
    unique_cuda_ptr<int> info(raw_info);
    CHECK_CUSOLVER(cusolverDnDpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    double* raw_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_work, lwork * sizeof(double)));
    unique_cuda_ptr<double> work(raw_work);
    CHECK_CUSOLVER(cusolverDnDpotrf(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work.get(), lwork, info.get()));
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    int* raw_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_info, 1 * sizeof(int)));
    unique_cuda_ptr<int> info(raw_info);
    CHECK_CUSOLVER(cusolverDnCpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex*>(A), lda, &lwork));
    cuComplex* raw_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_work, lwork * sizeof(cuComplex)));
    unique_cuda_ptr<cuComplex> work(raw_work);
    CHECK_CUSOLVER(cusolverDnCpotrf(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex*>(A), lda, work.get(), lwork, info.get()));
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    int* raw_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_info, 1 * sizeof(int)));
    unique_cuda_ptr<int> info(raw_info);
    CHECK_CUSOLVER(cusolverDnZpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex*>(A), lda, &lwork));
    cuDoubleComplex* raw_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_work, lwork * sizeof(cuDoubleComplex)));
    unique_cuda_ptr<cuDoubleComplex> work(raw_work);
    CHECK_CUSOLVER(cusolverDnZpotrf(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex*>(A), lda, work.get(), lwork, info.get()));
}


static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float * W)
{
    int lwork  = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnSsyevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, &lwork));
    float* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(float) * lwork));
    unique_cuda_ptr<float> d_work(raw_d_work);
    CHECK_CUSOLVER(cusolverDnSsyevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, d_work.get(), lwork, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
}
static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double * W)
{
    int lwork  = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnDsyevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, &lwork));
    double* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(double) * lwork));
    unique_cuda_ptr<double> d_work(raw_d_work);
    CHECK_CUSOLVER(cusolverDnDsyevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, d_work.get(), lwork, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
}
static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, float * W)
{
    int lwork  = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnCheevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, W, &lwork));
    cuComplex* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(cuComplex) * lwork));
    unique_cuda_ptr<cuComplex> d_work(raw_d_work);
    CHECK_CUSOLVER(cusolverDnCheevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, W, d_work.get(), lwork, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
}
static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, double* W)
{
    int lwork  = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnZheevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, W, &lwork));
    cuDoubleComplex* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(cuDoubleComplex) * lwork));
    unique_cuda_ptr<cuDoubleComplex> d_work(raw_d_work);
    CHECK_CUSOLVER(cusolverDnZheevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, W, d_work.get(), lwork, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
}

// =====================================================================================================
// heevdx: Compute eigenvalues and eigenvectors of symmetric/Hermitian matrix
// =====================================================================================================
// --- float ---
static inline
void heevdx(cusolverDnHandle_t& cusolver_handle,
    const int n,
    const int lda,
    float* d_A,
    const char jobz,
    const char uplo,
    const char range,
    const int il, const int iu,
    const float vl, const float vu,
    float* d_eigen_val,
    int* h_meig)
{
    int lwork = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);
    cusolverEigRange_t range_t = cublas_eig_range(range);

    CHECK_CUSOLVER(cusolverDnSsyevdx_bufferSize(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n, d_A, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    float* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(float) * lwork));
    unique_cuda_ptr<float> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnSsyevdx(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        d_A, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work.get(), lwork,
        d_info.get()
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevdx (float) failed with info = " + std::to_string(h_info));
    }
}

// --- double ---
static inline
void heevdx(cusolverDnHandle_t& cusolver_handle,
    const int n,
    const int lda,
    double* d_A,
    const char jobz,
    const char uplo,
    const char range,
    const int il, const int iu,
    const double vl, const double vu,
    double* d_eigen_val,
    int* h_meig)
{
    int lwork = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);
    cusolverEigRange_t range_t = cublas_eig_range(range);

    CHECK_CUSOLVER(cusolverDnDsyevdx_bufferSize(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n, d_A, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    double* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(double) * lwork));
    unique_cuda_ptr<double> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnDsyevdx(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        d_A, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work.get(), lwork,
        d_info.get()
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevdx (double) failed with info = " + std::to_string(h_info));
    }
}

// --- complex<float> ---
static inline
void heevdx(cusolverDnHandle_t& cusolver_handle,
    const int n,
    const int lda,
    std::complex<float>* d_A,
    const char jobz,
    const char uplo,
    const char range,
    const int il, const int iu,
    const float vl, const float vu,
    float* d_eigen_val,
    int* h_meig)
{
    int lwork = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);
    cusolverEigRange_t range_t = cublas_eig_range(range);

    CHECK_CUSOLVER(cusolverDnCheevdx_bufferSize(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        reinterpret_cast<cuComplex*>(d_A), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    cuComplex* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(cuComplex) * lwork));
    unique_cuda_ptr<cuComplex> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnCheevdx(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        reinterpret_cast<cuComplex*>(d_A), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work.get(), lwork,
        d_info.get()
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevdx (complex<float>) failed with info = " + std::to_string(h_info));
    }
}

// --- complex<double> ---
static inline
void heevdx(cusolverDnHandle_t& cusolver_handle,
    const int n,
    const int lda,
    std::complex<double>* d_A,
    const char jobz,
    const char uplo,
    const char range,
    const int il, const int iu,
    const double vl, const double vu,
    double* d_eigen_val,
    int* h_meig)
{
    int lwork = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);
    cusolverEigRange_t range_t = cublas_eig_range(range);

    CHECK_CUSOLVER(cusolverDnZheevdx_bufferSize(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        reinterpret_cast<cuDoubleComplex*>(d_A), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    cuDoubleComplex* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(cuDoubleComplex) * lwork));
    unique_cuda_ptr<cuDoubleComplex> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnZheevdx(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        reinterpret_cast<cuDoubleComplex*>(d_A), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work.get(), lwork,
        d_info.get()
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevdx (complex<double>) failed with info = " + std::to_string(h_info));
    }
}

static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float* B, const int& ldb, float * W)
{
    int lwork  = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnSsygvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, &lwork));
    float* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(float) * lwork));
    unique_cuda_ptr<float> d_work(raw_d_work);
    CHECK_CUSOLVER(cusolverDnSsygvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, d_work.get(), lwork, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
}
static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double* B, const int& ldb, double * W)
{
    int lwork  = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnDsygvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, &lwork));
    double* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(double) * lwork));
    unique_cuda_ptr<double> d_work(raw_d_work);
    CHECK_CUSOLVER(cusolverDnDsygvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, d_work.get(), lwork, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
}
static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, float* W)
{
    int lwork  = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnChegvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, reinterpret_cast<cuComplex*>(B), ldb, W, &lwork));
    cuComplex* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(cuComplex) * lwork));
    unique_cuda_ptr<cuComplex> d_work(raw_d_work);
    CHECK_CUSOLVER(cusolverDnChegvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, reinterpret_cast<cuComplex*>(B), ldb, W, d_work.get(), lwork, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
}
static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, double* W)
{
    int lwork  = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnZhegvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, reinterpret_cast<cuDoubleComplex*>(B), ldb, W, &lwork));
    cuDoubleComplex* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(cuDoubleComplex) * lwork));
    unique_cuda_ptr<cuDoubleComplex> d_work(raw_d_work);
    CHECK_CUSOLVER(cusolverDnZhegvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, reinterpret_cast<cuDoubleComplex*>(B), ldb, W, d_work.get(), lwork, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
}

// =====================================================================================================
// hegvd x: Compute selected eigenvalues and eigenvectors of generalized Hermitian-definite eigenproblem
//          A * x = lambda * B * x
// =====================================================================================================

// --- float ---
static inline
void hegvdx(
    cusolverDnHandle_t& cusolver_handle,
    const int itype,
    const char jobz,
    const char range,
    const char uplo,
    const int n,
    const int lda,
    float* d_A,
    float* d_B,
    const float vl,
    const float vu,
    const int il,
    const int iu,
    int* h_meig,
    float* d_eigen_val,
    float* d_eigen_vec
) {
    int lwork = 0;
    int *raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    float *raw_d_A_copy = nullptr, *raw_d_B_copy = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_A_copy, sizeof(float) * n * lda));
    unique_cuda_ptr<float> d_A_copy(raw_d_A_copy);
    CHECK_CUDA(cudaMalloc((void**)&raw_d_B_copy, sizeof(float) * n * lda));
    unique_cuda_ptr<float> d_B_copy(raw_d_B_copy);
    CHECK_CUDA(cudaMemcpy(d_A_copy.get(), d_A, sizeof(float) * n * lda, cudaMemcpyDeviceToDevice));
    CHECK_CUDA(cudaMemcpy(d_B_copy.get(), d_B, sizeof(float) * n * lda, cudaMemcpyDeviceToDevice));

    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    CHECK_CUSOLVER(cusolverDnSsygvdx_bufferSize(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy.get(), lda,
        d_B_copy.get(), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    float* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(float) * lwork));
    unique_cuda_ptr<float> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnSsygvdx(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy.get(), lda,
        d_B_copy.get(), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work.get(), lwork,
        d_info.get()
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info < 0) {
        throw std::runtime_error("hegvdx (float): illegal argument #" + std::to_string(-h_info));
    } else if (h_info > 0) {
        if (jobz_t == CUSOLVER_EIG_MODE_NOVECTOR && h_info <= n) {
            throw std::runtime_error("hegvdx (float): failed to converge, " + std::to_string(h_info) + " off-diagonal elements didn't converge");
        } else if (h_info > n) {
            throw std::runtime_error("hegvdx (float): leading minor of order " + std::to_string(h_info - n) + " of B is not positive definite");
        }
    }

    if (jobz == 'V') {
        const int m = (*h_meig);
        CHECK_CUDA(cudaMemcpy(d_eigen_vec, d_A_copy.get(), sizeof(float) * n * m, cudaMemcpyDeviceToDevice));
    }
}


// --- double ---
static inline
void hegvdx(
    cusolverDnHandle_t& cusolver_handle,
    const int itype,
    const char jobz,
    const char range,
    const char uplo,
    const int n,
    const int lda,
    double* d_A,
    double* d_B,
    const double vl,
    const double vu,
    const int il,
    const int iu,
    int* h_meig,
    double* d_eigen_val,
    double* d_eigen_vec
) {
    int lwork = 0;
    int *raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    double *raw_d_A_copy = nullptr, *raw_d_B_copy = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_A_copy, sizeof(double) * n * lda));
    unique_cuda_ptr<double> d_A_copy(raw_d_A_copy);
    CHECK_CUDA(cudaMalloc((void**)&raw_d_B_copy, sizeof(double) * n * lda));
    unique_cuda_ptr<double> d_B_copy(raw_d_B_copy);
    CHECK_CUDA(cudaMemcpy(d_A_copy.get(), d_A, sizeof(double) * n * lda, cudaMemcpyDeviceToDevice));
    CHECK_CUDA(cudaMemcpy(d_B_copy.get(), d_B, sizeof(double) * n * lda, cudaMemcpyDeviceToDevice));

    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    CHECK_CUSOLVER(cusolverDnDsygvdx_bufferSize(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy.get(), lda,
        d_B_copy.get(), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    double* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(double) * lwork));
    unique_cuda_ptr<double> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnDsygvdx(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy.get(), lda,
        d_B_copy.get(), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work.get(), lwork,
        d_info.get()
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info < 0) {
        throw std::runtime_error("hegvdx (double): illegal argument #" + std::to_string(-h_info));
    } else if (h_info > 0) {
        if (jobz_t == CUSOLVER_EIG_MODE_NOVECTOR && h_info <= n) {
            throw std::runtime_error("hegvdx (double): failed to converge, " + std::to_string(h_info) + " off-diagonal elements didn't converge");
        } else if (h_info > n) {
            throw std::runtime_error("hegvdx (double): leading minor of order " + std::to_string(h_info - n) + " of B is not positive definite");
        }
    }

    if (jobz == 'V') {
        const int m = (*h_meig);
        CHECK_CUDA(cudaMemcpy(d_eigen_vec, d_A_copy.get(), sizeof(double) * n * m, cudaMemcpyDeviceToDevice));
    }
}


// --- complex<float> ---
static inline
void hegvdx(
    cusolverDnHandle_t& cusolver_handle,
    const int itype,
    const char jobz,
    const char range,
    const char uplo,
    const int n,
    const int lda,
    std::complex<float>* d_A,
    std::complex<float>* d_B,
    const float vl,
    const float vu,
    const int il,
    const int iu,
    int* h_meig,
    float* d_eigen_val,
    std::complex<float>* d_eigen_vec
) {
    int lwork = 0;
    int *raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    cuComplex *raw_d_A_copy = nullptr, *raw_d_B_copy = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_A_copy, sizeof(cuComplex) * n * lda));
    unique_cuda_ptr<cuComplex> d_A_copy(raw_d_A_copy);
    CHECK_CUDA(cudaMalloc((void**)&raw_d_B_copy, sizeof(cuComplex) * n * lda));
    unique_cuda_ptr<cuComplex> d_B_copy(raw_d_B_copy);
    CHECK_CUDA(cudaMemcpy(d_A_copy.get(), reinterpret_cast<cuComplex*>(d_A), sizeof(cuComplex) * n * lda, cudaMemcpyDeviceToDevice));
    CHECK_CUDA(cudaMemcpy(d_B_copy.get(), reinterpret_cast<cuComplex*>(d_B), sizeof(cuComplex) * n * lda, cudaMemcpyDeviceToDevice));

    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    CHECK_CUSOLVER(cusolverDnChegvdx_bufferSize(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy.get(), lda,
        d_B_copy.get(), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    cuComplex* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(cuComplex) * lwork));
    unique_cuda_ptr<cuComplex> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnChegvdx(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy.get(), lda,
        d_B_copy.get(), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work.get(), lwork,
        d_info.get()
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info < 0) {
        throw std::runtime_error("hegvdx (complex<float>): illegal argument #" + std::to_string(-h_info));
    } else if (h_info > 0) {
        if (jobz_t == CUSOLVER_EIG_MODE_NOVECTOR && h_info <= n) {
            throw std::runtime_error("hegvdx (complex<float>): failed to converge, " + std::to_string(h_info) + " off-diagonal elements didn't converge");
        } else if (h_info > n) {
            throw std::runtime_error("hegvdx (complex<float>): leading minor of order " + std::to_string(h_info - n) + " of B is not positive definite");
        }
    }

    if (jobz == 'V') {
        const int m = (*h_meig);
        CHECK_CUDA(cudaMemcpy(reinterpret_cast<cuComplex*>(d_eigen_vec), d_A_copy.get(), sizeof(cuComplex) * n * m, cudaMemcpyDeviceToDevice));
    }
}


// --- complex<double> ---
static inline
void hegvdx(
    cusolverDnHandle_t& cusolver_handle,
    const int itype,
    const char jobz,
    const char range,
    const char uplo,
    const int n,
    const int lda,
    std::complex<double>* d_A,
    std::complex<double>* d_B,
    const double vl,
    const double vu,
    const int il,
    const int iu,
    int* h_meig,
    double* d_eigen_val,
    std::complex<double>* d_eigen_vec
) {
    int lwork = 0;
    int *raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    cuDoubleComplex *raw_d_A_copy = nullptr, *raw_d_B_copy = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_A_copy, sizeof(cuDoubleComplex) * n * lda));
    unique_cuda_ptr<cuDoubleComplex> d_A_copy(raw_d_A_copy);
    CHECK_CUDA(cudaMalloc((void**)&raw_d_B_copy, sizeof(cuDoubleComplex) * n * lda));
    unique_cuda_ptr<cuDoubleComplex> d_B_copy(raw_d_B_copy);
    CHECK_CUDA(cudaMemcpy(d_A_copy.get(), reinterpret_cast<cuDoubleComplex*>(d_A), sizeof(cuDoubleComplex) * n * lda, cudaMemcpyDeviceToDevice));
    CHECK_CUDA(cudaMemcpy(d_B_copy.get(), reinterpret_cast<cuDoubleComplex*>(d_B), sizeof(cuDoubleComplex) * n * lda, cudaMemcpyDeviceToDevice));

    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    CHECK_CUSOLVER(cusolverDnZhegvdx_bufferSize(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy.get(), lda,
        d_B_copy.get(), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    cuDoubleComplex* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(cuDoubleComplex) * lwork));
    unique_cuda_ptr<cuDoubleComplex> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnZhegvdx(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy.get(), lda,
        d_B_copy.get(), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work.get(), lwork,
        d_info.get()
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info < 0) {
        throw std::runtime_error("hegvdx (complex<double>): illegal argument #" + std::to_string(-h_info));
    } else if (h_info > 0) {
        if (jobz_t == CUSOLVER_EIG_MODE_NOVECTOR && h_info <= n) {
            throw std::runtime_error("hegvdx (complex<double>): failed to converge, " + std::to_string(h_info) + " off-diagonal elements didn't converge");
        } else if (h_info > n) {
            throw std::runtime_error("hegvdx (complex<double>): leading minor of order " + std::to_string(h_info - n) + " of B is not positive definite");
        }
    }

    if (jobz == 'V') {
        const int m = (*h_meig);
        CHECK_CUDA(cudaMemcpy(reinterpret_cast<cuDoubleComplex*>(d_eigen_vec), d_A_copy.get(), sizeof(cuDoubleComplex) * n * m, cudaMemcpyDeviceToDevice));
    }
}


// --- getrf
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, float* A, const int& lda, int* ipiv)
{
    int lwork = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnSgetrf_bufferSize(cusolver_handle, m, n, A, lda, &lwork));

    float* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(float) * lwork));
    unique_cuda_ptr<float> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnSgetrf(cusolver_handle, m, n, A, lda, d_work.get(), ipiv, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }
}
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, double* A, const int& lda, int* ipiv)
{
    int lwork = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnDgetrf_bufferSize(cusolver_handle, m, n, A, lda, &lwork));

    double* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(double) * lwork));
    unique_cuda_ptr<double> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnDgetrf(cusolver_handle, m, n, A, lda, d_work.get(), ipiv, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }
}
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, std::complex<float>* A, const int& lda, int* ipiv)
{
    int lwork = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnCgetrf_bufferSize(cusolver_handle, m, n, reinterpret_cast<cuComplex*>(A), lda, &lwork));

    cuComplex* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(cuComplex) * lwork));
    unique_cuda_ptr<cuComplex> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnCgetrf(cusolver_handle, m, n, reinterpret_cast<cuComplex*>(A), lda, d_work.get(), ipiv, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }
}
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, std::complex<double>* A, const int& lda, int* ipiv)
{
    int lwork = 0;
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnZgetrf_bufferSize(cusolver_handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda, &lwork));

    cuDoubleComplex* raw_d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_work, sizeof(cuDoubleComplex) * lwork));
    unique_cuda_ptr<cuDoubleComplex> d_work(raw_d_work);

    CHECK_CUSOLVER(cusolverDnZgetrf(cusolver_handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda, d_work.get(), ipiv, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }
}

static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, float* A, const int& lda, const int* ipiv, float* B, const int& ldb)
{
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnSgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, A, lda, ipiv, B, ldb, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }
}
static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, double* A, const int& lda, const int* ipiv, double* B, const int& ldb)
{
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnDgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, A, lda, ipiv, B, ldb, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }
}
static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, std::complex<float>* A, const int& lda, const int* ipiv, std::complex<float>* B, const int& ldb)
{
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnCgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, reinterpret_cast<cuComplex*>(A), lda, ipiv, reinterpret_cast<cuComplex*>(B), ldb, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }
}
static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, std::complex<double>* A, const int& lda, const int* ipiv, std::complex<double>* B, const int& ldb)
{
    int h_info = 0;
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&raw_d_info, sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnZgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, reinterpret_cast<cuDoubleComplex*>(A), lda, ipiv, reinterpret_cast<cuDoubleComplex*>(B), ldb, d_info.get()));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }
}

// QR decomposition
// geqrf, orgqr
// Note:
// there are two cusolver geqrf
// one is cusolverDn<t>geqrf
// one is cusolverDnXgeqrf
// which one is better?
//
// template<typename T>
// static inline void geqrf(
//     cusolverDnHandle_t& cusolver_handle,
//     const int64_t m,
//     const int64_t n,
//     T* d_A,           // device matrix A (m x n, column-major)
//     const int64_t lda,
//     T* d_tau         // output: scalar factors of elementary reflectors
// ) {
//     // query workspace size
//     int *d_info = nullptr;    /* error info */
//
//     size_t workspaceInBytesOnDevice = 0; /* size of workspace */
//     void *d_work = nullptr;              /* device workspace */
//     size_t workspaceInBytesOnHost = 0;   /* size of workspace */
//     void *h_work = nullptr;              /* host workspace */
//
//     CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&d_info), sizeof(int)));
//
//     cusolverDnParams_t params = NULL;
//     CHECK_CUSOLVER(cusolverDnCreateParams(&params));
//
//     CHECK_CUSOLVER(cusolverDnXgeqrf_bufferSize(
//         cusolver_handle,
//         params,
//         m, n,
//         traits<T>::cuda_data_type,
//         d_A,
//         lda,
//         traits<T>::cuda_data_type,
//         d_tau,
//         traits<T>::cuda_data_type,
//         &workspaceInBytesOnDevice,
//         &workspaceInBytesOnHost
//     ));
//
//     // allocate device workspace
//     CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&d_work), workspaceInBytesOnDevice));
//
//     // allocate host workspace
//     if (workspaceInBytesOnHost > 0) {
//         h_work = reinterpret_cast<void *>(malloc(workspaceInBytesOnHost));
//         if (h_work == nullptr) {
//             throw std::runtime_error("Error: h_work not allocated.");
//         }
//     }
//
//     // QR factorization
//     CHECK_CUSOLVER(cusolverDnXgeqrf(
//         cusolver_handle,
//         params,
//         m, n,
//         traits<T>::cuda_data_type,
//         d_A,
//         lda,
//         traits<T>::cuda_data_type,
//         d_tau,
//         traits<T>::cuda_data_type,
//         d_work,
//         workspaceInBytesOnDevice,
//         h_work,
//         workspaceInBytesOnHost,
//         d_info
//     ));
//
//     // check info
//     int h_info = 0;
//     CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
//     if (h_info != 0) {
//         // std::printf("%d-th parameter is wrong \n", -info);
//         // print error message
//         std::cout << -h_info << "th parameter is wrong" << std::endl;
//         throw std::runtime_error("geqrf: failed to compute QR decomposition");
//     }
//
//     // clean workspace
//     CHECK_CUDA(cudaFree(d_info));
//     CHECK_CUDA(cudaFree(d_work));
//     if (h_work) free(h_work);
//     CHECK_CUSOLVER(cusolverDnDestroyParams(params));
// }

// geqrf

// --- float ---
static inline void geqrf(
    cusolverDnHandle_t& cusolver_handle,
    const int m,
    const int n,
    float* d_A,
    const int lda,
    float* d_tau
) {
    int lwork = 0;
    CHECK_CUSOLVER(cusolverDnSgeqrf_bufferSize(
        cusolver_handle, m, n, d_A, lda, &lwork));

    unique_cuda_ptr<float> d_work;
    if (lwork > 0) {
        float* raw_d_work = nullptr;
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_work), sizeof(float) * lwork));
        d_work.reset(raw_d_work);
    }
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_info), sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnSgeqrf(
        cusolver_handle, m, n, d_A, lda, d_tau, d_work.get(), lwork, d_info.get()));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("geqrf (S): QR factorization failed");
    }
}

// --- double ---
static inline void geqrf(
    cusolverDnHandle_t& cusolver_handle,
    const int m,
    const int n,
    double* d_A,
    const int lda,
    double* d_tau
) {
    int lwork = 0;
    CHECK_CUSOLVER(cusolverDnDgeqrf_bufferSize(
        cusolver_handle, m, n, d_A, lda, &lwork));

    unique_cuda_ptr<double> d_work;
    if (lwork > 0) {
        double* raw_d_work = nullptr;
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_work), sizeof(double) * lwork));
        d_work.reset(raw_d_work);
    }
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_info), sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnDgeqrf(
        cusolver_handle, m, n, d_A, lda, d_tau, d_work.get(), lwork, d_info.get()));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("geqrf (D): QR factorization failed");
    }
}

// --- std::complex<float> ---
static inline void geqrf(
    cusolverDnHandle_t& cusolver_handle,
    const int m,
    const int n,
    std::complex<float>* d_A,
    const int lda,
    std::complex<float>* d_tau
) {
    int lwork = 0;
    CHECK_CUSOLVER(cusolverDnCgeqrf_bufferSize(
        cusolver_handle, m, n,
        reinterpret_cast<cuComplex*>(d_A),
        lda,
        &lwork
    ));

    unique_cuda_ptr<cuComplex> d_work;
    if (lwork > 0) {
        cuComplex* raw_d_work = nullptr;
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_work), sizeof(cuComplex) * lwork));
        d_work.reset(raw_d_work);
    }
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_info), sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnCgeqrf(
        cusolver_handle, m, n,
        reinterpret_cast<cuComplex*>(d_A),
        lda,
        reinterpret_cast<cuComplex*>(d_tau),
        d_work.get(), lwork, d_info.get()));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("geqrf (C): QR factorization failed");
    }
}

// --- std::complex<double> ---
static inline void geqrf(
    cusolverDnHandle_t& cusolver_handle,
    const int m,
    const int n,
    std::complex<double>* d_A,
    const int lda,
    std::complex<double>* d_tau
) {
    int lwork = 0;
    CHECK_CUSOLVER(cusolverDnZgeqrf_bufferSize(
        cusolver_handle, m, n,
        reinterpret_cast<cuDoubleComplex*>(d_A),
        lda,
        &lwork
    ));

    unique_cuda_ptr<cuDoubleComplex> d_work;
    if (lwork > 0) {
        cuDoubleComplex* raw_d_work = nullptr;
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_work), sizeof(cuDoubleComplex) * lwork));
        d_work.reset(raw_d_work);
    }
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_info), sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnZgeqrf(
        cusolver_handle, m, n,
        reinterpret_cast<cuDoubleComplex*>(d_A),
        lda,
        reinterpret_cast<cuDoubleComplex*>(d_tau),
        d_work.get(), lwork, d_info.get()));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("geqrf (Z): QR factorization failed");
    }
}


// --- float ---
static inline void orgqr(
    cusolverDnHandle_t& cusolver_handle,
    const int m,
    const int n,
    const int k,
    float* d_A,
    const int lda,
    float* d_tau
) {
    int lwork = 0;
    CHECK_CUSOLVER(cusolverDnSorgqr_bufferSize(
        cusolver_handle, m, n, k, d_A, lda, d_tau, &lwork));

    unique_cuda_ptr<float> d_work;
    if (lwork > 0) {
        float* raw_d_work = nullptr;
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_work), sizeof(float) * lwork));
        d_work.reset(raw_d_work);
    }
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_info), sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnSorgqr(
        cusolver_handle, m, n, k, d_A, lda, d_tau, d_work.get(), lwork, d_info.get()));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("orgqr (S): failed to generate Q matrix");
    }
}

// --- double ---
static inline void orgqr(
    cusolverDnHandle_t& cusolver_handle,
    const int m,
    const int n,
    const int k,
    double* d_A,
    const int lda,
    double* d_tau
) {
    int lwork = 0;
    CHECK_CUSOLVER(cusolverDnDorgqr_bufferSize(
        cusolver_handle, m, n, k, d_A, lda, d_tau, &lwork));

    unique_cuda_ptr<double> d_work;
    if (lwork > 0) {
        double* raw_d_work = nullptr;
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_work), sizeof(double) * lwork));
        d_work.reset(raw_d_work);
    }
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_info), sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnDorgqr(
        cusolver_handle, m, n, k, d_A, lda, d_tau, d_work.get(), lwork, d_info.get()));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("orgqr (D): failed to generate Q matrix");
    }
}

// --- std::complex<float> ---
static inline void orgqr(
    cusolverDnHandle_t& cusolver_handle,
    const int m,
    const int n,
    const int k,
    std::complex<float>* d_A,
    const int lda,
    std::complex<float>* d_tau
) {
    int lwork = 0;
    CHECK_CUSOLVER(cusolverDnCungqr_bufferSize(
        cusolver_handle, m, n, k,
        reinterpret_cast<cuComplex*>(d_A),
        lda,
        reinterpret_cast<cuComplex*>(d_tau),
        &lwork));

    unique_cuda_ptr<cuComplex> d_work;
    if (lwork > 0) {
        cuComplex* raw_d_work = nullptr;
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_work), sizeof(cuComplex) * lwork));
        d_work.reset(raw_d_work);
    }
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_info), sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnCungqr(
        cusolver_handle, m, n, k,
        reinterpret_cast<cuComplex*>(d_A),
        lda,
        reinterpret_cast<cuComplex*>(d_tau),
        d_work.get(), lwork, d_info.get()));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("orgqr (C): failed to generate Q matrix");
    }
}

// --- std::complex<double> ---
static inline void orgqr(
    cusolverDnHandle_t& cusolver_handle,
    const int m,
    const int n,
    const int k,
    std::complex<double>* d_A,
    const int lda,
    std::complex<double>* d_tau
) {
    int lwork = 0;
    CHECK_CUSOLVER(cusolverDnZungqr_bufferSize(
        cusolver_handle, m, n, k,
        reinterpret_cast<cuDoubleComplex*>(d_A),
        lda,
        reinterpret_cast<cuDoubleComplex*>(d_tau),
        &lwork));

    unique_cuda_ptr<cuDoubleComplex> d_work;
    if (lwork > 0) {
        cuDoubleComplex* raw_d_work = nullptr;
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_work), sizeof(cuDoubleComplex) * lwork));
        d_work.reset(raw_d_work);
    }
    int* raw_d_info = nullptr;
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&raw_d_info), sizeof(int)));
    unique_cuda_ptr<int> d_info(raw_d_info);

    CHECK_CUSOLVER(cusolverDnZungqr(
        cusolver_handle, m, n, k,
        reinterpret_cast<cuDoubleComplex*>(d_A),
        lda,
        reinterpret_cast<cuDoubleComplex*>(d_tau),
        d_work.get(), lwork, d_info.get()));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info.get(), sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("orgqr (Z): failed to generate Q matrix");
    }
}


} // namespace cuSolverConnector
} // namespace container

#endif // BASE_THIRD_PARTY_CUSOLVER_H_
