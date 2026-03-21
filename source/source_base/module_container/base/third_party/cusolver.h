#ifndef BASE_THIRD_PARTY_CUSOLVER_H_
#define BASE_THIRD_PARTY_CUSOLVER_H_

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <iostream>

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

template <typename T>
static inline
void trtri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, T* A, const int& lda)
{
    size_t d_lwork = 0, h_lwork = 0;
    using Type = typename GetTypeThrust<T>::type;
    CHECK_CUSOLVER(cusolverDnXtrtri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), cublas_diag_type(diag), n, GetTypeCuda<T>::cuda_data_type, reinterpret_cast<Type*>(A), lda, &d_lwork, &h_lwork));
    void* d_work = nullptr, *h_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_work, d_lwork));
    if (h_lwork) {
        h_work = malloc(h_lwork);
        if (h_work == nullptr) {
            throw std::bad_alloc();
        }
    }
    int h_info = 0;
    int* d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));
    // Perform Cholesky decomposition
    CHECK_CUSOLVER(cusolverDnXtrtri(cusolver_handle, cublas_fill_mode(uplo), cublas_diag_type(diag), n, GetTypeCuda<T>::cuda_data_type, reinterpret_cast<Type*>(A), n, d_work, d_lwork, h_work, h_lwork, d_info));
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("trtri: failed to invert matrix");
    }
    free(h_work);
    CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
}

static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, float * A, const int& lda)
{
    int lwork;
    CHECK_CUSOLVER(cusolverDnSpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    float* work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&work, lwork * sizeof(float)));
    // Perform Cholesky decomposition
    CHECK_CUSOLVER(cusolverDnSpotri(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, nullptr));
    CHECK_CUDA(cudaFree(work));
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, double * A, const int& lda)
{
    int lwork;
    CHECK_CUSOLVER(cusolverDnDpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    double* work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&work, lwork * sizeof(double)));
    // Perform Cholesky decomposition
    CHECK_CUSOLVER(cusolverDnDpotri(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, nullptr));
    CHECK_CUDA(cudaFree(work));
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    CHECK_CUSOLVER(cusolverDnCpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex *>(A), n, &lwork));
    cuComplex* work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&work, lwork * sizeof(cuComplex)));
    // Perform Cholesky decomposition
    CHECK_CUSOLVER(cusolverDnCpotri(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex *>(A), n, work, lwork, nullptr));
    CHECK_CUDA(cudaFree(work));
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    CHECK_CUSOLVER(cusolverDnZpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex *>(A), n, &lwork));
    cuDoubleComplex* work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&work, lwork * sizeof(cuDoubleComplex)));
    // Perform Cholesky decomposition
    CHECK_CUSOLVER(cusolverDnZpotri(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex *>(A), n, work, lwork, nullptr));
    CHECK_CUDA(cudaFree(work));
}


static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, float * A, const int& lda)
{
    int lwork;
    int *info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&info, 1 * sizeof(int)));
    CHECK_CUSOLVER(cusolverDnSpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    float* work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&work, lwork * sizeof(float)));
    // Perform Cholesky decomposition
    CHECK_CUSOLVER(cusolverDnSpotrf(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, info));
    CHECK_CUDA(cudaFree(work));
    CHECK_CUDA(cudaFree(info));
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, double * A, const int& lda)
{
    int lwork;
    int *info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&info, 1 * sizeof(int)));
    CHECK_CUSOLVER(cusolverDnDpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    double* work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&work, lwork * sizeof(double)));
    // Perform Cholesky decomposition
    CHECK_CUSOLVER(cusolverDnDpotrf(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, info));
    CHECK_CUDA(cudaFree(work));
    CHECK_CUDA(cudaFree(info));
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    int *info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&info, 1 * sizeof(int)));
    CHECK_CUSOLVER(cusolverDnCpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex*>(A), lda, &lwork));
    cuComplex* work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&work, lwork * sizeof(cuComplex)));
    // Perform Cholesky decomposition
    CHECK_CUSOLVER(cusolverDnCpotrf(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex*>(A), lda, work, lwork, info));
    CHECK_CUDA(cudaFree(work));
    CHECK_CUDA(cudaFree(info));
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    int *info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&info, 1 * sizeof(int)));
    CHECK_CUSOLVER(cusolverDnZpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex*>(A), lda, &lwork));
    cuDoubleComplex* work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&work, lwork * sizeof(cuDoubleComplex)));
    // Perform Cholesky decomposition
    CHECK_CUSOLVER(cusolverDnZpotrf(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex*>(A), lda, work, lwork, info));
    CHECK_CUDA(cudaFree(work));
    CHECK_CUDA(cudaFree(info));
}


static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float * W)
{
    // prepare some values for cusolverDnSsyevd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*   d_info = nullptr;
    float* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnSsyevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, &lwork));
    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(float) * lwork));
    // compute eigenvalues and eigenvectors.
    CHECK_CUSOLVER(cusolverDnSsyevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, d_work, lwork, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    CHECK_CUDA(cudaFree(d_info));
    CHECK_CUDA(cudaFree(d_work));
}
static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double * W)
{
    // prepare some values for cusolverDnDsyevd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*    d_info = nullptr;
    double* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnDsyevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, &lwork));
    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(double) * lwork));
    // compute eigenvalues and eigenvectors.
    CHECK_CUSOLVER(cusolverDnDsyevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, d_work, lwork, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    CHECK_CUDA(cudaFree(d_info));
    CHECK_CUDA(cudaFree(d_work));
}
static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, float * W)
{
    // prepare some values for cusolverDnCheevd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*    d_info = nullptr;
    cuComplex* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnCheevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, W, &lwork));
    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    CHECK_CUSOLVER(cusolverDnCheevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, W, d_work, lwork, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    CHECK_CUDA(cudaFree(d_info));
    CHECK_CUDA(cudaFree(d_work));
}
static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, double* W)
{
    // prepare some values for cusolverDnZheevd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*    d_info = nullptr;
    cuDoubleComplex* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnZheevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, W, &lwork));
    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    CHECK_CUSOLVER(cusolverDnZheevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, W, d_work, lwork, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    CHECK_CUDA(cudaFree(d_info));
    CHECK_CUDA(cudaFree(d_work));
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
    int* d_info = nullptr;
    float* d_work = nullptr;

    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);
    cusolverEigRange_t range_t = cublas_eig_range(range);

    CHECK_CUSOLVER(cusolverDnSsyevdx_bufferSize(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n, d_A, lda,
        vl, vu, il, iu,
        h_meig,         // ← int* output: number of eigenvalues found
        d_eigen_val,    // ← const float* W (used for query, can be dummy)
        &lwork          // ← int* lwork (output)
    ));

    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(float) * lwork));

    // Main call
    CHECK_CUSOLVER(cusolverDnSsyevdx(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        d_A, lda,
        vl, vu, il, iu,
        h_meig,         // ← int* output
        d_eigen_val,    // ← float* W: eigenvalues
        d_work, lwork,
        d_info
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        cudaFree(d_info); cudaFree(d_work);
        throw std::runtime_error("heevdx (float) failed with info = " + std::to_string(h_info));
    }

    cudaFree(d_info);
    cudaFree(d_work);
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
    int* d_info = nullptr;
    double* d_work = nullptr;

    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

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

    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(double) * lwork));

    CHECK_CUSOLVER(cusolverDnDsyevdx(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        d_A, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work, lwork,
        d_info
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        cudaFree(d_info); cudaFree(d_work);
        throw std::runtime_error("heevdx (double) failed with info = " + std::to_string(h_info));
    }

    cudaFree(d_info);
    cudaFree(d_work);
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
    int* d_info = nullptr;
    cuComplex* d_work = nullptr;

    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

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

    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork));

    CHECK_CUSOLVER(cusolverDnCheevdx(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        reinterpret_cast<cuComplex*>(d_A), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work, lwork,
        d_info
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        cudaFree(d_info); cudaFree(d_work);
        throw std::runtime_error("heevdx (complex<float>) failed with info = " + std::to_string(h_info));
    }

    cudaFree(d_info);
    cudaFree(d_work);
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
    int* d_info = nullptr;
    cuDoubleComplex* d_work = nullptr;

    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

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

    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));

    CHECK_CUSOLVER(cusolverDnZheevdx(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        reinterpret_cast<cuDoubleComplex*>(d_A), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work, lwork,
        d_info
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        cudaFree(d_info); cudaFree(d_work);
        throw std::runtime_error("heevdx (complex<double>) failed with info = " + std::to_string(h_info));
    }

    cudaFree(d_info);
    cudaFree(d_work);
}

static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float* B, const int& ldb, float * W)
{
    // prepare some values for cusolverDnSsygvd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*   d_info = nullptr;
    float* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnSsygvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, &lwork));
    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(float) * lwork));
    // compute eigenvalues and eigenvectors.
    CHECK_CUSOLVER(cusolverDnSsygvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, d_work, lwork, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    CHECK_CUDA(cudaFree(d_info));
    CHECK_CUDA(cudaFree(d_work));
}
static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double* B, const int& ldb, double * W)
{
    // prepare some values for cusolverDnDsygvd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*   d_info = nullptr;
    double* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnDsygvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, &lwork));
    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(double) * lwork));
    // compute eigenvalues and eigenvectors.
    CHECK_CUSOLVER(cusolverDnDsygvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, d_work, lwork, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    CHECK_CUDA(cudaFree(d_info));
    CHECK_CUDA(cudaFree(d_work));
}
static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, float* W)
{
    // prepare some values for cusolverDnChegvd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*   d_info = nullptr;
    cuComplex* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnChegvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, reinterpret_cast<cuComplex*>(B), ldb, W, &lwork));
    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    CHECK_CUSOLVER(cusolverDnChegvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, reinterpret_cast<cuComplex*>(B), ldb, W, d_work, lwork, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    CHECK_CUDA(cudaFree(d_info));
    CHECK_CUDA(cudaFree(d_work));
}
static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, double* W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*   d_info = nullptr;
    cuDoubleComplex* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnZhegvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, reinterpret_cast<cuDoubleComplex*>(B), ldb, W, &lwork));
    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    CHECK_CUSOLVER(cusolverDnZhegvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, reinterpret_cast<cuDoubleComplex*>(B), ldb, W, d_work, lwork, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    CHECK_CUDA(cudaFree(d_info));
    CHECK_CUDA(cudaFree(d_work));
}

// =====================================================================================================
// hegvd x: Compute selected eigenvalues and eigenvectors of generalized Hermitian-definite eigenproblem
//          A * x = lambda * B * x
// =====================================================================================================

// --- float ---
static inline
void hegvdx(
    cusolverDnHandle_t& cusolver_handle,
    const int itype,           // 1: A*x = lambda*B*x
    const char jobz,           // 'V' or 'N'
    const char range,          // 'I', 'V', 'A'
    const char uplo,           // 'U' or 'L'
    const int n,
    const int lda,
    float* d_A,                // Input matrix A (device)
    float* d_B,                // Input matrix B (device)
    const float vl,            // for RANGE='V'
    const float vu,
    const int il,              // for RANGE='I'
    const int iu,
    int* h_meig,               // output: number of eigenvalues found
    float* d_eigen_val,        // output: eigenvalues
    float* d_eigen_vec         // output: eigenvectors (if jobz='V'), size ldz × m
) {
    int lwork = 0;
    int *d_info = nullptr;
    float *d_work = nullptr;

    // Allocate device info
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // Copy A and B to temporary buffers since sygvdx may modify them
    float *d_A_copy = nullptr, *d_B_copy = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_A_copy, sizeof(float) * n * lda));
    CHECK_CUDA(cudaMalloc((void**)&d_B_copy, sizeof(float) * n * lda));
    CHECK_CUDA(cudaMemcpy(d_A_copy, d_A, sizeof(float) * n * lda, cudaMemcpyDeviceToDevice));
    CHECK_CUDA(cudaMemcpy(d_B_copy, d_B, sizeof(float) * n * lda, cudaMemcpyDeviceToDevice));

    // Set parameters
    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    // Query workspace size
    CHECK_CUSOLVER(cusolverDnSsygvdx_bufferSize(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy, lda,
        d_B_copy, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    // Allocate workspace
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(float) * lwork));

    // Main call
    CHECK_CUSOLVER(cusolverDnSsygvdx(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy, lda,
        d_B_copy, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work, lwork,
        d_info
    ));

    // Check result
    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info < 0) {
        throw std::runtime_error("hegvdx (float): illegal argument #" + std::to_string(-h_info));
    } else if (h_info > 0) {
        // If h_info <= n: convergence issue in tridiag solver (no vec) OR
        // If h_info > n: B's leading minor of order (h_info - n) is not positive definite
        if (jobz_t == CUSOLVER_EIG_MODE_NOVECTOR && h_info <= n) {
            throw std::runtime_error("hegvdx (float): failed to converge, " + std::to_string(h_info) + " off-diagonal elements didn't converge");
        } else if (h_info > n) {
            throw std::runtime_error("hegvdx (float): leading minor of order " + std::to_string(h_info - n) + " of B is not positive definite");
        }
    }

    // If jobz == 'V', copy eigenvectors from A (which now contains Z) to output
    if (jobz == 'V') {
        const int m = (*h_meig); // number of eigenvectors computed
        CHECK_CUDA(cudaMemcpy(d_eigen_vec, d_A_copy, sizeof(float) * n * m, cudaMemcpyDeviceToDevice));
    }

    // Cleanup
    cudaFree(d_info);
    cudaFree(d_work);
    cudaFree(d_A_copy);
    cudaFree(d_B_copy);
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
    int *d_info = nullptr;
    double *d_work = nullptr;

    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    double *d_A_copy = nullptr, *d_B_copy = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_A_copy, sizeof(double) * n * lda));
    CHECK_CUDA(cudaMalloc((void**)&d_B_copy, sizeof(double) * n * lda));
    CHECK_CUDA(cudaMemcpy(d_A_copy, d_A, sizeof(double) * n * lda, cudaMemcpyDeviceToDevice));
    CHECK_CUDA(cudaMemcpy(d_B_copy, d_B, sizeof(double) * n * lda, cudaMemcpyDeviceToDevice));

    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    CHECK_CUSOLVER(cusolverDnDsygvdx_bufferSize(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy, lda,
        d_B_copy, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(double) * lwork));

    CHECK_CUSOLVER(cusolverDnDsygvdx(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy, lda,
        d_B_copy, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work, lwork,
        d_info
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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
        CHECK_CUDA(cudaMemcpy(d_eigen_vec, d_A_copy, sizeof(double) * n * m, cudaMemcpyDeviceToDevice));
    }

    cudaFree(d_info);
    cudaFree(d_work);
    cudaFree(d_A_copy);
    cudaFree(d_B_copy);
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
    int *d_info = nullptr;
    cuComplex *d_work = nullptr;

    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    cuComplex *d_A_copy = nullptr, *d_B_copy = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_A_copy, sizeof(cuComplex) * n * lda));
    CHECK_CUDA(cudaMalloc((void**)&d_B_copy, sizeof(cuComplex) * n * lda));
    CHECK_CUDA(cudaMemcpy(d_A_copy, reinterpret_cast<cuComplex*>(d_A), sizeof(cuComplex) * n * lda, cudaMemcpyDeviceToDevice));
    CHECK_CUDA(cudaMemcpy(d_B_copy, reinterpret_cast<cuComplex*>(d_B), sizeof(cuComplex) * n * lda, cudaMemcpyDeviceToDevice));

    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    CHECK_CUSOLVER(cusolverDnChegvdx_bufferSize(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy, lda,
        d_B_copy, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork));

    CHECK_CUSOLVER(cusolverDnChegvdx(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy, lda,
        d_B_copy, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work, lwork,
        d_info
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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
        CHECK_CUDA(cudaMemcpy(reinterpret_cast<cuComplex*>(d_eigen_vec), d_A_copy, sizeof(cuComplex) * n * m, cudaMemcpyDeviceToDevice));
    }

    cudaFree(d_info);
    cudaFree(d_work);
    cudaFree(d_A_copy);
    cudaFree(d_B_copy);
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
    int *d_info = nullptr;
    cuDoubleComplex *d_work = nullptr;

    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    cuDoubleComplex *d_A_copy = nullptr, *d_B_copy = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_A_copy, sizeof(cuDoubleComplex) * n * lda));
    CHECK_CUDA(cudaMalloc((void**)&d_B_copy, sizeof(cuDoubleComplex) * n * lda));
    CHECK_CUDA(cudaMemcpy(d_A_copy, reinterpret_cast<cuDoubleComplex*>(d_A), sizeof(cuDoubleComplex) * n * lda, cudaMemcpyDeviceToDevice));
    CHECK_CUDA(cudaMemcpy(d_B_copy, reinterpret_cast<cuDoubleComplex*>(d_B), sizeof(cuDoubleComplex) * n * lda, cudaMemcpyDeviceToDevice));

    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    CHECK_CUSOLVER(cusolverDnZhegvdx_bufferSize(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy, lda,
        d_B_copy, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));

    CHECK_CUSOLVER(cusolverDnZhegvdx(
        cusolver_handle,
        itype_t, jobz_t, range_t, uplo_t,
        n,
        d_A_copy, lda,
        d_B_copy, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        d_work, lwork,
        d_info
    ));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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
        CHECK_CUDA(cudaMemcpy(reinterpret_cast<cuDoubleComplex*>(d_eigen_vec), d_A_copy, sizeof(cuDoubleComplex) * n * m, cudaMemcpyDeviceToDevice));
    }

    cudaFree(d_info);
    cudaFree(d_work);
    cudaFree(d_A_copy);
    cudaFree(d_B_copy);
}


// --- getrf
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, float* A, const int& lda, int* ipiv)
{
    // prepare some values for cusolverDnSgetrf_bufferSize
    int lwork = 0;
    int h_info = 0;
    int* d_info = nullptr;
    float* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnSgetrf_bufferSize(cusolver_handle, m, n, A, lda, &lwork));

    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(float) * lwork));

    // Perform LU decomposition
    CHECK_CUSOLVER(cusolverDnSgetrf(cusolver_handle, m, n, A, lda, d_work, ipiv, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }

    CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
}
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, double* A, const int& lda, int* ipiv)
{
    // prepare some values for cusolverDnDgetrf_bufferSize
    int lwork = 0;
    int h_info = 0;
    int* d_info = nullptr;
    double* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnDgetrf_bufferSize(cusolver_handle, m, n, A, lda, &lwork));

    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(double) * lwork));

    // Perform LU decomposition
    CHECK_CUSOLVER(cusolverDnDgetrf(cusolver_handle, m, n, A, lda, d_work, ipiv, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }

    CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
}
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, std::complex<float>* A, const int& lda, int* ipiv)
{
    // prepare some values for cusolverDnCgetrf_bufferSize
    int lwork = 0;
    int h_info = 0;
    int* d_info = nullptr;
    cuComplex* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnCgetrf_bufferSize(cusolver_handle, m, n, reinterpret_cast<cuComplex*>(A), lda, &lwork));

    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork));

    // Perform LU decomposition
    CHECK_CUSOLVER(cusolverDnCgetrf(cusolver_handle, m, n, reinterpret_cast<cuComplex*>(A), lda, d_work, ipiv, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }

    CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
}
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, std::complex<double>* A, const int& lda, int* ipiv)
{
    // prepare some values for cusolverDnZgetrf_bufferSize
    int lwork = 0;
    int h_info = 0;
    int* d_info = nullptr;
    cuDoubleComplex* d_work = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    CHECK_CUSOLVER(cusolverDnZgetrf_bufferSize(cusolver_handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda, &lwork));

    // allocate memory
    CHECK_CUDA(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));

    // Perform LU decomposition
    CHECK_CUSOLVER(cusolverDnZgetrf(cusolver_handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda, d_work, ipiv, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }

    CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
}

static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, float* A, const int& lda, const int* ipiv, float* B, const int& ldb)
{
    int h_info = 0;
    int* d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    CHECK_CUSOLVER(cusolverDnSgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, A, lda, ipiv, B, ldb, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }

    CHECK_CUDA(cudaFree(d_info));
}
static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, double* A, const int& lda, const int* ipiv, double* B, const int& ldb)
{
    int h_info = 0;
    int* d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    CHECK_CUSOLVER(cusolverDnDgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, A, lda, ipiv, B, ldb, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }

    CHECK_CUDA(cudaFree(d_info));
}
static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, std::complex<float>* A, const int& lda, const int* ipiv, std::complex<float>* B, const int& ldb)
{
    int h_info = 0;
    int* d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    CHECK_CUSOLVER(cusolverDnCgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, reinterpret_cast<cuComplex*>(A), lda, ipiv, reinterpret_cast<cuComplex*>(B), ldb, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }

    CHECK_CUDA(cudaFree(d_info));
}
static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, std::complex<double>* A, const int& lda, const int* ipiv, std::complex<double>* B, const int& ldb)
{
    int h_info = 0;
    int* d_info = nullptr;
    CHECK_CUDA(cudaMalloc((void**)&d_info, sizeof(int)));

    CHECK_CUSOLVER(cusolverDnZgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, reinterpret_cast<cuDoubleComplex*>(A), lda, ipiv, reinterpret_cast<cuDoubleComplex*>(B), ldb, d_info));

    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }

    CHECK_CUDA(cudaFree(d_info));
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

    float* d_work = nullptr;
    int*   d_info = nullptr;

    if (lwork > 0) {
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(float) * lwork));
    }
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(int)));

    CHECK_CUSOLVER(cusolverDnSgeqrf(
        cusolver_handle, m, n, d_A, lda, d_tau, d_work, lwork, d_info));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        std::cout << "geqrf (S): info = " << h_info << std::endl;
        if (d_work) CHECK_CUDA(cudaFree(d_work));
        CHECK_CUDA(cudaFree(d_info));
        throw std::runtime_error("geqrf (S): QR factorization failed");
    }

    if (d_work) CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
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

    double* d_work = nullptr;
    int*    d_info = nullptr;

    if (lwork > 0) {
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(double) * lwork));
    }
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(int)));

    CHECK_CUSOLVER(cusolverDnDgeqrf(
        cusolver_handle, m, n, d_A, lda, d_tau, d_work, lwork, d_info));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        std::cout << "geqrf (D): info = " << h_info << std::endl;
        if (d_work) CHECK_CUDA(cudaFree(d_work));
        CHECK_CUDA(cudaFree(d_info));
        throw std::runtime_error("geqrf (D): QR factorization failed");
    }

    if (d_work) CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
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
        &lwork  // ← 这里才是 lwork 的地址！
    ));

    cuComplex* d_work = nullptr;
    int*       d_info = nullptr;

    if (lwork > 0) {
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(cuComplex) * lwork));
    }
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(int)));

    CHECK_CUSOLVER(cusolverDnCgeqrf(
        cusolver_handle, m, n,
        reinterpret_cast<cuComplex*>(d_A),
        lda,
        reinterpret_cast<cuComplex*>(d_tau),  // ← 这里才是 d_tau
        d_work, lwork, d_info));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        std::cout << "geqrf (C): info = " << h_info << std::endl;
        if (d_work) CHECK_CUDA(cudaFree(d_work));
        CHECK_CUDA(cudaFree(d_info));
        throw std::runtime_error("geqrf (C): QR factorization failed");
    }

    if (d_work) CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
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

    cuDoubleComplex* d_work = nullptr;
    int*             d_info = nullptr;

    if (lwork > 0) {
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(cuDoubleComplex) * lwork));
    }
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(int)));

    CHECK_CUSOLVER(cusolverDnZgeqrf(
        cusolver_handle, m, n,
        reinterpret_cast<cuDoubleComplex*>(d_A),
        lda,
        reinterpret_cast<cuDoubleComplex*>(d_tau),
        d_work, lwork, d_info));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        std::cout << "geqrf (Z): info = " << h_info << std::endl;
        if (d_work) CHECK_CUDA(cudaFree(d_work));
        CHECK_CUDA(cudaFree(d_info));
        throw std::runtime_error("geqrf (Z): QR factorization failed");
    }

    if (d_work) CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
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

    float* d_work = nullptr;
    int*   d_info = nullptr;

    if (lwork > 0) {
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(float) * lwork));
    }
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(int)));

    CHECK_CUSOLVER(cusolverDnSorgqr(
        cusolver_handle, m, n, k, d_A, lda, d_tau, d_work, lwork, d_info));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        std::cout << "orgqr (S): info = " << h_info << " (failure at parameter " << -h_info << ")" << std::endl;
        if (d_work) CHECK_CUDA(cudaFree(d_work));
        CHECK_CUDA(cudaFree(d_info));
        throw std::runtime_error("orgqr (S): failed to generate Q matrix");
    }

    // clean workspace
    if (d_work) CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
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

    double* d_work = nullptr;
    int*    d_info = nullptr;

    if (lwork > 0) {
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(double) * lwork));
    }
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(int)));

    CHECK_CUSOLVER(cusolverDnDorgqr(
        cusolver_handle, m, n, k, d_A, lda, d_tau, d_work, lwork, d_info));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        std::cout << "orgqr (D): info = " << h_info << std::endl;
        if (d_work) CHECK_CUDA(cudaFree(d_work));
        CHECK_CUDA(cudaFree(d_info));
        throw std::runtime_error("orgqr (D): failed to generate Q matrix");
    }

    if (d_work) CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
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

    cuComplex* d_work = nullptr;
    int*       d_info = nullptr;

    if (lwork > 0) {
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(cuComplex) * lwork));
    }
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(int)));

    CHECK_CUSOLVER(cusolverDnCungqr(
        cusolver_handle, m, n, k,
        reinterpret_cast<cuComplex*>(d_A),
        lda,
        reinterpret_cast<cuComplex*>(d_tau),
        d_work, lwork, d_info));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        std::cout << "orgqr (C): info = " << h_info << std::endl;
        if (d_work) CHECK_CUDA(cudaFree(d_work));
        CHECK_CUDA(cudaFree(d_info));
        throw std::runtime_error("orgqr (C): failed to generate Q matrix");
    }

    if (d_work) CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
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

    cuDoubleComplex* d_work = nullptr;
    int*             d_info = nullptr;

    if (lwork > 0) {
        CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_work), sizeof(cuDoubleComplex) * lwork));
    }
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void**>(&d_info), sizeof(int)));

    CHECK_CUSOLVER(cusolverDnZungqr(
        cusolver_handle, m, n, k,
        reinterpret_cast<cuDoubleComplex*>(d_A),
        lda,
        reinterpret_cast<cuDoubleComplex*>(d_tau),
        d_work, lwork, d_info));

    int h_info = 0;
    CHECK_CUDA(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        std::cout << "orgqr (Z): info = " << h_info << std::endl;
        if (d_work) CHECK_CUDA(cudaFree(d_work));
        CHECK_CUDA(cudaFree(d_info));
        throw std::runtime_error("orgqr (Z): failed to generate Q matrix");
    }

    if (d_work) CHECK_CUDA(cudaFree(d_work));
    CHECK_CUDA(cudaFree(d_info));
}


} // namespace cuSolverConnector
} // namespace container

#endif // BASE_THIRD_PARTY_CUSOLVER_H_
