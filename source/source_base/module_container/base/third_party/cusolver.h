#ifndef BASE_THIRD_PARTY_CUSOLVER_H_
#define BASE_THIRD_PARTY_CUSOLVER_H_

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <base/macros/cuda.h>

namespace container {
namespace cuSolverConnector {

template <typename T>
static inline
void trtri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, T* A, const int& lda)
{
    size_t d_lwork = 0, h_lwork = 0;
    using Type = typename GetTypeThrust<T>::type;
    cusolverErrcheck(cusolverDnXtrtri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), cublas_diag_type(diag), n, GetTypeCuda<T>::cuda_data_type, reinterpret_cast<Type*>(A), lda, &d_lwork, &h_lwork));
    void* d_work = nullptr, *h_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_work, d_lwork));
    if (h_lwork) {
        h_work = malloc(h_lwork);
        if (h_work == nullptr) {
            throw std::bad_alloc();
        }
    }
    int h_info = 0;
    int* d_info = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));
    // Perform Cholesky decomposition
    cusolverErrcheck(cusolverDnXtrtri(cusolver_handle, cublas_fill_mode(uplo), cublas_diag_type(diag), n, GetTypeCuda<T>::cuda_data_type, reinterpret_cast<Type*>(A), n, d_work, d_lwork, h_work, h_lwork, d_info));
    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("trtri: failed to invert matrix");
    }
    free(h_work);
    cudaErrcheck(cudaFree(d_work));
    cudaErrcheck(cudaFree(d_info));
}

static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, float * A, const int& lda)
{
    int lwork;
    cusolverErrcheck(cusolverDnSpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    float* work;
    cudaErrcheck(cudaMalloc((void**)&work, lwork * sizeof(float)));
    // Perform Cholesky decomposition
    cusolverErrcheck(cusolverDnSpotri(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, nullptr));
    cudaErrcheck(cudaFree(work));
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, double * A, const int& lda)
{
    int lwork;
    cusolverErrcheck(cusolverDnDpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    double* work;
    cudaErrcheck(cudaMalloc((void**)&work, lwork * sizeof(double)));
    // Perform Cholesky decomposition
    cusolverErrcheck(cusolverDnDpotri(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, nullptr));
    cudaErrcheck(cudaFree(work));
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    cusolverErrcheck(cusolverDnCpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex *>(A), n, &lwork));
    cuComplex* work;
    cudaErrcheck(cudaMalloc((void**)&work, lwork * sizeof(cuComplex)));
    // Perform Cholesky decomposition
    cusolverErrcheck(cusolverDnCpotri(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex *>(A), n, work, lwork, nullptr));
    cudaErrcheck(cudaFree(work));
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    cusolverErrcheck(cusolverDnZpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex *>(A), n, &lwork));
    cuDoubleComplex* work;
    cudaErrcheck(cudaMalloc((void**)&work, lwork * sizeof(cuDoubleComplex)));
    // Perform Cholesky decomposition
    cusolverErrcheck(cusolverDnZpotri(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex *>(A), n, work, lwork, nullptr));
    cudaErrcheck(cudaFree(work));
}


static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, float * A, const int& lda)
{
    int lwork;
    int *info = nullptr;
    cudaErrcheck(cudaMalloc((void**)&info, 1 * sizeof(int)));
    cusolverErrcheck(cusolverDnSpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    float* work;
    cudaErrcheck(cudaMalloc((void**)&work, lwork * sizeof(float)));
    // Perform Cholesky decomposition
    cusolverErrcheck(cusolverDnSpotrf(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, info));
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(info));
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, double * A, const int& lda)
{
    int lwork;
    int *info = nullptr;
    cudaErrcheck(cudaMalloc((void**)&info, 1 * sizeof(int)));
    cusolverErrcheck(cusolverDnDpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork));
    double* work;
    cudaErrcheck(cudaMalloc((void**)&work, lwork * sizeof(double)));
    // Perform Cholesky decomposition
    cusolverErrcheck(cusolverDnDpotrf(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, info));
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(info));
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    int *info = nullptr;
    cudaErrcheck(cudaMalloc((void**)&info, 1 * sizeof(int)));
    cusolverErrcheck(cusolverDnCpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex*>(A), lda, &lwork));
    cuComplex* work;
    cudaErrcheck(cudaMalloc((void**)&work, lwork * sizeof(cuComplex)));
    // Perform Cholesky decomposition
    cusolverErrcheck(cusolverDnCpotrf(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex*>(A), lda, work, lwork, info));
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(info));
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    int *info = nullptr;
    cudaErrcheck(cudaMalloc((void**)&info, 1 * sizeof(int)));
    cusolverErrcheck(cusolverDnZpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex*>(A), lda, &lwork));
    cuDoubleComplex* work;
    cudaErrcheck(cudaMalloc((void**)&work, lwork * sizeof(cuDoubleComplex)));
    // Perform Cholesky decomposition
    cusolverErrcheck(cusolverDnZpotrf(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex*>(A), lda, work, lwork, info));
    cudaErrcheck(cudaFree(work));
    cudaErrcheck(cudaFree(info));
}


static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float * W)
{
    // prepare some values for cusolverDnSsyevd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*   d_info = nullptr;
    float* d_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnSsyevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, &lwork));
    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(float) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnSsyevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, d_work, lwork, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    cudaErrcheck(cudaFree(d_info));
    cudaErrcheck(cudaFree(d_work));
}
static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double * W)
{
    // prepare some values for cusolverDnDsyevd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*    d_info = nullptr;
    double* d_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnDsyevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, &lwork));
    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(double) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnDsyevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, d_work, lwork, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    cudaErrcheck(cudaFree(d_info));
    cudaErrcheck(cudaFree(d_work));
}
static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, float * W)
{
    // prepare some values for cusolverDnCheevd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*    d_info = nullptr;
    cuComplex* d_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnCheevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, W, &lwork));
    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnCheevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, W, d_work, lwork, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    cudaErrcheck(cudaFree(d_info));
    cudaErrcheck(cudaFree(d_work));
}
static inline
void heevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, double* W)
{
    // prepare some values for cusolverDnZheevd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*    d_info = nullptr;
    cuDoubleComplex* d_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZheevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, W, &lwork));
    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnZheevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, W, d_work, lwork, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    cudaErrcheck(cudaFree(d_info));
    cudaErrcheck(cudaFree(d_work));
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

    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);
    cusolverEigRange_t range_t = cublas_eig_range(range);

    cusolverErrcheck(cusolverDnSsyevdx_bufferSize(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n, d_A, lda,
        vl, vu, il, iu,
        h_meig,         // ← int* output: number of eigenvalues found
        d_eigen_val,    // ← const float* W (used for query, can be dummy)
        &lwork          // ← int* lwork (output)
    ));

    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(float) * lwork));

    // Main call
    cusolverErrcheck(cusolverDnSsyevdx(
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
    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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

    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);
    cusolverEigRange_t range_t = cublas_eig_range(range);

    cusolverErrcheck(cusolverDnDsyevdx_bufferSize(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n, d_A, lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(double) * lwork));

    cusolverErrcheck(cusolverDnDsyevdx(
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
    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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

    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);
    cusolverEigRange_t range_t = cublas_eig_range(range);

    cusolverErrcheck(cusolverDnCheevdx_bufferSize(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        reinterpret_cast<cuComplex*>(d_A), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork));

    cusolverErrcheck(cusolverDnCheevdx(
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
    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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

    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);
    cusolverEigRange_t range_t = cublas_eig_range(range);

    cusolverErrcheck(cusolverDnZheevdx_bufferSize(
        cusolver_handle,
        jobz_t, range_t, uplo_t,
        n,
        reinterpret_cast<cuDoubleComplex*>(d_A), lda,
        vl, vu, il, iu,
        h_meig,
        d_eigen_val,
        &lwork
    ));

    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));

    cusolverErrcheck(cusolverDnZheevdx(
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
    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnSsygvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, &lwork));
    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(float) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnSsygvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, d_work, lwork, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    cudaErrcheck(cudaFree(d_info));
    cudaErrcheck(cudaFree(d_work));
}
static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double* B, const int& ldb, double * W)
{
    // prepare some values for cusolverDnDsygvd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*   d_info = nullptr;
    double* d_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnDsygvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, &lwork));
    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(double) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnDsygvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, B, ldb, W, d_work, lwork, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    cudaErrcheck(cudaFree(d_info));
    cudaErrcheck(cudaFree(d_work));
}
static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, float* W)
{
    // prepare some values for cusolverDnChegvd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*   d_info = nullptr;
    cuComplex* d_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnChegvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, reinterpret_cast<cuComplex*>(B), ldb, W, &lwork));
    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnChegvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, reinterpret_cast<cuComplex*>(B), ldb, W, d_work, lwork, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    cudaErrcheck(cudaFree(d_info));
    cudaErrcheck(cudaFree(d_work));
}
static inline
void hegvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, double* W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork  = 0;
    int h_info = 0;
    int*   d_info = nullptr;
    cuDoubleComplex* d_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZhegvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, reinterpret_cast<cuDoubleComplex*>(B), ldb, W, &lwork));
    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    cusolverErrcheck(cusolverDnZhegvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, reinterpret_cast<cuDoubleComplex*>(B), ldb, W, d_work, lwork, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("heevd: failed to invert matrix");
    }
    cudaErrcheck(cudaFree(d_info));
    cudaErrcheck(cudaFree(d_work));
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
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // Copy A and B to temporary buffers since sygvdx may modify them
    float *d_A_copy = nullptr, *d_B_copy = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_A_copy, sizeof(float) * n * lda));
    cudaErrcheck(cudaMalloc((void**)&d_B_copy, sizeof(float) * n * lda));
    cudaErrcheck(cudaMemcpy(d_A_copy, d_A, sizeof(float) * n * lda, cudaMemcpyDeviceToDevice));
    cudaErrcheck(cudaMemcpy(d_B_copy, d_B, sizeof(float) * n * lda, cudaMemcpyDeviceToDevice));

    // Set parameters
    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    // Query workspace size
    cusolverErrcheck(cusolverDnSsygvdx_bufferSize(
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
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(float) * lwork));

    // Main call
    cusolverErrcheck(cusolverDnSsygvdx(
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
    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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
        cudaErrcheck(cudaMemcpy(d_eigen_vec, d_A_copy, sizeof(float) * n * m, cudaMemcpyDeviceToDevice));
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

    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    double *d_A_copy = nullptr, *d_B_copy = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_A_copy, sizeof(double) * n * lda));
    cudaErrcheck(cudaMalloc((void**)&d_B_copy, sizeof(double) * n * lda));
    cudaErrcheck(cudaMemcpy(d_A_copy, d_A, sizeof(double) * n * lda, cudaMemcpyDeviceToDevice));
    cudaErrcheck(cudaMemcpy(d_B_copy, d_B, sizeof(double) * n * lda, cudaMemcpyDeviceToDevice));

    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    cusolverErrcheck(cusolverDnDsygvdx_bufferSize(
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

    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(double) * lwork));

    cusolverErrcheck(cusolverDnDsygvdx(
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
    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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
        cudaErrcheck(cudaMemcpy(d_eigen_vec, d_A_copy, sizeof(double) * n * m, cudaMemcpyDeviceToDevice));
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

    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    cuComplex *d_A_copy = nullptr, *d_B_copy = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_A_copy, sizeof(cuComplex) * n * lda));
    cudaErrcheck(cudaMalloc((void**)&d_B_copy, sizeof(cuComplex) * n * lda));
    cudaErrcheck(cudaMemcpy(d_A_copy, reinterpret_cast<cuComplex*>(d_A), sizeof(cuComplex) * n * lda, cudaMemcpyDeviceToDevice));
    cudaErrcheck(cudaMemcpy(d_B_copy, reinterpret_cast<cuComplex*>(d_B), sizeof(cuComplex) * n * lda, cudaMemcpyDeviceToDevice));

    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    cusolverErrcheck(cusolverDnChegvdx_bufferSize(
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

    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork));

    cusolverErrcheck(cusolverDnChegvdx(
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
    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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
        cudaErrcheck(cudaMemcpy(reinterpret_cast<cuComplex*>(d_eigen_vec), d_A_copy, sizeof(cuComplex) * n * m, cudaMemcpyDeviceToDevice));
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

    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    cuDoubleComplex *d_A_copy = nullptr, *d_B_copy = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_A_copy, sizeof(cuDoubleComplex) * n * lda));
    cudaErrcheck(cudaMalloc((void**)&d_B_copy, sizeof(cuDoubleComplex) * n * lda));
    cudaErrcheck(cudaMemcpy(d_A_copy, reinterpret_cast<cuDoubleComplex*>(d_A), sizeof(cuDoubleComplex) * n * lda, cudaMemcpyDeviceToDevice));
    cudaErrcheck(cudaMemcpy(d_B_copy, reinterpret_cast<cuDoubleComplex*>(d_B), sizeof(cuDoubleComplex) * n * lda, cudaMemcpyDeviceToDevice));

    cusolverEigType_t itype_t = cublas_eig_type(itype);
    cusolverEigMode_t jobz_t = cublas_eig_mode(jobz);
    cusolverEigRange_t range_t = cublas_eig_range(range);
    cublasFillMode_t uplo_t = cublas_fill_mode(uplo);

    cusolverErrcheck(cusolverDnZhegvdx_bufferSize(
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

    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));

    cusolverErrcheck(cusolverDnZhegvdx(
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
    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
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
        cudaErrcheck(cudaMemcpy(reinterpret_cast<cuDoubleComplex*>(d_eigen_vec), d_A_copy, sizeof(cuDoubleComplex) * n * m, cudaMemcpyDeviceToDevice));
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
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnSgetrf_bufferSize(cusolver_handle, m, n, A, lda, &lwork));

    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(float) * lwork));

    // Perform LU decomposition
    cusolverErrcheck(cusolverDnSgetrf(cusolver_handle, m, n, A, lda, d_work, ipiv, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }

    cudaErrcheck(cudaFree(d_work));
    cudaErrcheck(cudaFree(d_info));
}
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, double* A, const int& lda, int* ipiv)
{
    // prepare some values for cusolverDnDgetrf_bufferSize
    int lwork = 0;
    int h_info = 0;
    int* d_info = nullptr;
    double* d_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnDgetrf_bufferSize(cusolver_handle, m, n, A, lda, &lwork));

    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(double) * lwork));

    // Perform LU decomposition
    cusolverErrcheck(cusolverDnDgetrf(cusolver_handle, m, n, A, lda, d_work, ipiv, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }

    cudaErrcheck(cudaFree(d_work));
    cudaErrcheck(cudaFree(d_info));
}
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, std::complex<float>* A, const int& lda, int* ipiv)
{
    // prepare some values for cusolverDnCgetrf_bufferSize
    int lwork = 0;
    int h_info = 0;
    int* d_info = nullptr;
    cuComplex* d_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnCgetrf_bufferSize(cusolver_handle, m, n, reinterpret_cast<cuComplex*>(A), lda, &lwork));

    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork));

    // Perform LU decomposition
    cusolverErrcheck(cusolverDnCgetrf(cusolver_handle, m, n, reinterpret_cast<cuComplex*>(A), lda, d_work, ipiv, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }

    cudaErrcheck(cudaFree(d_work));
    cudaErrcheck(cudaFree(d_info));
}
static inline
void getrf(cusolverDnHandle_t& cusolver_handle, const int& m, const int& n, std::complex<double>* A, const int& lda, int* ipiv)
{
    // prepare some values for cusolverDnZgetrf_bufferSize
    int lwork = 0;
    int h_info = 0;
    int* d_info = nullptr;
    cuDoubleComplex* d_work = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverErrcheck(cusolverDnZgetrf_bufferSize(cusolver_handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda, &lwork));

    // allocate memory
    cudaErrcheck(cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork));

    // Perform LU decomposition
    cusolverErrcheck(cusolverDnZgetrf(cusolver_handle, m, n, reinterpret_cast<cuDoubleComplex*>(A), lda, d_work, ipiv, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrf: failed to compute LU factorization");
    }

    cudaErrcheck(cudaFree(d_work));
    cudaErrcheck(cudaFree(d_info));
}

static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, float* A, const int& lda, const int* ipiv, float* B, const int& ldb)
{
    int h_info = 0;
    int* d_info = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    cusolverErrcheck(cusolverDnSgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, A, lda, ipiv, B, ldb, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }

    cudaErrcheck(cudaFree(d_info));
}
static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, double* A, const int& lda, const int* ipiv, double* B, const int& ldb)
{
    int h_info = 0;
    int* d_info = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    cusolverErrcheck(cusolverDnDgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, A, lda, ipiv, B, ldb, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }

    cudaErrcheck(cudaFree(d_info));
}
static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, std::complex<float>* A, const int& lda, const int* ipiv, std::complex<float>* B, const int& ldb)
{
    int h_info = 0;
    int* d_info = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    cusolverErrcheck(cusolverDnCgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, reinterpret_cast<cuComplex*>(A), lda, ipiv, reinterpret_cast<cuComplex*>(B), ldb, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }

    cudaErrcheck(cudaFree(d_info));
}
static inline
void getrs(cusolverDnHandle_t& cusolver_handle, const char& trans, const int& n, const int& nrhs, std::complex<double>* A, const int& lda, const int* ipiv, std::complex<double>* B, const int& ldb)
{
    int h_info = 0;
    int* d_info = nullptr;
    cudaErrcheck(cudaMalloc((void**)&d_info, sizeof(int)));

    cusolverErrcheck(cusolverDnZgetrs(cusolver_handle, GetCublasOperation(trans), n, nrhs, reinterpret_cast<cuDoubleComplex*>(A), lda, ipiv, reinterpret_cast<cuDoubleComplex*>(B), ldb, d_info));

    cudaErrcheck(cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("getrs: failed to solve the linear system");
    }

    cudaErrcheck(cudaFree(d_info));
}

} // namespace cuSolverConnector
} // namespace container

#endif // BASE_THIRD_PARTY_CUSOLVER_H_
