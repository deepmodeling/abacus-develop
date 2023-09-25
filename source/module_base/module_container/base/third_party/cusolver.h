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
    cusolverDnXtrtri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), cublas_diag_type(diag), n, GetTypeCuda<T>::cuda_data_type, reinterpret_cast<Type*>(A), lda, &d_lwork, &h_lwork);
    void* d_work = nullptr, *h_work = nullptr;
    cudaMalloc((void**)&d_work, d_lwork);
    if (h_lwork) {
        h_work = malloc(h_lwork);
        if (h_work == nullptr) {
            throw std::bad_alloc();
        }
    }
    int h_info = 0;
    int* d_info = nullptr;
    cudaMalloc((void**)&d_info, sizeof(int));
    // Perform Cholesky decomposition
    cusolverDnXtrtri(cusolver_handle, cublas_fill_mode(uplo), cublas_diag_type(diag), n, GetTypeCuda<T>::cuda_data_type, reinterpret_cast<Type*>(A), n, d_work, d_lwork, h_work, h_lwork, d_info);
    cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("trtri: failed to invert matrix");
    }
    free(h_work);
    cudaFree(d_work);
    cudaFree(d_info);
}

static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, float * A, const int& lda)
{
    int lwork;
    cusolverDnSpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork);
    float* work;
    cudaMalloc((void**)&work, lwork * sizeof(float));
    // Perform Cholesky decomposition
    cusolverDnSpotri(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, nullptr);
    cudaFree(work);
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, double * A, const int& lda)
{
    int lwork;
    cusolverDnDpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork);
    double* work;
    cudaMalloc((void**)&work, lwork * sizeof(double));
    // Perform Cholesky decomposition
    cusolverDnDpotri(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, nullptr);
    cudaFree(work);
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    cusolverDnCpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex *>(A), n, &lwork);
    cuComplex* work;
    cudaMalloc((void**)&work, lwork * sizeof(cuComplex));
    // Perform Cholesky decomposition
    cusolverDnCpotri(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex *>(A), n, work, lwork, nullptr);
    cudaFree(work);
}
static inline
void potri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    cusolverDnZpotri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex *>(A), n, &lwork);
    cuDoubleComplex* work;
    cudaMalloc((void**)&work, lwork * sizeof(cuDoubleComplex));
    // Perform Cholesky decomposition
    cusolverDnZpotri(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex *>(A), n, work, lwork, nullptr);
    cudaFree(work);
}


static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, float * A, const int& lda)
{
    int lwork;
    cusolverDnSpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork);
    float* work;
    cudaMalloc((void**)&work, lwork * sizeof(float));
    // Perform Cholesky decomposition
    cusolverDnSpotrf(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, nullptr);
    cudaFree(work);
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, double * A, const int& lda)
{
    int lwork;
    cusolverDnDpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, A, n, &lwork);
    double* work;
    cudaMalloc((void**)&work, lwork * sizeof(double));
    // Perform Cholesky decomposition
    cusolverDnDpotrf(cusolver_handle, cublas_fill_mode(uplo), n, A, n, work, lwork, nullptr);
    cudaFree(work);
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    cusolverDnCpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex*>(A), n, &lwork);
    cuComplex* work;
    cudaMalloc((void**)&work, lwork * sizeof(cuComplex));
    // Perform Cholesky decomposition
    cusolverDnCpotrf(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuComplex*>(A), n, work, lwork, nullptr);
    cudaFree(work);
}
static inline
void potrf (cusolverDnHandle_t& cusolver_handle, const char& uplo, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    cusolverDnZpotrf_bufferSize(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex*>(A), n, &lwork);
    cuDoubleComplex* work;
    cudaMalloc((void**)&work, lwork * sizeof(cuDoubleComplex));
    // Perform Cholesky decomposition
    cusolverDnZpotrf(cusolver_handle, cublas_fill_mode(uplo), n, reinterpret_cast<cuDoubleComplex*>(A), n, work, lwork, nullptr);
    cudaFree(work);
}


static inline
void dnevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    float* d_work = nullptr;
    cudaMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverDnSsyevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, A, lda, W, &lwork);
    // allocate memery
    cudaMalloc((void**)&d_work, sizeof(float) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverDnSsyevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, d_work, lwork, d_info);

    cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    cudaFree(d_info);
    cudaFree(d_work);
}
static inline
void dnevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*    d_info = nullptr;
    double* d_work = nullptr;
    cudaMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverDnDsyevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, A, lda, W, &lwork);
    // allocate memery
    cudaMalloc((void**)&d_work, sizeof(double) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverDnDsyevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, A, lda, W, d_work, lwork, d_info);

    cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    cudaFree(d_info);
    cudaFree(d_work);
}
static inline
void dnevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, float * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*    d_info = nullptr;
    cuComplex* d_work = nullptr;
    cudaMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverDnCheevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, reinterpret_cast<cuComplex*>(A), lda, W, &lwork);
    // allocate memery
    cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverDnCheevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuComplex*>(A), lda, W, d_work, lwork, d_info);

    cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    cudaFree(d_info);
    cudaFree(d_work);
}
static inline
void dnevd (cusolverDnHandle_t& cusolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, double* W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*    d_info = nullptr;
    cuDoubleComplex* d_work = nullptr;
    cudaMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverDnZheevd_bufferSize(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, W, &lwork);
    // allocate memery
    cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverDnZheevd(cusolver_handle, cublas_eig_mode(jobz), cublas_fill_mode(uplo),
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, W, d_work, lwork, d_info);

    cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    cudaFree(d_info);
    cudaFree(d_work);
}

static inline
void dngvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float* B, const int& ldb, float * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    float* d_work = nullptr;
    cudaMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverDnSsygvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, &lwork);
    // allocate memery
    cudaMalloc((void**)&d_work, sizeof(float) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverDnSsygvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, d_work, lwork, d_info);

    cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    cudaFree(d_info);
    cudaFree(d_work);
}
static inline
void dngvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double* B, const int& ldb, double * W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    double* d_work = nullptr;
    cudaMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverDnDsygvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, &lwork);
    // allocate memery
    cudaMalloc((void**)&d_work, sizeof(double) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverDnDsygvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, d_work, lwork, d_info);

    cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    cudaFree(d_info);
    cudaFree(d_work);
}
static inline
void dngvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, float* W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    cuComplex* d_work = nullptr;
    cudaMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverDnChegvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, reinterpret_cast<cuComplex*>(A), lda, reinterpret_cast<cuComplex*>(B), ldb, W, &lwork);
    // allocate memery
    cudaMalloc((void**)&d_work, sizeof(cuComplex) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverDnChegvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, reinterpret_cast<cuComplex*>(A), lda, reinterpret_cast<cuComplex*>(B), ldb, W, d_work, lwork, d_info);

    cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    cudaFree(d_info);
    cudaFree(d_work);
}
static inline
void dngvd (cusolverDnHandle_t& cusolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, double* W)
{
    // prepare some values for cusolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    cuDoubleComplex* d_work = nullptr;
    cudaMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    cusolverDnZhegvd_bufferSize(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, reinterpret_cast<cuDoubleComplex*>(B), ldb, W, &lwork);
    // allocate memery
    cudaMalloc((void**)&d_work, sizeof(cuDoubleComplex) * lwork);
    // compute eigenvalues and eigenvectors.
    cusolverDnZhegvd(cusolver_handle, cublas_eig_type(itype), cublas_eig_mode(jobz), cublas_fill_mode(uplo), 
                                n, reinterpret_cast<cuDoubleComplex*>(A), lda, reinterpret_cast<cuDoubleComplex*>(B), ldb, W, d_work, lwork, d_info);

    cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    cudaFree(d_info);
    cudaFree(d_work);
}

} // namespace cuSolverConnector
} // namespace container

#endif // BASE_THIRD_PARTY_CUSOLVER_H_