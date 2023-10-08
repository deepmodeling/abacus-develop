#ifndef BASE_THIRD_PARTY_HIPSOLVER_H_
#define BASE_THIRD_PARTY_HIPSOLVER_H_

#include <base/macros/rocm.h>

#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <hipsolver/hipsolver.h>

namespace container {
namespace hipSolverConnector {

template <typename T>
static inline
void trtri (hipsolverHandle_t& hipsolver_handle, const char& uplo, const char& diag, const int& n, T* A, const int& lda)
{
    size_t d_lwork = 0, h_lwork = 0;
    hipsolverDnXtrtri_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), hipblas_diag_type(diag), n, GetTypeRocm<T>::cuda_data_type, A, lda, &d_lwork, &h_lwork);
    void* d_work = nullptr, *h_work = nullptr;
    hipMalloc((void**)&d_work, d_lwork);
    if (h_lwork) {
        h_work = malloc(h_lwork);
        if (h_work == nullptr) {
            throw std::bad_alloc();
        }
    }
    int h_info = 0;
    int* d_info = nullptr;
    hipMalloc((void**)&d_info, sizeof(int));
    // Perform Cholesky decomposition
    hipsolverDnXtrtri(hipsolver_handle, hipsolver_fill_mode(uplo), hipblas_diag_type(diag), n, GetTypeRocm<T>::cuda_data_type, A, n, d_work, d_lwork, h_work, h_lwork, d_info);
    hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("trtri: failed to invert matrix");
    }
    free(h_work);
    hipFree(d_work);
    hipFree(d_info);
}

static inline
void potri (hipsolverHandle_t& hipsolver_handle, const char& uplo, const char& diag, const int& n, float * A, const int& lda)
{
    int lwork;
    hipsolverDnSpotri_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, &lwork);
    float* work;
    hipMalloc((void**)&work, lwork * sizeof(float));
    // Perform Cholesky decomposition
    hipsolverDnSpotri(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, work, lwork, nullptr);
    hipFree(work);
}
static inline
void potri (hipsolverHandle_t& hipsolver_handle, const char& uplo, const char& diag, const int& n, double * A, const int& lda)
{
    int lwork;
    hipsolverDnDpotri_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, &lwork);
    double* work;
    hipMalloc((void**)&work, lwork * sizeof(double));
    // Perform Cholesky decomposition
    hipsolverDnDpotri(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, work, lwork, nullptr);
    hipFree(work);
}
static inline
void potri (hipsolverHandle_t& hipsolver_handle, const char& uplo, const char& diag, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    hipsolverDnCpotri_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipFloatComplex *>(A), n, &lwork);
    hipFloatComplex* work;
    hipMalloc((void**)&work, lwork * sizeof(hipFloatComplex));
    // Perform Cholesky decomposition
    hipsolverDnCpotri(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipFloatComplex *>(A), n, work, lwork, nullptr);
    hipFree(work);
}
static inline
void potri (hipsolverHandle_t& hipsolver_handle, const char& uplo, const char& diag, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    hipsolverDnZpotri_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipDoubleComplex *>(A), n, &lwork);
    hipDoubleComplex* work;
    hipMalloc((void**)&work, lwork * sizeof(hipDoubleComplex));
    // Perform Cholesky decomposition
    hipsolverDnZpotri(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipDoubleComplex *>(A), n, work, lwork, nullptr);
    hipFree(work);
}


static inline
void potrf (hipsolverHandle_t& hipsolver_handle, const char& uplo, const int& n, float * A, const int& lda)
{
    int lwork;
    hipsolverDnSpotrf_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, &lwork);
    float* work;
    hipMalloc((void**)&work, lwork * sizeof(float));
    // Perform Cholesky decomposition
    hipsolverDnSpotrf(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, work, lwork, nullptr);
    hipFree(work);
}
static inline
void potrf (hipsolverHandle_t& hipsolver_handle, const char& uplo, const int& n, double * A, const int& lda)
{
    int lwork;
    hipsolverDnDpotrf_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, &lwork);
    double* work;
    hipMalloc((void**)&work, lwork * sizeof(double));
    // Perform Cholesky decomposition
    hipsolverDnDpotrf(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, work, lwork, nullptr);
    hipFree(work);
}
static inline
void potrf (hipsolverHandle_t& hipsolver_handle, const char& uplo, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    hipsolverDnCpotrf_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipFloatComplex*>(A), n, &lwork);
    hipFloatComplex* work;
    hipMalloc((void**)&work, lwork * sizeof(hipFloatComplex));
    // Perform Cholesky decomposition
    hipsolverDnCpotrf(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipFloatComplex*>(A), n, work, lwork, nullptr);
    hipFree(work);
}
static inline
void potrf (hipsolverHandle_t& hipsolver_handle, const char& uplo, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    hipsolverDnZpotrf_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipDoubleComplex*>(A), n, &lwork);
    hipDoubleComplex* work;
    hipMalloc((void**)&work, lwork * sizeof(hipDoubleComplex));
    // Perform Cholesky decomposition
    hipsolverDnZpotrf(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipDoubleComplex*>(A), n, work, lwork, nullptr);
    hipFree(work);
}


static inline
void dnevd (hipsolverHandle_t& hipsolver_handle, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float * W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    float* d_work = nullptr;
    hipMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverDnSsyevd_bufferSize(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, W, &lwork);
    // allocate memery
    hipMalloc((void**)&d_work, sizeof(float) * lwork);
    // compute eigenvalues and eigenvectors.
    hipsolverDnSsyevd(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo),
                                n, A, lda, W, d_work, lwork, d_info);

    hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipFree(d_info);
    hipFree(d_work);
}
static inline
void dnevd (hipsolverHandle_t& hipsolver_handle, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double * W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*    d_info = nullptr;
    double* d_work = nullptr;
    hipMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverDnDsyevd_bufferSize(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, W, &lwork);
    // allocate memery
    hipMalloc((void**)&d_work, sizeof(double) * lwork);
    // compute eigenvalues and eigenvectors.
    hipsolverDnDsyevd(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo),
                                n, A, lda, W, d_work, lwork, d_info);

    hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipFree(d_info);
    hipFree(d_work);
}
static inline
void dnevd (hipsolverHandle_t& hipsolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, float * W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*    d_info = nullptr;
    hipFloatComplex* d_work = nullptr;
    hipMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverDnCheevd_bufferSize(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipFloatComplex*>(A), lda, W, &lwork);
    // allocate memery
    hipMalloc((void**)&d_work, sizeof(hipFloatComplex) * lwork);
    // compute eigenvalues and eigenvectors.
    hipsolverDnCheevd(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo),
                                n, reinterpret_cast<hipFloatComplex*>(A), lda, W, d_work, lwork, d_info);

    hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipFree(d_info);
    hipFree(d_work);
}
static inline
void dnevd (hipsolverHandle_t& hipsolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, double* W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*    d_info = nullptr;
    hipDoubleComplex* d_work = nullptr;
    hipMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverDnZheevd_bufferSize(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipDoubleComplex*>(A), lda, W, &lwork);
    // allocate memery
    hipMalloc((void**)&d_work, sizeof(hipDoubleComplex) * lwork);
    // compute eigenvalues and eigenvectors.
    hipsolverDnZheevd(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo),
                                n, reinterpret_cast<hipDoubleComplex*>(A), lda, W, d_work, lwork, d_info);

    hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipFree(d_info);
    hipFree(d_work);
}

static inline
void dngvd (hipsolverHandle_t& hipsolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float* B, const int& ldb, float * W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    float* d_work = nullptr;
    hipMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverDnSsygvd_bufferSize(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, &lwork);
    // allocate memery
    hipMalloc((void**)&d_work, sizeof(float) * lwork);
    // compute eigenvalues and eigenvectors.
    hipsolverDnSsygvd(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, d_work, lwork, d_info);

    hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipFree(d_info);
    hipFree(d_work);
}
static inline
void dngvd (hipsolverHandle_t& hipsolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double* B, const int& ldb, double * W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    double* d_work = nullptr;
    hipMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverDnDsygvd_bufferSize(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, &lwork);
    // allocate memery
    hipMalloc((void**)&d_work, sizeof(double) * lwork);
    // compute eigenvalues and eigenvectors.
    hipsolverDnDsygvd(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, d_work, lwork, d_info);

    hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipFree(d_info);
    hipFree(d_work);
}
static inline
void dngvd (hipsolverHandle_t& hipsolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, float* W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    hipFloatComplex* d_work = nullptr;
    hipMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverDnChegvd_bufferSize(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipFloatComplex*>(A), lda, reinterpret_cast<hipFloatComplex*>(B), ldb, W, &lwork);
    // allocate memery
    hipMalloc((void**)&d_work, sizeof(hipFloatComplex) * lwork);
    // compute eigenvalues and eigenvectors.
    hipsolverDnChegvd(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipFloatComplex*>(A), lda, reinterpret_cast<hipFloatComplex*>(B), ldb, W, d_work, lwork, d_info);

    hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipFree(d_info);
    hipFree(d_work);
}
static inline
void dngvd (hipsolverHandle_t& hipsolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, double* W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    hipDoubleComplex* d_work = nullptr;
    hipMalloc((void**)&d_info, sizeof(int));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverDnZhegvd_bufferSize(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipDoubleComplex*>(A), lda, reinterpret_cast<hipDoubleComplex*>(B), ldb, W, &lwork);
    // allocate memery
    hipMalloc((void**)&d_work, sizeof(hipDoubleComplex) * lwork);
    // compute eigenvalues and eigenvectors.
    hipsolverDnZhegvd(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipDoubleComplex*>(A), lda, reinterpret_cast<hipDoubleComplex*>(B), ldb, W, d_work, lwork, d_info);

    hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost);
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipFree(d_info);
    hipFree(d_work);
}

} // namespace hipSolverConnector
} // namespace container

#endif // BASE_THIRD_PARTY_HIPSOLVER_H_