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
    hipsolverErrcheck(hipsolverDnXtrtri_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), hipblas_diag_type(diag), n, GetTypeRocm<T>::cuda_data_type, A, lda, &d_lwork, &h_lwork));
    void* d_work = nullptr, *h_work = nullptr;
    hipErrcheck(hipMalloc((void**)&d_work, d_lwork));
    if (h_lwork) {
        h_work = malloc(h_lwork);
        if (h_work == nullptr) {
            throw std::bad_alloc();
        }
    }
    int h_info = 0;
    int* d_info = nullptr;
    hipErrcheck(hipMalloc((void**)&d_info, sizeof(int)));
    // Perform Cholesky decomposition
    hipsolverErrcheck(hipsolverDnXtrtri(hipsolver_handle, hipsolver_fill_mode(uplo), hipblas_diag_type(diag), n, GetTypeRocm<T>::cuda_data_type, A, n, d_work, d_lwork, h_work, h_lwork, d_info));
    hipErrcheck(hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("trtri: failed to invert matrix");
    }
    free(h_work);
    hipErrcheck(hipFree(d_work));
    hipErrcheck(hipFree(d_info));
}

static inline
void potri (hipsolverHandle_t& hipsolver_handle, const char& uplo, const char& diag, const int& n, float * A, const int& lda)
{
    int lwork;
    hipsolverErrcheck(hipsolverDnSpotri_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, &lwork));
    float* work;
    hipErrcheck(hipMalloc((void**)&work, lwork * sizeof(float)));
    // Perform Cholesky decomposition
    hipsolverErrcheck(hipsolverDnSpotri(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, work, lwork, nullptr));
    hipErrcheck(hipFree(work));
}
static inline
void potri (hipsolverHandle_t& hipsolver_handle, const char& uplo, const char& diag, const int& n, double * A, const int& lda)
{
    int lwork;
    hipsolverErrcheck(hipsolverDnDpotri_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, &lwork));
    double* work;
    hipErrcheck(hipMalloc((void**)&work, lwork * sizeof(double)));
    // Perform Cholesky decomposition
    hipsolverErrcheck(hipsolverDnDpotri(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, work, lwork, nullptr));
    hipErrcheck(hipFree(work));
}
static inline
void potri (hipsolverHandle_t& hipsolver_handle, const char& uplo, const char& diag, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    hipsolverErrcheck(hipsolverDnCpotri_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipFloatComplex *>(A), n, &lwork));
    hipFloatComplex* work;
    hipErrcheck(hipMalloc((void**)&work, lwork * sizeof(hipFloatComplex)));
    // Perform Cholesky decomposition
    hipsolverErrcheck(hipsolverDnCpotri(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipFloatComplex *>(A), n, work, lwork, nullptr));
    hipErrcheck(hipFree(work));
}
static inline
void potri (hipsolverHandle_t& hipsolver_handle, const char& uplo, const char& diag, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    hipsolverErrcheck(hipsolverDnZpotri_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipDoubleComplex *>(A), n, &lwork));
    hipDoubleComplex* work;
    hipErrcheck(hipMalloc((void**)&work, lwork * sizeof(hipDoubleComplex)));
    // Perform Cholesky decomposition
    hipsolverErrcheck(hipsolverDnZpotri(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipDoubleComplex *>(A), n, work, lwork, nullptr));
    hipErrcheck(hipFree(work));
}


static inline
void potrf (hipsolverHandle_t& hipsolver_handle, const char& uplo, const int& n, float * A, const int& lda)
{
    int lwork;
    hipsolverErrcheck(hipsolverDnSpotrf_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, &lwork));
    float* work;
    hipErrcheck(hipMalloc((void**)&work, lwork * sizeof(float)));
    // Perform Cholesky decomposition
    hipsolverErrcheck(hipsolverDnSpotrf(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, work, lwork, nullptr));
    hipErrcheck(hipFree(work));
}
static inline
void potrf (hipsolverHandle_t& hipsolver_handle, const char& uplo, const int& n, double * A, const int& lda)
{
    int lwork;
    hipsolverErrcheck(hipsolverDnDpotrf_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, &lwork));
    double* work;
    hipErrcheck(hipMalloc((void**)&work, lwork * sizeof(double)));
    // Perform Cholesky decomposition
    hipsolverErrcheck(hipsolverDnDpotrf(hipsolver_handle, hipsolver_fill_mode(uplo), n, A, n, work, lwork, nullptr));
    hipErrcheck(hipFree(work));
}
static inline
void potrf (hipsolverHandle_t& hipsolver_handle, const char& uplo, const int& n, std::complex<float> * A, const int& lda)
{
    int lwork;
    hipsolverErrcheck(hipsolverDnCpotrf_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipFloatComplex*>(A), n, &lwork));
    hipFloatComplex* work;
    hipErrcheck(hipMalloc((void**)&work, lwork * sizeof(hipFloatComplex)));
    // Perform Cholesky decomposition
    hipsolverErrcheck(hipsolverDnCpotrf(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipFloatComplex*>(A), n, work, lwork, nullptr));
    hipErrcheck(hipFree(work));
}
static inline
void potrf (hipsolverHandle_t& hipsolver_handle, const char& uplo, const int& n, std::complex<double> * A, const int& lda)
{
    int lwork;
    hipsolverErrcheck(hipsolverDnZpotrf_bufferSize(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipDoubleComplex*>(A), n, &lwork));
    hipDoubleComplex* work;
    hipErrcheck(hipMalloc((void**)&work, lwork * sizeof(hipDoubleComplex)));
    // Perform Cholesky decomposition
    hipsolverErrcheck(hipsolverDnZpotrf(hipsolver_handle, hipsolver_fill_mode(uplo), n, reinterpret_cast<hipDoubleComplex*>(A), n, work, lwork, nullptr));
    hipErrcheck(hipFree(work));
}


static inline
void dnevd (hipsolverHandle_t& hipsolver_handle, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float * W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    float* d_work = nullptr;
    hipErrcheck(hipMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverErrcheck(hipsolverDnSsyevd_bufferSize(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, W, &lwork));
    // allocate memery
    hipErrcheck(hipMalloc((void**)&d_work, sizeof(float) * lwork));
    // compute eigenvalues and eigenvectors.
    hipsolverErrcheck(hipsolverDnSsyevd(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo),
                                n, A, lda, W, d_work, lwork, d_info));

    hipErrcheck(hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipErrcheck(hipFree(d_info));
    hipErrcheck(hipFree(d_work));
}
static inline
void dnevd (hipsolverHandle_t& hipsolver_handle, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double * W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*    d_info = nullptr;
    double* d_work = nullptr;
    hipErrcheck(hipMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverErrcheck(hipsolverDnDsyevd_bufferSize(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, W, &lwork));
    // allocate memery
    hipErrcheck(hipMalloc((void**)&d_work, sizeof(double) * lwork));
    // compute eigenvalues and eigenvectors.
    hipsolverErrcheck(hipsolverDnDsyevd(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo),
                                n, A, lda, W, d_work, lwork, d_info));

    hipErrcheck(hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipErrcheck(hipFree(d_info));
    hipErrcheck(hipFree(d_work));
}
static inline
void dnevd (hipsolverHandle_t& hipsolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, float * W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*    d_info = nullptr;
    hipFloatComplex* d_work = nullptr;
    hipErrcheck(hipMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverErrcheck(hipsolverDnCheevd_bufferSize(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipFloatComplex*>(A), lda, W, &lwork));
    // allocate memery
    hipErrcheck(hipMalloc((void**)&d_work, sizeof(hipFloatComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    hipsolverErrcheck(hipsolverDnCheevd(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo),
                                n, reinterpret_cast<hipFloatComplex*>(A), lda, W, d_work, lwork, d_info));

    hipErrcheck(hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipErrcheck(hipFree(d_info));
    hipErrcheck(hipFree(d_work));
}
static inline
void dnevd (hipsolverHandle_t& hipsolver_handle, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, double* W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*    d_info = nullptr;
    hipDoubleComplex* d_work = nullptr;
    hipErrcheck(hipMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverErrcheck(hipsolverDnZheevd_bufferSize(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipDoubleComplex*>(A), lda, W, &lwork));
    // allocate memery
    hipErrcheck(hipMalloc((void**)&d_work, sizeof(hipDoubleComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    hipsolverErrcheck(hipsolverDnZheevd(hipsolver_handle, hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo),
                                n, reinterpret_cast<hipDoubleComplex*>(A), lda, W, d_work, lwork, d_info));

    hipErrcheck(hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipErrcheck(hipFree(d_info));
    hipErrcheck(hipFree(d_work));
}

static inline
void dngvd (hipsolverHandle_t& hipsolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, float* A, const int& lda, float* B, const int& ldb, float * W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    float* d_work = nullptr;
    hipErrcheck(hipMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverErrcheck(hipsolverDnSsygvd_bufferSize(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, &lwork));
    // allocate memery
    hipErrcheck(hipMalloc((void**)&d_work, sizeof(float) * lwork));
    // compute eigenvalues and eigenvectors.
    hipsolverErrcheck(hipsolverDnSsygvd(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, d_work, lwork, d_info));

    hipErrcheck(hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipErrcheck(hipFree(d_info));
    hipErrcheck(hipFree(d_work));
}
static inline
void dngvd (hipsolverHandle_t& hipsolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, double* A, const int& lda, double* B, const int& ldb, double * W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    double* d_work = nullptr;
    hipErrcheck(hipMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverErrcheck(hipsolverDnDsygvd_bufferSize(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, &lwork));
    // allocate memery
    hipErrcheck(hipMalloc((void**)&d_work, sizeof(double) * lwork));
    // compute eigenvalues and eigenvectors.
    hipsolverErrcheck(hipsolverDnDsygvd(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, A, lda, B, ldb, W, d_work, lwork, d_info));

    hipErrcheck(hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipErrcheck(hipFree(d_info));
    hipErrcheck(hipFree(d_work));
}
static inline
void dngvd (hipsolverHandle_t& hipsolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<float>* A, const int& lda, std::complex<float>* B, const int& ldb, float* W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    hipFloatComplex* d_work = nullptr;
    hipErrcheck(hipMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverErrcheck(hipsolverDnChegvd_bufferSize(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipFloatComplex*>(A), lda, reinterpret_cast<hipFloatComplex*>(B), ldb, W, &lwork));
    // allocate memery
    hipErrcheck(hipMalloc((void**)&d_work, sizeof(hipFloatComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    hipsolverErrcheck(hipsolverDnChegvd(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipFloatComplex*>(A), lda, reinterpret_cast<hipFloatComplex*>(B), ldb, W, d_work, lwork, d_info));

    hipErrcheck(hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipErrcheck(hipFree(d_info));
    hipErrcheck(hipFree(d_work));
}
static inline
void dngvd (hipsolverHandle_t& hipsolver_handle, const int& itype, const char& jobz, const char& uplo, const int& n, std::complex<double>* A, const int& lda, std::complex<double>* B, const int& ldb, double* W)
{
    // prepare some values for hipsolverDnZhegvd_bufferSize
    int lwork  = 0; 
    int h_info = 0; 
    int*   d_info = nullptr;
    hipDoubleComplex* d_work = nullptr;
    hipErrcheck(hipMalloc((void**)&d_info, sizeof(int)));

    // calculate the sizes needed for pre-allocated buffer.
    hipsolverErrcheck(hipsolverDnZhegvd_bufferSize(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipDoubleComplex*>(A), lda, reinterpret_cast<hipDoubleComplex*>(B), ldb, W, &lwork));
    // allocate memery
    hipErrcheck(hipMalloc((void**)&d_work, sizeof(hipDoubleComplex) * lwork));
    // compute eigenvalues and eigenvectors.
    hipsolverErrcheck(hipsolverDnZhegvd(hipsolver_handle, hipblas_eig_type(itype), hipblas_eig_mode(jobz), hipsolver_fill_mode(uplo), 
                                n, reinterpret_cast<hipDoubleComplex*>(A), lda, reinterpret_cast<hipDoubleComplex*>(B), ldb, W, d_work, lwork, d_info));

    hipErrcheck(hipMemcpy(&h_info, d_info, sizeof(int), hipMemcpyDeviceToHost));
    if (h_info != 0) {
        throw std::runtime_error("dnevd: failed to invert matrix");
    }
    hipErrcheck(hipFree(d_info));
    hipErrcheck(hipFree(d_work));
}

} // namespace hipSolverConnector
} // namespace container

#endif // BASE_THIRD_PARTY_HIPSOLVER_H_