#ifndef BASE_THIRD_PARTY_LAPACK_H_
#define BASE_THIRD_PARTY_LAPACK_H_

#include <complex>

#if __CUDA || __ROCM
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <base/macros/cuda.h>
#endif // __CUDA || __ROCM

//Naming convention of lapack subroutines : ammxxx, where
//"a" specifies the data type:
//  - d stands for double
//  - z stands for complex double
//"mm" specifies the type of matrix, for example:
//  - he stands for hermitian
//  - sy stands for symmetric
//"xxx" specifies the type of problem, for example:
//  - gv stands for generalized eigenvalue

extern "C"
{
int ilaenv_(int* ispec,const char* name,const char* opts,
            const int* n1,const int* n2,const int* n3,const int* n4);


// solve the generalized eigenproblem Ax=eBx, where A is Hermitian and complex couble
// zhegv_ & zhegvd_ returns all eigenvalues while zhegvx_ returns selected ones
void ssygvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             float* a, const int* lda,
             const float* b, const int* ldb, float* w,
             float* work, int* lwork,
             int* iwork, int* liwork, int* info);

void dsygvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             double* a, const int* lda,
             const double* b, const int* ldb, double* w,
             double* work, int* lwork,
             int* iwork, int* liwork, int* info);

void chegvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             std::complex<float>* a, const int* lda,
             const std::complex<float>* b, const int* ldb, float* w,
             std::complex<float>* work, int* lwork, float* rwork, int* lrwork,
             int* iwork, int* liwork, int* info);

void zhegvd_(const int* itype, const char* jobz, const char* uplo, const int* n,
             std::complex<double>* a, const int* lda,
             const std::complex<double>* b, const int* ldb, double* w,
             std::complex<double>* work, int* lwork, double* rwork, int* lrwork,
             int* iwork, int* liwork, int* info);

void ssyevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             float *a, const int* lda,
             const float* vl, const float* vu, const int* il, const int* iu, const float* abstol,
             const int* m, float* w, float *z, const int *ldz,
             float *work, const int* lwork, float* rwork, int* iwork, int* ifail, int* info);
void dsyevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             double *a, const int* lda,
             const double* vl, const double* vu, const int* il, const int* iu, const double* abstol,
             const int* m, double* w, double *z, const int *ldz,
             double *work, const int* lwork, double* rwork, int* iwork, int* ifail, int* info);
void cheevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             std::complex<float> *a, const int* lda,
             const float* vl, const float* vu, const int* il, const int* iu, const float* abstol,
             const int* m, float* w, std::complex<float> *z, const int *ldz,
             std::complex<float> *work, const int* lwork, float* rwork, int* iwork, int* ifail, int* info);
void zheevx_(const char* jobz, const char* range, const char* uplo, const int* n,
             std::complex<double> *a, const int* lda,
             const double* vl, const double* vu, const int* il, const int* iu, const double* abstol,
             const int* m, double* w, std::complex<double> *z, const int *ldz,
             std::complex<double> *work, const int* lwork, double* rwork, int* iwork, int* ifail, int* info);

void ssyevd_(const char *jobz, const char *uplo, const int *n,
         float *a, const int *lda, float *w,
         float *work, int *lwork,
         int *iwork, int *liwork, int *info);
void dsyevd_(const char *jobz, const char *uplo, const int *n,
             double *a, const int *lda, double *w,
             double *work, int *lwork,
             int *iwork, int *liwork, int *info);
void cheevd_(const char *jobz, const char *uplo, const int *n,
         std::complex<float> *a, const int *lda, float *w,
         std::complex<float> *work, int *lwork, float *rwork, int *lrwork,
         int *iwork, int *liwork, int *info);
void zheevd_(const char *jobz, const char *uplo, const int *n,
             std::complex<double> *a, const int *lda, double *w,
             std::complex<double> *work, int *lwork, double *rwork, int *lrwork,
             int *iwork, int *liwork, int *info);

void spotrf_(const char*const uplo, const int*const n, float*const A, const int*const lda, int*const info);
void dpotrf_(const char*const uplo, const int*const n, double*const A, const int*const lda, int*const info);
void cpotrf_(const char*const uplo, const int*const n, std::complex<float>*const A, const int*const lda, int*const info);
void zpotrf_(const char*const uplo, const int*const n, std::complex<double>*const A, const int*const lda, int*const info);

void spotri_(const char*const uplo, const int*const n, float*const A, const int*const lda, int*const info);
void dpotri_(const char*const uplo, const int*const n, double*const A, const int*const lda, int*const info);
void cpotri_(const char*const uplo, const int*const n, std::complex<float>*const A, const int*const lda, int*const info);
void zpotri_(const char*const uplo, const int*const n, std::complex<double>*const A, const int*const lda, int*const info);

void strtri_(const char* uplo, const char* diag, const int* n, float* a, const int* lda, int* info);
void dtrtri_(const char* uplo, const char* diag, const int* n, double* a, const int* lda, int* info);
void ctrtri_(const char* uplo, const char* diag, const int* n, std::complex<float>* a, const int* lda, int* info);
void ztrtri_(const char* uplo, const char* diag, const int* n, std::complex<double>* a, const int* lda, int* info);

}

// Class LapackConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	LapackConnector::functionname(parameter list).
namespace container {
namespace lapackConnector
{
static inline
int ilaenv( int ispec, const char *name,const char *opts,const int n1,const int n2,
            const int n3,const int n4)
{
    const int nb = ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4);
    return nb;
}
// wrap function of fortran lapack routine zhegvd. (pointer version)
static inline
void dngvd(const int itype, const char jobz, const char uplo, const int n,
            float* a, const int lda,
            const float* b, const int ldb, float* w,
            float* work, int lwork, float* rwork, int lrwork,
            int* iwork, int liwork, int info)
{
    // call the fortran routine
    ssygvd_(&itype, &jobz, &uplo, &n,
            a, &lda, b, &ldb, w,
            work, &lwork,
            iwork, &liwork, &info);
}
// wrap function of fortran lapack routine zhegvd.
static inline
void dngvd(const int itype, const char jobz, const char uplo, const int n,
            double* a, const int lda,
            const double* b, const int ldb, double* w,
            double* work, int lwork, double* rwork, int lrwork,
            int* iwork, int liwork, int info)
{
    // call the fortran routine
    dsygvd_(&itype, &jobz, &uplo, &n,
            a, &lda, b, &ldb, w,
            work, &lwork,
            iwork, &liwork, &info);
}
static inline
void dngvd(const int itype, const char jobz, const char uplo, const int n,
            std::complex<float>* a, const int lda,
            const std::complex<float>* b, const int ldb, float* w,
            std::complex<float>* work, int lwork, float* rwork, int lrwork,
            int* iwork, int liwork, int info)
{
    // call the fortran routine
    chegvd_(&itype, &jobz, &uplo, &n,
            a, &lda, b, &ldb, w,
            work, &lwork, rwork, &lrwork,
            iwork, &liwork, &info);
}
// wrap function of fortran lapack routine zhegvd.
static inline
void dngvd(const int itype, const char jobz, const char uplo, const int n,
            std::complex<double>* a, const int lda,
            const std::complex<double>* b, const int ldb, double* w,
            std::complex<double>* work, int lwork, double* rwork, int lrwork,
            int* iwork, int liwork, int info)
{
    // call the fortran routine
    zhegvd_(&itype, &jobz, &uplo, &n,
            a, &lda, b, &ldb, w,
            work, &lwork, rwork, &lrwork,
            iwork, &liwork, &info);
}

// wrap function of fortran lapack routine zheevx.
static inline
void dnevx( const int itype, const char jobz, const char range, const char uplo, const int n,
             float* a, const int lda,
             const float vl, const float vu, const int il, const int iu, const float abstol,
             const int m, float* w, float* z, const int ldz,
             float* work, const int lwork, float* rwork, int* iwork, int* ifail, int info)
{
    ssyevx_(&jobz, &range, &uplo, &n,
            a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
}
// wrap function of fortran lapack routine zheevx.
static inline
void dnevx( const int itype, const char jobz, const char range, const char uplo, const int n,
             double* a, const int lda,
             const double vl, const double vu, const int il, const int iu, const double abstol,
             const int m, double* w, double* z, const int ldz,
             double* work, const int lwork, double* rwork, int* iwork, int* ifail, int info)
{
    dsyevx_(&jobz, &range, &uplo, &n,
            a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
}
static inline
void dnevx( const int itype, const char jobz, const char range, const char uplo, const int n,
             std::complex<float>* a, const int lda,
             const float vl, const float vu, const int il, const int iu, const float abstol,
             const int m, float* w, std::complex<float>* z, const int ldz,
             std::complex<float>* work, const int lwork, float* rwork, int* iwork, int* ifail, int info)
{
    cheevx_(&jobz, &range, &uplo, &n,
            a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
}
// wrap function of fortran lapack routine zheevx.
static inline
void dnevx( const int itype, const char jobz, const char range, const char uplo, const int n,
             std::complex<double>* a, const int lda,
             const double vl, const double vu, const int il, const int iu, const double abstol,
             const int m, double* w, std::complex<double>* z, const int ldz,
             std::complex<double>* work, const int lwork, double* rwork, int* iwork, int* ifail, int info)
{
    zheevx_(&jobz, &range, &uplo, &n,
            a, &lda, &vl, &vu, &il, &iu,
            &abstol, &m, w, z, &ldz,
            work, &lwork, rwork, iwork, ifail, &info);
}

static inline
void dnevd(const char jobz, const char uplo, const int n,
            float* a, const int lda, float* w,
            float* work, int lwork, float* rwork, int lrwork,
            int* iwork, int liwork, int& info)
{
    // call the fortran routine
    ssyevd_( &jobz, &uplo, &n,
             a, &lda, w,
             work, &lwork,
             iwork, &liwork, &info);
}
// wrap function of fortran lapack routine zhegvd.
static inline
void dnevd(const char jobz, const char uplo, const int n,
            double* a, const int lda, double* w,
            double* work, int lwork, double* rwork, int lrwork,
            int* iwork, int liwork, int& info)
{
    // call the fortran routine
    dsyevd_( &jobz, &uplo, &n,
             a, &lda, w,
             work, &lwork,
             iwork, &liwork, &info);
}
static inline
void dnevd(const char jobz, const char uplo, const int n,
            std::complex<float>* a, const int lda, float* w,
            std::complex<float>* work, int lwork, float* rwork, int lrwork,
            int* iwork, int liwork, int& info)
{
    // call the fortran routine
    cheevd_( &jobz, &uplo, &n,
             a, &lda, w,
             work, &lwork, rwork, &lrwork,
             iwork, &liwork, &info);
}
// wrap function of fortran lapack routine zhegvd.
static inline
void dnevd(const char jobz, const char uplo, const int n,
            std::complex<double>* a, const int lda, double* w,
            std::complex<double>* work, int lwork, double* rwork, int lrwork,
            int* iwork, int liwork, int& info)
{
    // call the fortran routine
    zheevd_( &jobz, &uplo, &n,
             a, &lda, w,
             work, &lwork, rwork, &lrwork,
             iwork, &liwork, &info);
}

static inline
void potrf( const char &uplo, const int &n, float* A, const int &lda, int &info )
{
	spotrf_(&uplo, &n, A, &lda, &info );
}	
static inline
void potrf( const char &uplo, const int &n, double* A, const int &lda, int &info )
{
	dpotrf_(&uplo, &n, A, &lda, &info );
}	
static inline
void potrf( const char &uplo, const int &n, std::complex<float>* A, const int &lda, int &info )
{
	cpotrf_(&uplo, &n, A, &lda, &info );
}	
static inline
void potrf( const char &uplo, const int &n, std::complex<double>* A, const int &lda, int &info )
{
	zpotrf_( &uplo, &n, A, &lda, &info );
}

static inline
void trtri( const char &uplo, const char &diag, const int &n, float* A, const int &lda, int &info )
{
    strtri_( &uplo, &diag, &n, A, &lda, &info);
}
static inline
void trtri( const char &uplo, const char &diag, const int &n, double* A, const int &lda, int &info)
{
    dtrtri_( &uplo, &diag, &n, A, &lda, &info);
}
static inline
void trtri( const char &uplo, const char &diag, const int &n, std::complex<float>* A, const int &lda, int &info )
{
    ctrtri_( &uplo, &diag, &n, A, &lda, &info);
}
static inline
void trtri( const char &uplo, const char &diag, const int &n, std::complex<double>* A, const int &lda, int &info)
{
    ztrtri_( &uplo, &diag, &n, A, &lda, &info);
}

} // namespace lapackConnector

#if __CUDA || __ROCM
namespace cuSolverConnector {

template <typename T>
static inline
void trtri (cusolverDnHandle_t& cusolver_handle, const char& uplo, const char& diag, const int& n, T* A, const int& lda)
{
    size_t d_lwork = 0, h_lwork = 0;
    using Type = typename PossibleStdComplexToThrustComplex<T>::type;
    cusolverDnXtrtri_bufferSize(cusolver_handle, cublas_fill_mode(uplo), cublas_diag_type(diag), n, DataTypeToCudaType<T>::cuda_data_type, reinterpret_cast<Type*>(A), lda, &d_lwork, &h_lwork);
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
    cusolverDnXtrtri(cusolver_handle, cublas_fill_mode(uplo), cublas_diag_type(diag), n, DataTypeToCudaType<T>::cuda_data_type, reinterpret_cast<Type*>(A), n, d_work, d_lwork, h_work, h_lwork, d_info);
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
#endif

}
#endif  // BASE_THIRD_PARTY_LAPACK_H_
