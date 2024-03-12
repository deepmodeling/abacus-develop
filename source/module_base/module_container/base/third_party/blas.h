#ifndef BASE_THIRD_PARTY_BLAS_H_
#define BASE_THIRD_PARTY_BLAS_H_

#include <complex>

#if defined(__CUDA)
#include <base/third_party/cublas.h>
#elif defined(__ROCM)
#include <base/third_party/hipblas.h>
#endif

extern "C"
{
// level 1: std::vector-std::vector operations, O(n) data and O(n) work.

// Peize Lin add ?scal 2016-08-04, to compute x=a*x
void sscal_(const int *N, const float *alpha, float *x, const int *incx);
void dscal_(const int *N, const double *alpha, double *x, const int *incx);
void cscal_(const int *N, const std::complex<float> *alpha, std::complex<float> *x, const int *incx);
void zscal_(const int *N, const std::complex<double> *alpha, std::complex<double> *x, const int *incx);

// Peize Lin add ?axpy 2016-08-04, to compute y=a*x+y
void saxpy_(const int *N, const float *alpha, const float *x, const int *incx, float *y, const int *incy);
void daxpy_(const int *N, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
void caxpy_(const int *N, const std::complex<float> *alpha, const std::complex<float> *x, const int *incx, std::complex<float> *y, const int *incy);
void zaxpy_(const int *N, const std::complex<double> *alpha, const std::complex<double> *x, const int *incx, std::complex<double> *y, const int *incy);

void dcopy_(long const *n, const double *a, int const *incx, double *b, int const *incy);
void zcopy_(long const *n, const std::complex<double> *a, int const *incx, std::complex<double> *b, int const *incy);

//reason for passing results as argument instead of returning it:
//see https://www.numbercrunch.de/blog/2014/07/lost-in-translation/
void cdotc_(const int *n, const std::complex<float> *zx, const int *incx, 
            const std::complex<float> *zy, const int *incy, std::complex<float> *result);
void zdotc_(const int *n, const std::complex<double> *zx, const int *incx, 
            const std::complex<double> *zy, const int *incy, std::complex<double> *result);
// Peize Lin add ?dot 2017-10-27, to compute d=x*y
float sdot_(const int *N, const float *x, const int *incx, const float *y, const int *incy);
double ddot_(const int *N, const double *x, const int *incx, const double *y, const int *incy);

// Peize Lin add ?nrm2 2018-06-12, to compute out = ||x||_2 = \sqrt{ \sum_i x_i**2 }
float snrm2_( const int *n, const float *x, const int *incx );
double dnrm2_( const int *n, const double *x, const int *incx );
double dznrm2_( const int *n, const std::complex<double> *x, const int *incx );

// level 2: matrix-std::vector operations, O(n^2) data and O(n^2) work.
void sgemv_(const char*const transa, const int*const m, const int*const n,
            const float*const alpha, const float*const a, const int*const lda, const float*const x, const int*const incx,
            const float*const eta, float*const y, const int*const incy);
void dgemv_(const char*const transa, const int*const m, const int*const n,
            const double*const alpha, const double*const a, const int*const lda, const double*const x, const int*const incx,
            const double*const beta, double*const y, const int*const incy);

void cgemv_(const char *trans, const int *m, const int *n, const std::complex<float> *alpha,
            const std::complex<float> *a, const int *lda, const std::complex<float> *x, const int *incx,
            const std::complex<float> *beta, std::complex<float> *y, const int *incy);

void zgemv_(const char *trans, const int *m, const int *n, const std::complex<double> *alpha,
            const std::complex<double> *a, const int *lda, const std::complex<double> *x, const int *incx,
            const std::complex<double> *beta, std::complex<double> *y, const int *incy);

void dsymv_(const char *uplo, const int *n,
            const double *alpha, const double *a, const int *lda,
            const double *x, const int *incx,
            const double *beta, double *y, const int *incy);

// A := alpha x * y.T + A
void dger_(const int* m,
           const int* n,
           const double* alpha,
           const double* x,
           const int* incx,
           const double* y,
           const int* incy,
           double* a,
           const int* lda);
void zgerc_(const int* m,
            const int* n,
            const std::complex<double>* alpha,
            const std::complex<double>* x,
            const int* incx,
            const std::complex<double>* y,
            const int* incy,
            std::complex<double>* a,
            const int* lda);

// level 3: matrix-matrix operations, O(n^2) data and O(n^3) work.

// Peize Lin add ?gemm 2017-10-27, to compute C = a * A.? * B.? + b * C
// A is general
void sgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
            const float *alpha, const float *a, const int *lda, const float *b, const int *ldb,
            const float *beta, float *c, const int *ldc);
void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
            const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
            const double *beta, double *c, const int *ldc);
void cgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
            const std::complex<float> *alpha, const std::complex<float> *a, const int *lda, const std::complex<float> *b, const int *ldb,
            const std::complex<float> *beta, std::complex<float> *c, const int *ldc);
void zgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
            const std::complex<double> *alpha, const std::complex<double> *a, const int *lda, const std::complex<double> *b, const int *ldb,
            const std::complex<double> *beta, std::complex<double> *c, const int *ldc);


//a is symmetric
void dsymm_(const char *side, const char *uplo, const int *m, const int *n,
            const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
            const double *beta, double *c, const int *ldc);
//a is hermitian
void zhemm_(char *side, char *uplo, int *m, int *n,std::complex<double> *alpha,
            std::complex<double> *a,  int *lda,  std::complex<double> *b, int *ldb, std::complex<double> *beta, std::complex<double> *c, int *ldc);

//solving triangular matrix with multiple right hand sides
void dtrsm_(char *side, char* uplo, char *transa, char *diag, int *m, int *n,
            double* alpha, double* a, int *lda, double*b, int *ldb);
void ztrsm_(char *side, char* uplo, char *transa, char *diag, int *m, int *n,
            std::complex<double>* alpha, std::complex<double>* a, int *lda, std::complex<double>*b, int *ldb);

}

namespace container {

// Class BlasConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	BlasConnector::functionname(parameter list).
namespace BlasConnector {

static inline
void axpy( const int& n, const float& alpha, const float *x, const int& incx, float *y, const int& incy)
{
    saxpy_(&n, &alpha, x, &incx, y, &incy);
}
static inline
void axpy( const int& n, const double& alpha, const double *x, const int& incx, double *y, const int& incy)
{
    daxpy_(&n, &alpha, x, &incx, y, &incy);
}
static inline
void axpy( const int& n, const std::complex<float>& alpha, const std::complex<float> *x, const int& incx, std::complex<float> *y, const int& incy)
{
    caxpy_(&n, &alpha, x, &incx, y, &incy);
}
static inline
void axpy( const int& n, const std::complex<double>& alpha, const std::complex<double> *x, const int& incx, std::complex<double> *y, const int& incy)
{
    zaxpy_(&n, &alpha, x, &incx, y, &incy);
}

// Peize Lin add 2016-08-04
// x=a*x
static inline
void scal( const int& n,  const float& alpha, float *x, const int& incx)
{
    sscal_(&n, &alpha, x, &incx);
}
static inline
void scal( const int& n, const double& alpha, double *x, const int& incx)
{
    dscal_(&n, &alpha, x, &incx);
}
static inline
void scal( const int& n, const std::complex<float>& alpha, std::complex<float> *x, const int& incx)
{
    cscal_(&n, &alpha, x, &incx);
}
static inline
void scal( const int& n, const std::complex<double>& alpha, std::complex<double> *x, const int& incx)
{
    zscal_(&n, &alpha, x, &incx);
}

// Peize Lin add 2017-10-27
// d=x*y
static inline
float dot( const int& n, const float *x, const int& incx, const float *y, const int& incy)
{
    return sdot_(&n, x, &incx, y, &incy);
}
static inline
double dot( const int& n, const double *x, const int& incx, const double *y, const int& incy)
{
    return ddot_(&n, x, &incx, y, &incy);
}
// Denghui Lu add 2023-8-01
static inline
std::complex<float> dot(const int& n, const std::complex<float> *x, const int& incx, const std::complex<float> *y, const int& incy)
{
    std::complex<float> result = {0, 0};
    // cdotc_(&n, x, &incx, y, &incy, &result);
    for (int ii = 0; ii < n; ii++) {
        result += std::conj(x[ii * incx]) * y[ii * incy];
    }
    return result;
}
static inline
std::complex<double> dot(const int& n, const std::complex<double> *x, const int& incx, const std::complex<double> *y, const int& incy)
{
    std::complex<double> result = {0, 0};
    // zdotc_(&n, x, &incx, y, &incy, &result);
    for (int ii = 0; ii < n; ii++) {
        result += std::conj(x[ii * incx]) * y[ii * incy];
    }
    return result;
}

// Peize Lin add 2017-10-27, fix bug trans 2019-01-17
// C = a * A.? * B.? + b * C
static inline
void gemm(const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const float& alpha, const float* A, const int& lda, const float* B, const int& ldb,
          const float& beta, float* C, const int& ldc)
{
    sgemm_(&transa, &transb, &m, &n, &k,
           &alpha, A, &lda, B, &ldb,
           &beta, C, &ldc);
}
static inline
void gemm(const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const double& alpha, const double* A, const int& lda, const double* B, const int& ldb,
          const double& beta, double* C, const int& ldc)
{
    dgemm_(&transa, &transb, &m, &n, &k,
           &alpha, A, &lda, B, &ldb,
           &beta, C, &ldc);
}
static inline
void gemm(const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<float>& alpha, const std::complex<float>* A, const int& lda, const std::complex<float>* B, const int& ldb,
          const std::complex<float>& beta, std::complex<float>* C, const int& ldc)
{
    cgemm_(&transa, &transb, &m, &n, &k,
           &alpha, A, &lda, B, &ldb,
           &beta, C, &ldc);
}
static inline
void gemm(const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const std::complex<double>& alpha, const std::complex<double>* A, const int& lda, const std::complex<double>* B, const int& ldb,
          const std::complex<double>& beta, std::complex<double>* C, const int& ldc)
{
    zgemm_(&transa, &transb, &m, &n, &k,
           &alpha, A, &lda, B, &ldb,
           &beta, C, &ldc);
}

template <typename T>
static inline
void gemm_batched(const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const T& alpha, T** A, const int& lda, T** B, const int& ldb,
          const T& beta, T** C, const int& ldc, const int& batch_size)
{
    for (int ii = 0; ii < batch_size; ++ii) {
        // Call the single GEMV for each pair of matrix A[ii] and vector x[ii]
        BlasConnector::gemm(transa, transb, m, n, k, alpha, A[ii], lda, B[ii], ldb, beta, C[ii], ldc);
    }
}

template <typename T>
static inline
void gemm_batched_strided(const char& transa, const char& transb, const int& m, const int& n, const int& k,
          const T& alpha, const T* A, const int& lda, const int& stride_a, const T* B, const int& ldb, const int& stride_b,
          const T& beta, T* C, const int& ldc, const int& stride_c, const int& batch_size)
{
    for (int ii = 0; ii < batch_size; ii++) {
        // Call the single GEMV for each pair of matrix A[ii] and vector x[ii]
        BlasConnector::gemm(transa, transb, m, n, k, alpha, A + ii * stride_a, lda, B + ii * stride_b, ldb, beta, C + ii * stride_c, ldc);
    }
}

static inline
void gemv(const char& trans, const int& m, const int& n,
          const float& alpha, const float *A, const int& lda, const float *x, const int& incx,
          const float& beta, float *y, const int& incy)
{
    sgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}
static inline
void gemv(const char& trans, const int& m, const int& n,
          const double& alpha, const double *A, const int& lda, const double *x, const int& incx,
          const double& beta, double *y, const int& incy)
{
    dgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}
static inline
void gemv(const char& trans, const int& m, const int& n,
          const std::complex<float>& alpha, const std::complex<float> *A, const int& lda, const std::complex<float> *x, const int& incx,
          const std::complex<float>& beta, std::complex<float> *y, const int& incy)
{
    cgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}
static inline
void gemv(const char& trans, const int& m, const int& n,
          const std::complex<double>& alpha, const std::complex<double> *A, const int& lda, const std::complex<double> *x, const int& incx,
          const std::complex<double>& beta, std::complex<double> *y, const int& incy)
{
    zgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

template <typename T>
static inline
void gemv_batched(const char& trans, const int& m, const int& n,
          const T& alpha, T** A, const int& lda, T** x, const int& incx,
          const T& beta, T** y, const int& incy, const int& batch_size)
{
    for (int ii = 0; ii < batch_size; ++ii) {
        // Call the single GEMV for each pair of matrix A[ii] and vector x[ii]
        BlasConnector::gemv(trans, m, n, alpha, A[ii], lda, x[ii], incy, beta, y[ii], incy);
    }
}

template <typename T>
static inline
void gemv_batched_strided(const char& transa, const int& m, const int& n,
          const T& alpha, const T* A, const int& lda, const int& stride_a, const T* x, const int& incx, const int& stride_x,
          const T& beta, T* y, const int& incy, const int& stride_y, const int& batch_size)
{
    for (int ii = 0; ii < batch_size; ii++) {
        // Call the single GEMV for each pair of matrix A[ii] and vector x[ii]
        BlasConnector::gemv(transa, m, n, alpha, A + ii * stride_a, lda, x + ii * stride_x, incx, beta, y + ii * stride_y, incy);
    }
}

// Peize Lin add 2018-06-12
// out = ||x||_2
static inline
float nrm2( const int n, const float *x, const int incx )
{
    return snrm2_( &n, x, &incx );
}
static inline
double nrm2( const int n, const double *x, const int incx )
{
    return dnrm2_( &n, x, &incx );
}
static inline
double nrm2( const int n, const std::complex<double> *x, const int incx )
{
    return dznrm2_( &n, x, &incx );
}

// copies a into b
static inline
void copy(const long n, const double *a, const int incx, double *b, const int incy)
{
    dcopy_(&n, a, &incx, b, &incy);
}
static inline
void copy(const long n, const std::complex<double> *a, const int incx, std::complex<double> *b, const int incy)
{
    zcopy_(&n, a, &incx, b, &incy);
}

} // namespace BlasConnector
} // namespace container

#endif // BASE_THIRD_PARTY_BLAS_H_
