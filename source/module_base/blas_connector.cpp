#include "module_base/blas_connector.h"

inline
void BlasConnector::axpy( const int n, const float alpha, const float *X, const int incX, float *Y, const int incY)
{
	saxpy_(&n, &alpha, X, &incX, Y, &incY);
}

inline
void BlasConnector::axpy( const int n, const double alpha, const double *X, const int incX, double *Y, const int incY)
{
	daxpy_(&n, &alpha, X, &incX, Y, &incY);
}

inline
void BlasConnector::axpy( const int n, const std::complex<float> alpha, const std::complex<float> *X, const int incX, std::complex<float> *Y, const int incY)
{
	caxpy_(&n, &alpha, X, &incX, Y, &incY);
}

inline
void BlasConnector::axpy( const int n, const std::complex<double> alpha, const std::complex<double> *X, const int incX, std::complex<double> *Y, const int incY)
{
	zaxpy_(&n, &alpha, X, &incX, Y, &incY);
}


// x=a*x
inline
void BlasConnector::scal( const int n,  const float alpha, float *X, const int incX)
{
	sscal_(&n, &alpha, X, &incX);
}

inline
void BlasConnector::scal( const int n, const double alpha, double *X, const int incX)
{
	dscal_(&n, &alpha, X, &incX);
}
	
inline
void BlasConnector::scal( const int n, const std::complex<float> alpha, std::complex<float> *X, const int incX)
{
	cscal_(&n, &alpha, X, &incX);
}

inline
void BlasConnector::scal( const int n, const std::complex<double> alpha, std::complex<double> *X, const int incX)
{
	zscal_(&n, &alpha, X, &incX);
}


// d=x*y
inline
float BlasConnector::dot( const int n, const float *X, const int incX, const float *Y, const int incY)
{
	return sdot_(&n, X, &incX, Y, &incY);
}
inline
double BlasConnector::dot( const int n, const double *X, const int incX, const double *Y, const int incY)
{
	return ddot_(&n, X, &incX, Y, &incY);
}

// C = a * A.? * B.? + b * C
inline
void BlasConnector::gemm(const char transa, const char transb, const int m, const int n, const int k,
	const float alpha, const float *a, const int lda, const float *b, const int ldb,
	const float beta, float *c, const int ldc)
{
	sgemm_(&transb, &transa, &n, &m, &k,
		&alpha, b, &ldb, a, &lda,
		&beta, c, &ldc);
}

inline
void BlasConnector::gemm(const char transa, const char transb, const int m, const int n, const int k,
	const double alpha, const double *a, const int lda, const double *b, const int ldb,
	const double beta, double *c, const int ldc)
{
	dgemm_(&transb, &transa, &n, &m, &k,
		&alpha, b, &ldb, a, &lda,
		&beta, c, &ldc);
}

inline
void BlasConnector::gemm(const char transa, const char transb, const int m, const int n, const int k,
    const std::complex<float> alpha, const std::complex<float> *a, const int lda, const std::complex<float> *b, const int ldb,
    const std::complex<float> beta, std::complex<float> *c, const int ldc)
{
    cgemm_(&transb, &transa, &n, &m, &k,
            &alpha, b, &ldb, a, &lda,
            &beta, c, &ldc);
}

inline
void BlasConnector::gemm(const char transa, const char transb, const int m, const int n, const int k,
	const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb,
	const std::complex<double> beta, std::complex<double> *c, const int ldc)
{
	zgemm_(&transb, &transa, &n, &m, &k,
		&alpha, b, &ldb, a, &lda,
		&beta, c, &ldc);
}

inline
void BlasConnector::gemv(const char trans, const int m, const int n,
    const double alpha, const double* A, const int lda, const double* X, const int incx,
    const double beta, double* Y, const int incy)
{
    dgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
}

inline
void BlasConnector::gemv(const char trans, const int m, const int n,
    const std::complex<float> alpha, const std::complex<float> *A, const int lda, const std::complex<float> *X, const int incx,
    const std::complex<float> beta, std::complex<float> *Y, const int incy)
{
    cgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
}

inline
void BlasConnector::gemv(const char trans, const int m, const int n,
    const std::complex<double> alpha, const std::complex<double> *A, const int lda, const std::complex<double> *X, const int incx,
    const std::complex<double> beta, std::complex<double> *Y, const int incy)
{
    zgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
}


// out = ||x||_2
inline
float BlasConnector::nrm2( const int n, const float *X, const int incX )
{
	return snrm2_( &n, X, &incX );
}

inline
double BlasConnector::nrm2( const int n, const double *X, const int incX )
{
	return dnrm2_( &n, X, &incX );
}

inline
double BlasConnector::nrm2( const int n, const std::complex<double> *X, const int incX )
{
	return dznrm2_( &n, X, &incX );
}

// copies a into b
inline
void BlasConnector::copy(const long n, const double *a, const int incx, double *b, const int incy)
{
	dcopy_(&n, a, &incx, b, &incy);
}

inline
void BlasConnector::copy(const long n, const std::complex<double> *a, const int incx, std::complex<double> *b, const int incy)
{
	zcopy_(&n, a, &incx, b, &incy);
}