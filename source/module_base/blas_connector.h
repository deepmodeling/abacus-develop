#ifndef BLAS_CONNECTOR_H
#define BLAS_CONNECTOR_H

#include <complex>

// Class BlasConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	BlasConnector::functionname(parameter list).
class BlasConnector
{
public:

	// Peize Lin add 2016-08-04
	// y=a*x+y
	static inline
	void axpy( const int n, const float alpha, const float *X, const int incX, float *Y, const int incY);

	static inline
	void axpy( const int n, const double alpha, const double *X, const int incX, double *Y, const int incY);

	static inline
	void axpy( const int n, const std::complex<float> alpha, const std::complex<float> *X, const int incX, std::complex<float> *Y, const int incY);

	static inline
	void axpy( const int n, const std::complex<double> alpha, const std::complex<double> *X, const int incX, std::complex<double> *Y, const int incY);


	// Peize Lin add 2016-08-04
	// x=a*x
	static inline
	void scal( const int n,  const float alpha, float *X, const int incX);

	static inline
	void scal( const int n, const double alpha, double *X, const int incX);

	static inline
	void scal( const int n, const std::complex<float> alpha, std::complex<float> *X, const int incX);

	static inline
	void scal( const int n, const std::complex<double> alpha, std::complex<double> *X, const int incX);


	// Peize Lin add 2017-10-27
	// d=x*y
	static inline
	float dot( const int n, const float *X, const int incX, const float *Y, const int incY);

	static inline
	double dot( const int n, const double *X, const int incX, const double *Y, const int incY);


	// Peize Lin add 2017-10-27, fix bug trans 2019-01-17
	// C = a * A.? * B.? + b * C
	static inline
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const float alpha, const float *a, const int lda, const float *b, const int ldb,
		const float beta, float *c, const int ldc);

	static inline
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const double alpha, const double *a, const int lda, const double *b, const int ldb,
		const double beta, double *c, const int ldc);

    static inline
    void gemm(const char transa, const char transb, const int m, const int n, const int k,
              const std::complex<float> alpha, const std::complex<float> *a, const int lda, const std::complex<float> *b, const int ldb,
              const std::complex<float> beta, std::complex<float> *c, const int ldc);

	static inline
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb,
		const std::complex<double> beta, std::complex<double> *c, const int ldc);

    static inline
    void gemv(const char trans, const int m, const int n,
        const double alpha, const double* A, const int lda, const double* X, const int incx,
        const double beta, double* Y, const int incy);

    static inline
    void gemv(const char trans, const int m, const int n,
          const std::complex<float> alpha, const std::complex<float> *A, const int lda, const std::complex<float> *X, const int incx,
          const std::complex<float> beta, std::complex<float> *Y, const int incy);

    static inline
    void gemv(const char trans, const int m, const int n,
              const std::complex<double> alpha, const std::complex<double> *A, const int lda, const std::complex<double> *X, const int incx,
              const std::complex<double> beta, std::complex<double> *Y, const int incy);
 

	// Peize Lin add 2018-06-12
	// out = ||x||_2
	static inline
	float nrm2( const int n, const float *X, const int incX );

	static inline
	double nrm2( const int n, const double *X, const int incX );

	static inline
	double nrm2( const int n, const std::complex<double> *X, const int incX );


	// copies a into b
	static inline
	void copy(const long n, const double *a, const int incx, double *b, const int incy);

	static inline
	void copy(const long n, const std::complex<double> *a, const int incx, std::complex<double> *b, const int incy);
};

// If GATHER_INFO is defined, the original function is replaced with a "i" suffix,
// preventing changes on the original code.
// The real function call is at gather_math_lib_info.cpp
#ifdef GATHER_INFO

#define zgemm_ zgemm_i
void zgemm_i(const char *transa,
             const char *transb,
             const int *m,
             const int *n,
             const int *k,
             const std::complex<double> *alpha,
             const std::complex<double> *a,
             const int *lda,
             const std::complex<double> *b,
             const int *ldb,
             const std::complex<double> *beta,
             std::complex<double> *c,
             const int *ldc);

#define zaxpy_  zaxpy_i
void zaxpy_i(const int *N,
            const std::complex<double> *alpha,
            const std::complex<double> *X,
            const int *incX,
            std::complex<double> *Y,
            const int *incY);

/*
#define zgemv_ zgemv_i

void zgemv_i(const char *trans,
             const int *m,
             const int *n,
             const std::complex<double> *alpha,
             const std::complex<double> *a,
             const int *lda,
             const std::complex<double> *x,
             const int *incx,
             const std::complex<double> *beta,
             std::complex<double> *y,
             const int *incy);
*/

#endif // GATHER_INFO
#endif // BLAS_CONNECTOR_H
