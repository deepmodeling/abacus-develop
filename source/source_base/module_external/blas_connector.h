#ifndef BLAS_CONNECTOR_H
#define BLAS_CONNECTOR_H

#include <complex>
#include "source_base/module_device/types.h"
#include "../macros.h"

// These still need to be linked in the header file
// Because quite a lot of code will directly use the original cblas kernels.

// The following declarations cover only a subset of BLAS routines.
// If you need a BLAS function that is not included here, feel free to add its declaration as needed.
extern "C"
{
// Level 1 BLAS
void sscal_(const int *N, const float *alpha, float *X, const int *incX);
void dscal_(const int *N, const double *alpha, double *X, const int *incX);
void cscal_(const int *N, const std::complex<float> *alpha, std::complex<float> *X, const int *incX);
void zscal_(const int *N, const std::complex<double> *alpha, std::complex<double> *X, const int *incX);

void saxpy_(const int *N, const float *alpha, const float *X, const int *incX, float *Y, const int *incY);
void daxpy_(const int *N, const double *alpha, const double *X, const int *incX, double *Y, const int *incY);
void caxpy_(const int *N, const std::complex<float> *alpha, const std::complex<float> *X, const int *incX, std::complex<float> *Y, const int *incY);
void zaxpy_(const int *N, const std::complex<double> *alpha, const std::complex<double> *X, const int *incX, std::complex<double> *Y, const int *incY);

void scopy_(const int *n, const float *a, const int *incx, float *b, const int *incy);
void dcopy_(const int *n, const double *a, const int *incx, double *b, const int *incy);
void ccopy_(const int *n, const std::complex<float> *a, const int *incx, std::complex<float> *b, const int *incy);
void zcopy_(const int *n, const std::complex<double> *a, const int *incx, std::complex<double> *b, const int *incy);

// Note: sdot_/ddot_ return by value
float sdot_(const int *N, const float *X, const int *incX, const float *Y, const int *incY);
double ddot_(const int *N, const double *X, const int *incX, const double *Y, const int *incY);

float snrm2_(const int *n, const float *X, const int *incX);
double dnrm2_(const int *n, const double *X, const int *incX);
float scnrm2_(const int *n, const std::complex<float> *X, const int *incX);
double dznrm2_(const int *n, const std::complex<double> *X, const int *incX);

// Level 2 BLAS

void sgemv_(const char *transa, const int *m, const int *n,
            const float *alpha, const float *a, const int *lda,
            const float *x, const int *incx,
            const float *beta, float *y, const int *incy);

void dgemv_(const char *transa, const int *m, const int *n,
            const double *alpha, const double *a, const int *lda,
            const double *x, const int *incx,
            const double *beta, double *y, const int *incy);

void cgemv_(const char *trans, const int *m, const int *n,
            const std::complex<float> *alpha,
            const std::complex<float> *a, const int *lda,
            const std::complex<float> *x, const int *incx,
            const std::complex<float> *beta,
            std::complex<float> *y, const int *incy);

void zgemv_(const char *trans, const int *m, const int *n,
            const std::complex<double> *alpha,
            const std::complex<double> *a, const int *lda,
            const std::complex<double> *x, const int *incx,
            const std::complex<double> *beta,
            std::complex<double> *y, const int *incy);

void dsymv_(const char *uplo, const int *n,
            const double *alpha, const double *a, const int *lda,
            const double *x, const int *incx,
            const double *beta, double *y, const int *incy);

void dger_(const int *m, const int *n,
           const double *alpha,
           const double *x, const int *incx,
           const double *y, const int *incy,
           double *a, const int *lda);

void zgerc_(const int *m, const int *n,
            const std::complex<double> *alpha,
            const std::complex<double> *x, const int *incx,
            const std::complex<double> *y, const int *incy,
            std::complex<double> *a, const int *lda);

// Level 3 BLAS

void sgemm_(const char *transa, const char *transb,
            const int *m, const int *n, const int *k,
            const float *alpha,
            const float *a, const int *lda,
            const float *b, const int *ldb,
            const float *beta,
            float *c, const int *ldc);

void dgemm_(const char *transa, const char *transb,
            const int *m, const int *n, const int *k,
            const double *alpha,
            const double *a, const int *lda,
            const double *b, const int *ldb,
            const double *beta,
            double *c, const int *ldc);

void cgemm_(const char *transa, const char *transb,
            const int *m, const int *n, const int *k,
            const std::complex<float> *alpha,
            const std::complex<float> *a, const int *lda,
            const std::complex<float> *b, const int *ldb,
            const std::complex<float> *beta,
            std::complex<float> *c, const int *ldc);

void zgemm_(const char *transa, const char *transb,
            const int *m, const int *n, const int *k,
            const std::complex<double> *alpha,
            const std::complex<double> *a, const int *lda,
            const std::complex<double> *b, const int *ldb,
            const std::complex<double> *beta,
            std::complex<double> *c, const int *ldc);

void ssymm_(const char *side, const char *uplo,
            const int *m, const int *n,
            const float *alpha,
            const float *a, const int *lda,
            const float *b, const int *ldb,
            const float *beta,
            float *c, const int *ldc);

void dsymm_(const char *side, const char *uplo,
            const int *m, const int *n,
            const double *alpha,
            const double *a, const int *lda,
            const double *b, const int *ldb,
            const double *beta,
            double *c, const int *ldc);

void csymm_(const char *side, const char *uplo,
            const int *m, const int *n,
            const std::complex<float> *alpha,
            const std::complex<float> *a, const int *lda,
            const std::complex<float> *b, const int *ldb,
            const std::complex<float> *beta,
            std::complex<float> *c, const int *ldc);

void zsymm_(const char *side, const char *uplo,
            const int *m, const int *n,
            const std::complex<double> *alpha,
            const std::complex<double> *a, const int *lda,
            const std::complex<double> *b, const int *ldb,
            const std::complex<double> *beta,
            std::complex<double> *c, const int *ldc);

void chemm_(const char *side, const char *uplo,
            const int *m, const int *n,
            const std::complex<float> *alpha,
            const std::complex<float> *a, const int *lda,
            const std::complex<float> *b, const int *ldb,
            const std::complex<float> *beta,
            std::complex<float> *c, const int *ldc);

void zhemm_(const char *side, const char *uplo,
            const int *m, const int *n,
            const std::complex<double> *alpha,
            const std::complex<double> *a, const int *lda,
            const std::complex<double> *b, const int *ldb,
            const std::complex<double> *beta,
            std::complex<double> *c, const int *ldc);

void dtrsm_(const char *side, const char *uplo, const char *transa, const char *diag,
            const int *m, const int *n,
            const double *alpha,
            const double *a, const int *lda,
            double *b, const int *ldb);

void ztrsm_(const char *side, const char *uplo, const char *transa, const char *diag,
            const int *m, const int *n,
            const std::complex<double> *alpha,
            const std::complex<double> *a, const int *lda,
            std::complex<double> *b, const int *ldb);

// === Hermitian rank-k update ===
void cherk_(const char* uplo, const char* trans, const int* n, const int* k,
            const float* alpha,
            const std::complex<float>* a, const int* lda,
            const float* beta,
            std::complex<float>* c, const int* ldc);

void zherk_(const char* uplo, const char* trans, const int* n, const int* k,
            const double* alpha,
            const std::complex<double>* a, const int* lda,
            const double* beta,
            std::complex<double>* c, const int* ldc);

// === Symmetric rank-k update ===
void dsyrk_(const char* uplo, const char* trans, const int* n, const int* k,
            const double* alpha,
            const double* a, const int* lda,
            const double* beta, 
            double* c,
            const int* ldc);
}

// Class BlasConnector provide the connector to fortran lapack routine.
// The entire function in this class are static and inline function.
// Usage example:	BlasConnector::functionname(parameter list).
class BlasConnector
{
public:

	// Peize Lin add 2016-08-04
	// y=a*x+y
	static
	void axpy( const int n, const float alpha, const float *X, const int incX, float *Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void axpy( const int n, const double alpha, const double *X, const int incX, double *Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void axpy( const int n, const std::complex<float> alpha, const std::complex<float> *X, const int incX, std::complex<float> *Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void axpy( const int n, const std::complex<double> alpha, const std::complex<double> *X, const int incX, std::complex<double> *Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);


	// Peize Lin add 2016-08-04
	// x=a*x
	static
	void scal( const int n,  const float alpha, float *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void scal( const int n, const double alpha, double *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void scal( const int n, const std::complex<float> alpha, std::complex<float> *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void scal( const int n, const std::complex<double> alpha, std::complex<double> *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);


	// Peize Lin add 2017-10-27
	// d=x*y
	static
	float dot( const int n, const float*const X, const int incX, const float*const Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	double dot( const int n, const double*const X, const int incX, const double*const Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	// d=x*y
	static
	float dotu( const int n, const float*const X, const int incX, const float*const Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	double dotu( const int n, const double*const X, const int incX, const double*const Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	std::complex<float> dotu( const int n, const std::complex<float>*const X, const int incX, const std::complex<float>*const Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	std::complex<double> dotu( const int n, const std::complex<double>*const X, const int incX, const std::complex<double>*const Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	// d=x.conj()*y
	static
	float dotc( const int n, const float*const X, const int incX, const float*const Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	double dotc( const int n, const double*const X, const int incX, const double*const Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	std::complex<float> dotc( const int n, const std::complex<float>*const X, const int incX, const std::complex<float>*const Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	std::complex<double> dotc( const int n, const std::complex<double>*const X, const int incX, const std::complex<double>*const Y, const int incY, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	// Peize Lin add 2017-10-27, fix bug trans 2019-01-17
	// C = a * A.? * B.? + b * C
	// Row Major by default
	static
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const float alpha, const float *a, const int lda, const float *b, const int ldb,
		const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const double alpha, const double *a, const int lda, const double *b, const int ldb,
		const double beta, double *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

    static
    void gemm(const char transa, const char transb, const int m, const int n, const int k,
              const std::complex<float> alpha, const std::complex<float> *a, const int lda, const std::complex<float> *b, const int ldb,
              const std::complex<float> beta, std::complex<float> *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void gemm(const char transa, const char transb, const int m, const int n, const int k,
		const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb,
		const std::complex<double> beta, std::complex<double> *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	// Col-Major if you need to use it

	static
	void gemm_cm(const char transa, const char transb, const int m, const int n, const int k,
		const float alpha, const float *a, const int lda, const float *b, const int ldb,
		const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void gemm_cm(const char transa, const char transb, const int m, const int n, const int k,
		const double alpha, const double *a, const int lda, const double *b, const int ldb,
		const double beta, double *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

    static
    void gemm_cm(const char transa, const char transb, const int m, const int n, const int k,
              const std::complex<float> alpha, const std::complex<float> *a, const int lda, const std::complex<float> *b, const int ldb,
              const std::complex<float> beta, std::complex<float> *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void gemm_cm(const char transa, const char transb, const int m, const int n, const int k,
		const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb,
		const std::complex<double> beta, std::complex<double> *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	// side=='L': C = a * A * B + b * C.
	// side=='R': C = a * B * A + b * C.
	// 		A == A^T
	// Because you cannot pack symm or hemm into a row-major kernel by exchanging parameters, so only col-major functions are provided.
	static
	void symm_cm(const char side, const char uplo, const int m, const int n,
		const float alpha, const float *a, const int lda, const float *b, const int ldb,
		const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void symm_cm(const char side, const char uplo, const int m, const int n,
		const double alpha, const double *a, const int lda, const double *b, const int ldb,
		const double beta, double *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

    static
    void symm_cm(const char side, const char uplo, const int m, const int n,
              const std::complex<float> alpha, const std::complex<float> *a, const int lda, const std::complex<float> *b, const int ldb,
              const std::complex<float> beta, std::complex<float> *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void symm_cm(const char side, const char uplo, const int m, const int n,
		const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb,
		const std::complex<double> beta, std::complex<double> *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	// side=='L': C = a * A * B + b * C.
	// side=='R': C = a * B * A + b * C.
	// 		A == A^H
	static
	void hemm_cm(const char side, const char uplo, const int m, const int n,
		const float alpha, const float *a, const int lda, const float *b, const int ldb,
		const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void hemm_cm(const char side, const char uplo, const int m, const int n,
		const double alpha, const double *a, const int lda, const double *b, const int ldb,
		const double beta, double *c, const int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
    void hemm_cm(char side, char uplo, int m, int n,
            std::complex<float> alpha, std::complex<float> *a, int lda, std::complex<float> *b, int ldb,
            std::complex<float> beta, std::complex<float> *c, int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void hemm_cm(char side, char uplo, int m, int n,
		std::complex<double> alpha, std::complex<double> *a, int lda, std::complex<double> *b, int ldb,
		std::complex<double> beta, std::complex<double> *c, int ldc, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	// y = A*x + beta*y
	static
	void gemv(const char trans, const int m, const int n,
        const float alpha, const float* A, const int lda, const float* X, const int incx,
        const float beta, float* Y, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

    static
    void gemv(const char trans, const int m, const int n,
        const double alpha, const double* A, const int lda, const double* X, const int incx,
        const double beta, double* Y, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

    static
    void gemv(const char trans, const int m, const int n,
          const std::complex<float> alpha, const std::complex<float> *A, const int lda, const std::complex<float> *X, const int incx,
          const std::complex<float> beta, std::complex<float> *Y, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

    static
    void gemv(const char trans, const int m, const int n,
              const std::complex<double> alpha, const std::complex<double> *A, const int lda, const std::complex<double> *X, const int incx,
              const std::complex<double> beta, std::complex<double> *Y, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	// Peize Lin add 2018-06-12
	// out = ||x||_2
	static
	float nrm2( const int n, const float *X, const int, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice );

	static
	double nrm2( const int n, const double *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice );

	static
	double nrm2( const int n, const std::complex<double> *X, const int incX, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice );


	// copies a into b
	static
	void copy(const int n, const double *a, const int incx, double *b, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void copy(const int n, const float *a, const int incx, float *b, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void copy(const int n, const std::complex<float> *a, const int incx, std::complex<float> *b, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void copy(const int n, const std::complex<double> *a, const int incx, std::complex<double> *b, const int incy, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	// There is some other operators needed, so implemented manually here
	template <typename T>
	static
	void vector_mul_vector(const int& dim, T* result, const T* vector1, const T* vector2, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	template <typename T>
	static
	void vector_div_vector(const int& dim, T* result, const T* vector1, const T* vector2, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	// y = alpha * x + beta * y
	static
	void vector_add_vector(const int& dim, float *result, const float *vector1, const float constant1, const float *vector2, const float constant2, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void vector_add_vector(const int& dim, double *result, const double *vector1, const double constant1, const double *vector2, const double constant2, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void vector_add_vector(const int& dim, std::complex<float> *result, const std::complex<float> *vector1, const float constant1, const std::complex<float> *vector2, const float constant2, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);

	static
	void vector_add_vector(const int& dim, std::complex<double> *result, const std::complex<double> *vector1, const double constant1, const std::complex<double> *vector2, const double constant2, base_device::AbacusDevice_t device_type = base_device::AbacusDevice_t::CpuDevice);
};

#ifdef __CUDA

#include <cuda_runtime.h>
#include "cublas_v2.h"

// If you want to use cublas, you need these functions to create and destroy the cublas/hipblas handle.
// You also need to use these functions to translate the transpose parameter into cublas/hipblas datatype.

namespace BlasUtils{

	static cublasHandle_t cublas_handle = nullptr;

	void createGpuBlasHandle(); // Create a cublas/hipblas handle.

	void destoryBLAShandle(); // Destroy the cublas/hipblas handle. Do this when the software is about to end.

	cublasOperation_t judge_trans(bool is_complex, const char& trans, const char* name); // Translate a normal transpose parameter to a cublas/hipblas type.

	cublasSideMode_t judge_side(const char& trans); // Translate a normal side parameter to a cublas/hipblas type.

	cublasFillMode_t judge_fill(const char& trans); // Translate a normal fill parameter to a cublas/hipblas type.

}

#endif

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
