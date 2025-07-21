/* level 2: matrix-std::vector operations, O(n^2) data and O(n^2) work.
 * This file contains the implementation of the BLAS level 2 operations.
 * These operations include matrix-vector multiplication and related operations.
 */
#include "blas_connector.h"
#include "macros.h"
#include <cblas.h>

#ifdef __DSP
#include "source_base/kernels/dsp/dsp_connector.h"
#include "source_base/global_variable.h"
#endif

#ifdef __CUDA
#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#endif

extern "C"
{
    // level 2: matrix-std::vector operations, O(n^2) data and O(n^2) work.
	void sgemv_(const char*const transa, const int*const m, const int*const n,
		const float*const alpha, const float*const a, const int*const lda, const float*const x, const int*const incx,
		const float*const beta, float*const y, const int*const incy);
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
}

void BlasConnector::gemv_cm(const char trans, const int m, const int n,
    const float alpha, const float* A, const int lda, const float* X, const int incx,
    const float beta, float* Y, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	sgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasOperation_t cutransA = BlasUtils::judge_trans(false, trans, "gemv_op");
		cublasErrcheck(cublasSgemv(BlasUtils::cublas_handle, cutransA, m, n, &alpha, A, lda, X, incx, &beta, Y, incy));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::gemv_cm(const char trans, const int m, const int n,
    const double alpha, const double* A, const int lda, const double* X, const int incx,
    const double beta, double* Y, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	dgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasOperation_t cutransA = BlasUtils::judge_trans(false, trans, "gemv_op");
		cublasErrcheck(cublasDgemv(BlasUtils::cublas_handle, cutransA, m, n, &alpha, A, lda, X, incx, &beta, Y, incy));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::gemv_cm(const char trans, const int m, const int n,
    const std::complex<float> alpha, const std::complex<float> *A, const int lda, const std::complex<float> *X, const int incx,
    const std::complex<float> beta, std::complex<float> *Y, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	cgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cuFloatComplex alpha_cu = make_cuFloatComplex(alpha.real(), alpha.imag());
    	cuFloatComplex beta_cu = make_cuFloatComplex(beta.real(), beta.imag());
		cublasOperation_t cutransA = BlasUtils::judge_trans(true, trans, "gemv_op");
		cublasErrcheck(cublasCgemv(BlasUtils::cublas_handle, cutransA, m, n, &alpha_cu, (cuFloatComplex*)A, lda, (cuFloatComplex*)X, incx, &beta_cu, (cuFloatComplex*)Y, incy));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::gemv_cm(const char trans, const int m, const int n,
    const std::complex<double> alpha, const std::complex<double> *A, const int lda, const std::complex<double> *X, const int incx,
    const std::complex<double> beta, std::complex<double> *Y, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	zgemv_(&trans, &m, &n, &alpha, A, &lda, X, &incx, &beta, Y, &incy);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cuDoubleComplex alpha_cu = make_cuDoubleComplex(alpha.real(), alpha.imag());
    	cuDoubleComplex beta_cu = make_cuDoubleComplex(beta.real(), beta.imag());
		cublasOperation_t cutransA = BlasUtils::judge_trans(true, trans, "gemv_op");
		cublasErrcheck(cublasZgemv(BlasUtils::cublas_handle, cutransA, m, n, &alpha_cu, (cuDoubleComplex*)A, lda, (cuDoubleComplex*)X, incx, &beta_cu, (cuDoubleComplex*)Y, incy));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::ger_cm(int m, int n, double alpha, const double* x,
    int incx, const double* y, const int incy, double a, int lda, base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice) {
        dger_(&m, &n, &alpha, x, &incx, y, &incy, &a, &lda);
    }
    else {
        throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
    }
}