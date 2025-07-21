/* level 3: matrix-matrix operations, O(n^2) data and O(n^3) work.
 * This file contains the implementation of the BLAS level 3 operations.
 * These operations include matrix-matrix multiplication and related operations.
 */
#include "blas_connector.h"
#include "../macros.h"

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

	// A is symmetric. C = a * A.? * B.? + b * C
	void ssymm_(const char *side, const char *uplo, const int *m, const int *n,
		const float *alpha, const float *a, const int *lda, const float *b, const int *ldb,
		const float *beta, float *c, const int *ldc);
	void dsymm_(const char *side, const char *uplo, const int *m, const int *n,
		const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);
	void csymm_(const char *side, const char *uplo, const int *m, const int *n,
		const std::complex<float> *alpha, const std::complex<float> *a, const int *lda, const std::complex<float> *b, const int *ldb,
		const std::complex<float> *beta, std::complex<float> *c, const int *ldc);
	void zsymm_(const char *side, const char *uplo, const int *m, const int *n,
		const std::complex<double> *alpha, const std::complex<double> *a, const int *lda, const std::complex<double> *b, const int *ldb,
		const std::complex<double> *beta, std::complex<double> *c, const int *ldc);

	// A is hermitian. C = a * A.? * B.? + b * C
	void chemm_(char *side, char *uplo, int *m, int *n,std::complex<float> *alpha,
		std::complex<float> *a,  int *lda,  std::complex<float> *b, int *ldb, std::complex<float> *beta, std::complex<float> *c, int *ldc);
	void zhemm_(char *side, char *uplo, int *m, int *n,std::complex<double> *alpha,
		std::complex<double> *a,  int *lda,  std::complex<double> *b, int *ldb, std::complex<double> *beta, std::complex<double> *c, int *ldc);

    // symmetric rank-k update
    void dsyrk_(
        const char* uplo,
        const char* trans,
        const int* n,
        const int* k,
        const double* alpha,
        const double* a,
        const int* lda,
        const double* beta,
        double* c,
        const int* ldc
    );
}

// C = a * A.? * B.? + b * C
// Row-Major part
void BlasConnector::gemm(const char transa, const char transb, const int m, const int n, const int k,
	const float alpha, const float *a, const int lda, const float *b, const int ldb,
	const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		sgemm_(&transb, &transa, &n, &m, &k,
		&alpha, b, &ldb, a, &lda,
		&beta, c, &ldc);
	}
#ifdef __DSP
	else if (device_type == base_device::AbacusDevice_t::DspDevice){
		mtfunc::sgemm_mth_(&transb, &transa, &n, &m, &k,
		&alpha, b, &ldb, a, &lda,
		&beta, c, &ldc, GlobalV::MY_RANK);
	}
#endif
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasOperation_t cutransA = BlasUtils::judge_trans(false, transa, "gemm_op");
		cublasOperation_t cutransB = BlasUtils::judge_trans(false, transb, "gemm_op");
		cublasErrcheck(cublasSgemm(BlasUtils::cublas_handle, cutransA, cutransB, n, m, k, &alpha, b, ldb, a, lda, &beta, c, ldc));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::gemm(const char transa,
                         const char transb,
                         const int m,
                         const int n,
                         const int k,
                         const double alpha,
                         const double* a,
                         const int lda,
                         const double* b,
                         const int ldb,
                         const double beta,
                         double* c,
                         const int ldc,
                         base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        dgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
    }
#ifdef __DSP
    else if (device_type == base_device::AbacusDevice_t::DspDevice)
    {
        mtfunc::dgemm_mth_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc, GlobalV::MY_RANK);
    }
#endif
    else if (device_type == base_device::AbacusDevice_t::GpuDevice)
    {
#ifdef __CUDA
        cublasOperation_t cutransA = BlasUtils::judge_trans(false, transa, "gemm_op");
        cublasOperation_t cutransB = BlasUtils::judge_trans(false, transb, "gemm_op");
        cublasErrcheck(
            cublasDgemm(BlasUtils::cublas_handle, cutransA, cutransB, n, m, k, &alpha, b, ldb, a, lda, &beta, c, ldc));
#endif
    }
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::gemm(const char transa,
                         const char transb,
                         const int m,
                         const int n,
                         const int k,
                         const std::complex<float> alpha,
                         const std::complex<float>* a,
                         const int lda,
                         const std::complex<float>* b,
                         const int ldb,
                         const std::complex<float> beta,
                         std::complex<float>* c,
                         const int ldc,
                         base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        cgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
    }
#ifdef __DSP
    else if (device_type == base_device::AbacusDevice_t::DspDevice)
    {
        mtfunc::cgemm_mth_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc, GlobalV::MY_RANK);
    }
#endif
    else if (device_type == base_device::AbacusDevice_t::GpuDevice)
    {
#ifdef __CUDA
        cublasOperation_t cutransA = BlasUtils::judge_trans(false, transa, "gemm_op");
        cublasOperation_t cutransB = BlasUtils::judge_trans(false, transb, "gemm_op");
        cublasErrcheck(cublasCgemm(BlasUtils::cublas_handle,
                                   cutransA,
                                   cutransB,
                                   n,
                                   m,
                                   k,
                                   (float2*)&alpha,
                                   (float2*)b,
                                   ldb,
                                   (float2*)a,
                                   lda,
                                   (float2*)&beta,
                                   (float2*)c,
                                   ldc));
#endif
    }
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::gemm(const char transa,
                         const char transb,
                         const int m,
                         const int n,
                         const int k,
                         const std::complex<double> alpha,
                         const std::complex<double>* a,
                         const int lda,
                         const std::complex<double>* b,
                         const int ldb,
                         const std::complex<double> beta,
                         std::complex<double>* c,
                         const int ldc,
                         base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        zgemm_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc);
    }
#ifdef __DSP
    else if (device_type == base_device::AbacusDevice_t::DspDevice)
    {
        mtfunc::zgemm_mth_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc, GlobalV::MY_RANK);
    }
#endif
    else if (device_type == base_device::AbacusDevice_t::GpuDevice)
    {
#ifdef __CUDA
        cublasOperation_t cutransA = BlasUtils::judge_trans(false, transa, "gemm_op");
        cublasOperation_t cutransB = BlasUtils::judge_trans(false, transb, "gemm_op");
        cublasErrcheck(cublasZgemm(BlasUtils::cublas_handle,
                                   cutransA,
                                   cutransB,
                                   n,
                                   m,
                                   k,
                                   (double2*)&alpha,
                                   (double2*)b,
                                   ldb,
                                   (double2*)a,
                                   lda,
                                   (double2*)&beta,
                                   (double2*)c,
                                   ldc));
#endif
    }
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

// Col-Major part
void BlasConnector::gemm_cm(const char transa, const char transb, const int m, const int n, const int k,
	const float alpha, const float *a, const int lda, const float *b, const int ldb,
	const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		sgemm_(&transa, &transb, &m, &n, &k,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
#ifdef __DSP
	else if (device_type == base_device::AbacusDevice_t::DspDevice){
		mtfunc::sgemm_mth_(&transb, &transa, &m, &n, &k,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc, GlobalV::MY_RANK);
	}
#endif
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasOperation_t cutransA = BlasUtils::judge_trans(false, transa, "gemm_op");
		cublasOperation_t cutransB = BlasUtils::judge_trans(false, transb, "gemm_op");
		cublasErrcheck(cublasSgemm(BlasUtils::cublas_handle, cutransA, cutransB, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::gemm_cm(const char transa,
                            const char transb,
                            const int m,
                            const int n,
                            const int k,
                            const double alpha,
                            const double* a,
                            const int lda,
                            const double* b,
                            const int ldb,
                            const double beta,
                            double* c,
                            const int ldc,
                            base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    }
#ifdef __DSP
    else if (device_type == base_device::AbacusDevice_t::DspDevice)
    {
        mtfunc::dgemm_mth_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, GlobalV::MY_RANK);
    }
#endif
#ifdef __CUDA
    else if (device_type == base_device::AbacusDevice_t::GpuDevice)
    {
        cublasOperation_t cutransA = BlasUtils::judge_trans(false, transa, "gemm_op");
        cublasOperation_t cutransB = BlasUtils::judge_trans(false, transb, "gemm_op");
        cublasErrcheck(
            cublasDgemm(BlasUtils::cublas_handle, cutransA, cutransB, m, n, k, &alpha, a, lda, b, ldb, &beta, c, ldc));
    }
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::gemm_cm(const char transa,
                            const char transb,
                            const int m,
                            const int n,
                            const int k,
                            const std::complex<float> alpha,
                            const std::complex<float>* a,
                            const int lda,
                            const std::complex<float>* b,
                            const int ldb,
                            const std::complex<float> beta,
                            std::complex<float>* c,
                            const int ldc,
                            base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        cgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    }
#ifdef __DSP
    else if (device_type == base_device::AbacusDevice_t::DspDevice)
    {
        mtfunc::cgemm_mth_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, GlobalV::MY_RANK);
    }
#endif
#ifdef __CUDA
    else if (device_type == base_device::AbacusDevice_t::GpuDevice)
    {
        cublasOperation_t cutransA = BlasUtils::judge_trans(false, transa, "gemm_op");
        cublasOperation_t cutransB = BlasUtils::judge_trans(false, transb, "gemm_op");
        cublasErrcheck(cublasCgemm(BlasUtils::cublas_handle,
                                   cutransA,
                                   cutransB,
                                   m,
                                   n,
                                   k,
                                   (float2*)&alpha,
                                   (float2*)a,
                                   lda,
                                   (float2*)b,
                                   ldb,
                                   (float2*)&beta,
                                   (float2*)c,
                                   ldc));
    }
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::gemm_cm(const char transa,
                            const char transb,
                            const int m,
                            const int n,
                            const int k,
                            const std::complex<double> alpha,
                            const std::complex<double>* a,
                            const int lda,
                            const std::complex<double>* b,
                            const int ldb,
                            const std::complex<double> beta,
                            std::complex<double>* c,
                            const int ldc,
                            base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        zgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
    }
#ifdef __DSP
    else if (device_type == base_device::AbacusDevice_t::DspDevice)
    {
        mtfunc::zgemm_mth_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, GlobalV::MY_RANK);
    }
#endif
#ifdef __CUDA
    else if (device_type == base_device::AbacusDevice_t::GpuDevice)
    {
        cublasOperation_t cutransA = BlasUtils::judge_trans(false, transa, "gemm_op");
        cublasOperation_t cutransB = BlasUtils::judge_trans(false, transb, "gemm_op");
        cublasErrcheck(cublasZgemm(BlasUtils::cublas_handle,
                                   cutransA,
                                   cutransB,
                                   m,
                                   n,
                                   k,
                                   (double2*)&alpha,
                                   (double2*)a,
                                   lda,
                                   (double2*)b,
                                   ldb,
                                   (double2*)&beta,
                                   (double2*)c,
                                   ldc));
    }
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

// Symm and Hemm part. Only col-major is supported.

void BlasConnector::symm_cm(const char side, const char uplo, const int m, const int n,
	const float alpha, const float *a, const int lda, const float *b, const int ldb,
	const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		ssymm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasSideMode_t sideMode = BlasUtils::judge_side(side);
		cublasFillMode_t fillMode = BlasUtils::judge_fill(uplo);
		cublasErrcheck(cublasSsymm(BlasUtils::cublas_handle, sideMode, fillMode, m, n, &alpha, a, lda, b, ldb, &beta, c, ldc));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::symm_cm(const char side, const char uplo, const int m, const int n,
	const double alpha, const double *a, const int lda, const double *b, const int ldb,
	const double beta, double *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		dsymm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasSideMode_t sideMode = BlasUtils::judge_side(side);
		cublasFillMode_t fillMode = BlasUtils::judge_fill(uplo);
		cublasErrcheck(cublasDsymm(BlasUtils::cublas_handle, sideMode, fillMode, m, n, &alpha, a, lda, b, ldb, &beta, c, ldc));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::symm_cm(const char side, const char uplo, const int m, const int n,
    const std::complex<float> alpha, const std::complex<float> *a, const int lda, const std::complex<float> *b, const int ldb,
    const std::complex<float> beta, std::complex<float> *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	csymm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasSideMode_t sideMode = BlasUtils::judge_side(side);
		cublasFillMode_t fillMode = BlasUtils::judge_fill(uplo);
		cublasErrcheck(cublasCsymm(BlasUtils::cublas_handle, sideMode, fillMode, m, n, (float2*)&alpha, (float2*)a, lda, (float2*)b, ldb, (float2*)&beta, (float2*)c, ldc));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::symm_cm(const char side, const char uplo, const int m, const int n,
	const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb,
	const std::complex<double> beta, std::complex<double> *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zsymm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasSideMode_t sideMode = BlasUtils::judge_side(side);
		cublasFillMode_t fillMode = BlasUtils::judge_fill(uplo);
		cublasErrcheck(cublasZsymm(BlasUtils::cublas_handle, sideMode, fillMode, m, n, (double2*)&alpha, (double2*)a, lda, (double2*)b, ldb, (double2*)&beta, (double2*)c, ldc));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::hemm_cm(const char side, const char uplo, const int m, const int n,
	const float alpha, const float *a, const int lda, const float *b, const int ldb,
	const float beta, float *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	symm_cm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, device_type);
}

void BlasConnector::hemm_cm(const char side, const char uplo, const int m, const int n,
	const double alpha, const double *a, const int lda, const double *b, const int ldb,
	const double beta, double *c, const int ldc, base_device::AbacusDevice_t device_type)
{
	symm_cm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc, device_type);
}

void BlasConnector::hemm_cm(char side, char uplo, int m, int n,
    std::complex<float> alpha, std::complex<float> *a, int lda, std::complex<float> *b, int ldb,
    std::complex<float> beta, std::complex<float> *c, int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
    	chemm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasSideMode_t sideMode = BlasUtils::judge_side(side);
		cublasFillMode_t fillMode = BlasUtils::judge_fill(uplo);
		cublasErrcheck(cublasChemm(BlasUtils::cublas_handle, sideMode, fillMode, m, n, (float2*)&alpha, (float2*)a, lda, (float2*)b, ldb, (float2*)&beta, (float2*)c, ldc));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::hemm_cm(char side, char uplo, int m, int n,
	std::complex<double> alpha, std::complex<double> *a, int lda, std::complex<double> *b, int ldb,
	std::complex<double> beta, std::complex<double> *c, int ldc, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zhemm_(&side, &uplo, &m, &n,
		&alpha, a, &lda, b, &ldb,
		&beta, c, &ldc);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasSideMode_t sideMode = BlasUtils::judge_side(side);
		cublasFillMode_t fillMode = BlasUtils::judge_fill(uplo);
		cublasErrcheck(cublasZhemm(BlasUtils::cublas_handle, sideMode, fillMode, m, n, (double2*)&alpha, (double2*)a, lda, (double2*)b, ldb, (double2*)&beta, (double2*)c, ldc));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::syrk(char uplo, char trans, int n, int k,
    double alpha, const double* a, int lda, double beta, double* c, int ldc, base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        dsyrk_(&uplo, &trans, &n, &k, &alpha, a, &lda, &beta, c, &ldc);
    }
    else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::herk(char uplo, char trans, int n, int k, float alpha,
    const std::complex<float> *A, int lda, float beta, std::complex<float> *C, int ldc, base_device::AbacusDevice_t device_type)
{
    auto cblas_uplo = BlasUtils::toCblasUplo(uplo);
    auto cblas_trans = BlasUtils::toCblasTrans(trans);
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        cblas_cherk(CblasRowMajor, cblas_uplo, cblas_trans, n, k, alpha, A, lda, beta, C, ldc);
    }
    else
    {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
    }
}

void BlasConnector::herk(char uplo, char trans, int n, int k, double alpha,
    const std::complex<double> *A, int lda, double beta, std::complex<double> *C, int ldc, base_device::AbacusDevice_t device_type)
{
    auto cblas_uplo = BlasUtils::toCblasUplo(uplo);
    auto cblas_trans = BlasUtils::toCblasTrans(trans);
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        cblas_zherk(CblasRowMajor, cblas_uplo, cblas_trans, n, k, alpha, A, lda, beta, C, ldc);
    }
    else
    {
       	throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
    }
} 