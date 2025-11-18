#include "blas_connector.h"
#include "../macros.h"

#ifdef __DSP
#include "source_base/kernels/dsp/dsp_connector.h"
#include "source_base/global_variable.h"
#include "source_io/module_parameter/parameter.h"
#endif

#ifdef __CUDA
#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "source_base/kernels/math_kernel_op.h"
#include "source_base/module_device/memory_op.h"
#endif


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
		&beta, c, &ldc, GlobalV::MY_RANK % PARAM.inp.dsp_count);
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
        mtfunc::dgemm_mth_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc, GlobalV::MY_RANK % PARAM.inp.dsp_count);
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
        mtfunc::cgemm_mth_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc, GlobalV::MY_RANK % PARAM.inp.dsp_count);
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
        mtfunc::zgemm_mth_(&transb, &transa, &n, &m, &k, &alpha, b, &ldb, a, &lda, &beta, c, &ldc, GlobalV::MY_RANK % PARAM.inp.dsp_count);
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
		&beta, c, &ldc, GlobalV::MY_RANK % PARAM.inp.dsp_count);
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
        mtfunc::dgemm_mth_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, GlobalV::MY_RANK % PARAM.inp.dsp_count);
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
        mtfunc::cgemm_mth_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, GlobalV::MY_RANK % PARAM.inp.dsp_count);
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
        mtfunc::zgemm_mth_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, GlobalV::MY_RANK % PARAM.inp.dsp_count);
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

void BlasConnector::gemv(const char trans, const int m, const int n,
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

void BlasConnector::gemv(const char trans, const int m, const int n,
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

void BlasConnector::gemv(const char trans, const int m, const int n,
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

void BlasConnector::gemv(const char trans, const int m, const int n,
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