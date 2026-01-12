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


void BlasConnector::axpy( const int n, const float alpha, const float *X, const int incX, float *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		saxpy_(&n, &alpha, X, &incX, Y, &incY);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasErrcheck(cublasSaxpy(BlasUtils::cublas_handle, n, &alpha, X, incX, Y, incY));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::axpy( const int n, const double alpha, const double *X, const int incX, double *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		daxpy_(&n, &alpha, X, &incX, Y, &incY);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasErrcheck(cublasDaxpy(BlasUtils::cublas_handle, n, &alpha, X, incX, Y, incY));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::axpy( const int n, const std::complex<float> alpha, const std::complex<float> *X, const int incX, std::complex<float> *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		caxpy_(&n, &alpha, X, &incX, Y, &incY);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasErrcheck(cublasCaxpy(BlasUtils::cublas_handle, n, (float2*)&alpha, (float2*)X, incX, (float2*)Y, incY));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::axpy( const int n, const std::complex<double> alpha, const std::complex<double> *X, const int incX, std::complex<double> *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zaxpy_(&n, &alpha, X, &incX, Y, &incY);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasErrcheck(cublasZaxpy(BlasUtils::cublas_handle, n, (double2*)&alpha, (double2*)X, incX, (double2*)Y, incY));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}


// x=a*x
void BlasConnector::scal( const int n,  const float alpha, float *X, const int incX, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		sscal_(&n, &alpha, X, &incX);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasErrcheck(cublasSscal(BlasUtils::cublas_handle, n, &alpha, X, incX));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::scal( const int n, const double alpha, double *X, const int incX, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		dscal_(&n, &alpha, X, &incX);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasErrcheck(cublasDscal(BlasUtils::cublas_handle, n, &alpha, X, incX));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::scal( const int n, const std::complex<float> alpha, std::complex<float> *X, const int incX, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		cscal_(&n, &alpha, X, &incX);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasErrcheck(cublasCscal(BlasUtils::cublas_handle, n, (float2*)&alpha, (float2*)X, incX));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::scal( const int n, const std::complex<double> alpha, std::complex<double> *X, const int incX, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zscal_(&n, &alpha, X, &incX);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		cublasErrcheck(cublasZscal(BlasUtils::cublas_handle, n, (double2*)&alpha, (double2*)X, incX));
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}


// d=x*y
float BlasConnector::dot( const int n, const float*const X, const int incX, const float*const Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		return sdot_(&n, X, &incX, Y, &incY);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice){
		float result = 0.0;
		cublasErrcheck(cublasSdot(BlasUtils::cublas_handle, n, X, incX, Y, incY, &result));
		return result;
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

double BlasConnector::dot( const int n, const double*const X, const int incX, const double*const Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		return ddot_(&n, X, &incX, Y, &incY);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice){
		double result = 0.0;
		cublasErrcheck(cublasDdot(BlasUtils::cublas_handle, n, X, incX, Y, incY, &result));
		return result;
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

// d=x*y
float BlasConnector::dotu(const int n, const float*const X, const int incX, const float*const Y, const int incY, base_device::AbacusDevice_t device_type)
{
	return BlasConnector::dot(n, X, incX, Y, incY, device_type);
}

double BlasConnector::dotu(const int n, const double*const X, const int incX, const double*const Y, const int incY, base_device::AbacusDevice_t device_type)
{
	return BlasConnector::dot(n, X, incX, Y, incY, device_type);
}

std::complex<float> BlasConnector::dotu(const int n, const std::complex<float>*const X, const int incX, const std::complex<float>*const Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		const int incX2 = 2 * incX;
		const int incY2 = 2 * incY;
		const float*const x = reinterpret_cast<const float*const>(X);
		const float*const y = reinterpret_cast<const float*const>(Y);
		//Re(result)=Re(x)*Re(y)-Im(x)*Im(y)
		//Im(result)=Re(x)*Im(y)+Im(x)*Re(y)
		return std::complex<float>(
			BlasConnector::dot(n, x, incX2, y,   incY2, device_type) - dot(n, x+1, incX2, y+1, incY2, device_type),
			BlasConnector::dot(n, x, incX2, y+1, incY2, device_type) + dot(n, x+1, incX2, y,   incY2, device_type));
	}
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

std::complex<double> BlasConnector::dotu(const int n, const std::complex<double>*const X, const int incX, const std::complex<double>*const Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		const int incX2 = 2 * incX;
		const int incY2 = 2 * incY;
		const double*const x = reinterpret_cast<const double*const>(X);
		const double*const y = reinterpret_cast<const double*const>(Y);
		//Re(result)=Re(x)*Re(y)-Im(x)*Im(y)
		//Im(result)=Re(x)*Im(y)+Im(x)*Re(y)
		return std::complex<double>(
			BlasConnector::dot(n, x, incX2, y,   incY2, device_type) - dot(n, x+1, incX2, y+1, incY2, device_type),
			BlasConnector::dot(n, x, incX2, y+1, incY2, device_type) + dot(n, x+1, incX2, y,   incY2, device_type));
	}
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

// d = x.conj() * Vy
float BlasConnector::dotc(const int n, const float*const X, const int incX, const float*const Y, const int incY, base_device::AbacusDevice_t device_type)
{
	return BlasConnector::dot(n, X, incX, Y, incY, device_type);
}

double BlasConnector::dotc(const int n, const double*const X, const int incX, const double*const Y, const int incY, base_device::AbacusDevice_t device_type)
{
	return BlasConnector::dot(n, X, incX, Y, incY, device_type);
}

std::complex<float> BlasConnector::dotc(const int n, const std::complex<float>*const X, const int incX, const std::complex<float>*const Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		const int incX2 = 2 * incX;
		const int incY2 = 2 * incY;
		const float*const x = reinterpret_cast<const float*const>(X);
		const float*const y = reinterpret_cast<const float*const>(Y);
		// Re(result)=Re(X)*Re(Y)+Im(X)*Im(Y)
		// Im(result)=Re(X)*Im(Y)-Im(X)*Re(Y)
		return std::complex<float>(
			BlasConnector::dot(n, x, incX2, y,   incY2, device_type) + dot(n, x+1, incX2, y+1, incY2, device_type),
			BlasConnector::dot(n, x, incX2, y+1, incY2, device_type) - dot(n, x+1, incX2, y,   incY2, device_type));
	}
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

std::complex<double> BlasConnector::dotc(const int n, const std::complex<double>*const X, const int incX, const std::complex<double>*const Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		const int incX2 = 2 * incX;
		const int incY2 = 2 * incY;
		const double*const x = reinterpret_cast<const double*const>(X);
		const double*const y = reinterpret_cast<const double*const>(Y);
		// Re(result)=Re(X)*Re(Y)+Im(X)*Im(Y)
		// Im(result)=Re(X)*Im(Y)-Im(X)*Re(Y)
		return std::complex<double>(
			BlasConnector::dot(n, x, incX2, y,   incY2, device_type) + dot(n, x+1, incX2, y+1, incY2, device_type),
			BlasConnector::dot(n, x, incX2, y+1, incY2, device_type) - dot(n, x+1, incX2, y,   incY2, device_type));
	}
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

// out = ||x||_2
float BlasConnector::nrm2( const int n, const float *X, const int incX, base_device::AbacusDevice_t device_type )
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		return snrm2_( &n, X, &incX );
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice){
		float result = 0.0;
		cublasErrcheck(cublasSnrm2(BlasUtils::cublas_handle, n, X, incX, &result));
		return result;
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}


double BlasConnector::nrm2( const int n, const double *X, const int incX, base_device::AbacusDevice_t device_type )
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		return dnrm2_( &n, X, &incX );
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice){
		double result = 0.0;
		cublasErrcheck(cublasDnrm2(BlasUtils::cublas_handle, n, X, incX, &result));
		return result;
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}


double BlasConnector::nrm2( const int n, const std::complex<double> *X, const int incX, base_device::AbacusDevice_t device_type )
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		return dznrm2_( &n, X, &incX );
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice){
		double result = 0.0;
		cublasErrcheck(cublasDznrm2(BlasUtils::cublas_handle, n, (double2*)X, incX, &result));
		return result;
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

// copies a into b
void BlasConnector::copy(const int n, const float *a, const int incx, float *b, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		scopy_(&n, a, &incx, b, &incy);
	}
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::copy(const int n, const double *a, const int incx, double *b, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		dcopy_(&n, a, &incx, b, &incy);
	}
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::copy(const int n, const std::complex<float> *a, const int incx, std::complex<float> *b, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		ccopy_(&n, a, &incx, b, &incy);
	}
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void BlasConnector::copy(const int n, const std::complex<double> *a, const int incx, std::complex<double> *b, const int incy, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		zcopy_(&n, a, &incx, b, &incy);
	}
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}


template <typename T>
void vector_mul_vector(const int& dim, T* result, const T* vector1, const T* vector2, base_device::AbacusDevice_t device_type){
	using Real = typename GetTypeReal<T>::type;
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(Real))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * vector2[i];
        }
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		ModuleBase::vector_mul_vector_op<T, base_device::DEVICE_GPU>()(dim, result, vector1, vector2);
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}


template <typename T>
void vector_div_vector(const int& dim, T* result, const T* vector1, const T* vector2, base_device::AbacusDevice_t device_type){
	using Real = typename GetTypeReal<T>::type;
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096 / sizeof(Real))
#endif
		for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] / vector2[i];
        }
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		ModuleBase::vector_div_vector_op<T, base_device::DEVICE_GPU>()(dim, result, vector1, vector2);
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void vector_add_vector(const int& dim, float *result, const float *vector1, const float constant1, const float *vector2, const float constant2, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::CpuDevice){
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 8192 / sizeof(float))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * constant1 + vector2[i] * constant2;
        }
	}
#ifdef __CUDA
	else if (device_type == base_device::GpuDevice) {
		ModuleBase::vector_add_vector_op<float, base_device::DEVICE_GPU>()(dim, result, vector1, constant1, vector2, constant2);
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void vector_add_vector(const int& dim, double *result, const double *vector1, const double constant1, const double *vector2, const double constant2, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::CpuDevice){
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 8192 / sizeof(double))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * constant1 + vector2[i] * constant2;
        }
	}
#ifdef __CUDA
	else if (device_type == base_device::GpuDevice) {
		ModuleBase::vector_add_vector_op<double, base_device::DEVICE_GPU>()(dim, result, vector1, constant1, vector2, constant2);
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void vector_add_vector(const int& dim, std::complex<float> *result, const std::complex<float> *vector1, const float constant1, const std::complex<float> *vector2, const float constant2, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::CpuDevice){
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 8192 / sizeof(std::complex<float>))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * constant1 + vector2[i] * constant2;
        }
	}
#ifdef __CUDA
	else if (device_type == base_device::GpuDevice) {
		ModuleBase::vector_add_vector_op<std::complex<float>, base_device::DEVICE_GPU>()(dim, result, vector1, constant1, vector2, constant2);
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}

void vector_add_vector(const int& dim, std::complex<double> *result, const std::complex<double> *vector1, const double constant1, const std::complex<double> *vector2, const double constant2, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::CpuDevice){
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 8192 / sizeof(std::complex<double>))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * constant1 + vector2[i] * constant2;
        }
	}
#ifdef __CUDA
	else if (device_type == base_device::GpuDevice) {
		ModuleBase::vector_add_vector_op<std::complex<double>, base_device::DEVICE_GPU>()(dim, result, vector1, constant1, vector2, constant2);
	}
#endif
	else {
		throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__) + " line " + std::to_string(__LINE__));
	}
}