#include "blas_connector.h"
#include "../macros.h"
#include "source_base/kernels/math_kernel_op.h"

#ifdef __DSP
#include "source_base/kernels/dsp/dsp_connector.h"
#include "source_base/global_variable.h"
#endif

#ifdef __CUDA
#include <base/macros/macros.h>
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include "source_base/module_device/memory_op.h"
#endif


void BlasConnector::axpy( const int n, const float alpha, const float *X, const int incX, float *Y, const int incY, base_device::AbacusDevice_t device_type)
{
	if (device_type == base_device::AbacusDevice_t::CpuDevice) {
		saxpy_(&n, &alpha, X, &incX, Y, &incY);
	}
#ifdef __CUDA
	else if (device_type == base_device::AbacusDevice_t::GpuDevice) {
		CHECK_CUBLAS(cublasSaxpy(BlasUtils::cublas_handle, n, &alpha, X, incX, Y, incY));
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
		CHECK_CUBLAS(cublasDaxpy(BlasUtils::cublas_handle, n, &alpha, X, incX, Y, incY));
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
		CHECK_CUBLAS(cublasCaxpy(BlasUtils::cublas_handle, n, (float2*)&alpha, (float2*)X, incX, (float2*)Y, incY));
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
		CHECK_CUBLAS(cublasZaxpy(BlasUtils::cublas_handle, n, (double2*)&alpha, (double2*)X, incX, (double2*)Y, incY));
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
		CHECK_CUBLAS(cublasSscal(BlasUtils::cublas_handle, n, &alpha, X, incX));
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
		CHECK_CUBLAS(cublasDscal(BlasUtils::cublas_handle, n, &alpha, X, incX));
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
		CHECK_CUBLAS(cublasCscal(BlasUtils::cublas_handle, n, (float2*)&alpha, (float2*)X, incX));
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
		CHECK_CUBLAS(cublasZscal(BlasUtils::cublas_handle, n, (double2*)&alpha, (double2*)X, incX));
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
		CHECK_CUBLAS(cublasSdot(BlasUtils::cublas_handle, n, X, incX, Y, incY, &result));
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
		CHECK_CUBLAS(cublasDdot(BlasUtils::cublas_handle, n, X, incX, Y, incY, &result));
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
		CHECK_CUBLAS(cublasSnrm2(BlasUtils::cublas_handle, n, X, incX, &result));
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
		CHECK_CUBLAS(cublasDnrm2(BlasUtils::cublas_handle, n, X, incX, &result));
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
		CHECK_CUBLAS(cublasDznrm2(BlasUtils::cublas_handle, n, (double2*)X, incX, &result));
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


template <typename T, typename Operand>
void BlasConnector::vector_mul_vector(const int& dim,
                                      T* result,
                                      const T* vector1,
                                      const Operand* vector2,
                                      base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        ModuleBase::vector_mul_vector_op<T, base_device::DEVICE_CPU, Operand>()(
            dim, result, vector1, vector2, false);
    }
#if defined(__CUDA) || defined(__ROCM)
    else if (device_type == base_device::AbacusDevice_t::GpuDevice)
    {
        ModuleBase::vector_mul_vector_op<T, base_device::DEVICE_GPU, Operand>()(
            dim, result, vector1, vector2, false);
    }
#endif
    else
    {
        throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__)
                                    + " line " + std::to_string(__LINE__));
    }
}

template <typename T, typename Operand>
void BlasConnector::vector_div_vector(const int& dim,
                                      T* result,
                                      const T* vector1,
                                      const Operand* vector2,
                                      base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        ModuleBase::vector_div_vector_op<T, base_device::DEVICE_CPU, Operand>()(dim, result, vector1, vector2);
    }
#if defined(__CUDA) || defined(__ROCM)
    else if (device_type == base_device::AbacusDevice_t::GpuDevice)
    {
        ModuleBase::vector_div_vector_op<T, base_device::DEVICE_GPU, Operand>()(dim, result, vector1, vector2);
    }
#endif
    else
    {
        throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__)
                                    + " line " + std::to_string(__LINE__));
    }
}

template <typename T, typename Scalar>
void BlasConnector::vector_add_vector(const int& dim,
                                      T* result,
                                      const T* vector1,
                                      const Scalar constant1,
                                      const T* vector2,
                                      const Scalar constant2,
                                      base_device::AbacusDevice_t device_type)
{
    if (device_type == base_device::AbacusDevice_t::CpuDevice)
    {
        ModuleBase::vector_add_vector_op<T, base_device::DEVICE_CPU, Scalar>()(
            dim, result, vector1, constant1, vector2, constant2);
    }
#if defined(__CUDA) || defined(__ROCM)
    else if (device_type == base_device::AbacusDevice_t::GpuDevice)
    {
        ModuleBase::vector_add_vector_op<T, base_device::DEVICE_GPU, Scalar>()(
            dim, result, vector1, constant1, vector2, constant2);
    }
#endif
    else
    {
        throw std::invalid_argument("device_type = " + std::to_string(device_type) + " in " + std::string(__FILE__)
                                    + " line " + std::to_string(__LINE__));
    }
}

template void BlasConnector::vector_mul_vector<float, float>(
    const int&, float*, const float*, const float*, base_device::AbacusDevice_t);
template void BlasConnector::vector_mul_vector<double, double>(
    const int&, double*, const double*, const double*, base_device::AbacusDevice_t);
template void BlasConnector::vector_mul_vector<std::complex<float>, float>(
    const int&, std::complex<float>*, const std::complex<float>*, const float*, base_device::AbacusDevice_t);
template void BlasConnector::vector_mul_vector<std::complex<double>, double>(
    const int&, std::complex<double>*, const std::complex<double>*, const double*, base_device::AbacusDevice_t);
template void BlasConnector::vector_mul_vector<std::complex<float>, std::complex<float>>(
    const int&,
    std::complex<float>*,
    const std::complex<float>*,
    const std::complex<float>*,
    base_device::AbacusDevice_t);
template void BlasConnector::vector_mul_vector<std::complex<double>, std::complex<double>>(
    const int&,
    std::complex<double>*,
    const std::complex<double>*,
    const std::complex<double>*,
    base_device::AbacusDevice_t);

template void BlasConnector::vector_div_vector<float, float>(
    const int&, float*, const float*, const float*, base_device::AbacusDevice_t);
template void BlasConnector::vector_div_vector<double, double>(
    const int&, double*, const double*, const double*, base_device::AbacusDevice_t);
template void BlasConnector::vector_div_vector<std::complex<float>, float>(
    const int&, std::complex<float>*, const std::complex<float>*, const float*, base_device::AbacusDevice_t);
template void BlasConnector::vector_div_vector<std::complex<double>, double>(
    const int&, std::complex<double>*, const std::complex<double>*, const double*, base_device::AbacusDevice_t);
template void BlasConnector::vector_div_vector<std::complex<float>, std::complex<float>>(
    const int&,
    std::complex<float>*,
    const std::complex<float>*,
    const std::complex<float>*,
    base_device::AbacusDevice_t);
template void BlasConnector::vector_div_vector<std::complex<double>, std::complex<double>>(
    const int&,
    std::complex<double>*,
    const std::complex<double>*,
    const std::complex<double>*,
    base_device::AbacusDevice_t);

template void BlasConnector::vector_add_vector<float, float>(
    const int&, float*, const float*, const float, const float*, const float, base_device::AbacusDevice_t);
template void BlasConnector::vector_add_vector<double, double>(
    const int&, double*, const double*, const double, const double*, const double, base_device::AbacusDevice_t);
template void BlasConnector::vector_add_vector<std::complex<float>, float>(const int&,
                                                                           std::complex<float>*,
                                                                           const std::complex<float>*,
                                                                           const float,
                                                                           const std::complex<float>*,
                                                                           const float,
                                                                           base_device::AbacusDevice_t);
template void BlasConnector::vector_add_vector<std::complex<double>, double>(const int&,
                                                                             std::complex<double>*,
                                                                             const std::complex<double>*,
                                                                             const double,
                                                                             const std::complex<double>*,
                                                                             const double,
                                                                             base_device::AbacusDevice_t);
template void BlasConnector::vector_add_vector<std::complex<float>, std::complex<float>>(
    const int&,
    std::complex<float>*,
    const std::complex<float>*,
    const std::complex<float>,
    const std::complex<float>*,
    const std::complex<float>,
    base_device::AbacusDevice_t);
template void BlasConnector::vector_add_vector<std::complex<double>, std::complex<double>>(
    const int&,
    std::complex<double>*,
    const std::complex<double>*,
    const std::complex<double>,
    const std::complex<double>*,
    const std::complex<double>,
    base_device::AbacusDevice_t);
