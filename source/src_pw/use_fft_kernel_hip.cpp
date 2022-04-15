#include "global.h"
#include "hip/hip_runtime.h"
#include "hipfft.h"
#include "use_fft.h"
using namespace HipCheck;

template <class T2> __global__ void kernel_set(int size, T2 *dst, const T2 *src, const int *index_list)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int p = index_list[idx];
	if (idx < size)
	{
		dst[p].x = src[idx].x;
		dst[p].y = src[idx].y;
	}
}

template <class T, class T2> __global__ void kernel_roundtrip(int size, T2 *dst, const T *src)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		dst[idx].x *= src[idx];
		dst[idx].y *= src[idx];
	}
}

template <class T, class T2> __global__ void kernel_normalization(int size, T2 *data, T norm)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		data[idx].x /= norm;
		data[idx].y /= norm;
	}
}

__global__ void kernel_reorder(double2 *dst, double2 *src, int size)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	if (idx < size)
	{
		dst[idx].x = src[idx].x;
		dst[idx].y = src[idx].y;
	}
}

void RoundTrip_kernel(const hipblasComplex *psi, const float *vr, const int *fft_index, hipblasComplex *psic)
{
	// (1) set value
	int thread = 512;
	int block = (GlobalC::wf.npw + thread - 1) / thread;
	int block2 = (GlobalC::pw.nrxx + thread - 1) / thread;
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_set<float2>),
					   dim3(block),
					   dim3(thread),
					   0,
					   0,
					   GlobalC::wf.npw,
					   reinterpret_cast<float2 *>(psic),
					   reinterpret_cast<const float2 *>(psi),
					   fft_index);

	CHECK_CUFFT(hipfftExecC2C(GlobalC::UFFT.fft_handle,
							  reinterpret_cast<hipfftComplex *>(psic),
							  reinterpret_cast<hipfftComplex *>(psic),
							  HIPFFT_BACKWARD));
	hipDeviceSynchronize();
	// CHECK_CUFFT(hipfftDestroy(cufftplan_gpu));

	// int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
	// hipLaunchKernelGGL(kernel_normalization, dim3(block3), dim3(thread), 0, 0, GlobalC::pw.nrxx, psic,
	// (double)(GlobalC::pw.nrxx));

	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_roundtrip<float, float2>),
					   dim3(block2),
					   dim3(thread),
					   0,
					   0,
					   GlobalC::pw.nrxx,
					   reinterpret_cast<float2 *>(psic),
					   vr);

	// hipfftHandle cufftplan_gpu2;
	// CHECK_CUFFT(hipfftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_Z2Z));
	CHECK_CUFFT(hipfftExecC2C(GlobalC::UFFT.fft_handle,
							  reinterpret_cast<hipfftComplex *>(psic),
							  reinterpret_cast<hipfftComplex *>(psic),
							  HIPFFT_FORWARD));
	hipDeviceSynchronize();
	// CHECK_CUFFT(hipfftDestroy(cufftplan_gpu));

	// Reorder_psi_minus(psic, ordered_psi);

	int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_normalization<float, float2>),
					   dim3(block3),
					   dim3(thread),
					   0,
					   0,
					   GlobalC::pw.nrxx,
					   reinterpret_cast<float2 *>(psic),
					   (double)(GlobalC::pw.nrxx));

	// CHECK_CUDA(hipFree(ordered_psi));

	return;
}

void RoundTrip_kernel(const hipblasDoubleComplex *psi,
					  const double *vr,
					  const int *fft_index,
					  hipblasDoubleComplex *psic)
{
	// (1) set value
	int thread = 512;
	int block = (GlobalC::wf.npw + thread - 1) / thread;
	int block2 = (GlobalC::pw.nrxx + thread - 1) / thread;
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_set<double2>),
					   dim3(block),
					   dim3(thread),
					   0,
					   0,
					   GlobalC::wf.npw,
					   reinterpret_cast<double2 *>(psic),
					   reinterpret_cast<const double2 *>(psi),
					   fft_index);
	CHECK_CUFFT(hipfftExecZ2Z(GlobalC::UFFT.fft_handle,
							  (hipfftDoubleComplex *)(psic),
							  (hipfftDoubleComplex *)(psic),
							  HIPFFT_BACKWARD));
	hipDeviceSynchronize();

	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_roundtrip<double, double2>),
					   dim3(block2),
					   dim3(thread),
					   0,
					   0,
					   GlobalC::pw.nrxx,
					   reinterpret_cast<double2 *>(psic),
					   vr);

	// hipfftHandle cufftplan_gpu2;
	// CHECK_CUFFT(hipfftPlan3d(&cufftplan_gpu, GlobalC::pw.nx, GlobalC::pw.ny, GlobalC::pw.nz, HIPFFT_Z2Z));
	CHECK_CUFFT(hipfftExecZ2Z(GlobalC::UFFT.fft_handle,
							  reinterpret_cast<hipfftDoubleComplex *>(psic),
							  reinterpret_cast<hipfftDoubleComplex *>(psic),
							  HIPFFT_FORWARD));
	hipDeviceSynchronize();
	// CHECK_CUFFT(hipfftDestroy(cufftplan_gpu));

	// Reorder_psi_minus(psic, ordered_psi);

	int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
	hipLaunchKernelGGL(HIP_KERNEL_NAME(kernel_normalization<double, double2>),
					   dim3(block3),
					   dim3(thread),
					   0,
					   0,
					   GlobalC::pw.nrxx,
					   reinterpret_cast<double2 *>(psic),
					   (double)(GlobalC::pw.nrxx));

	// CHECK_CUDA(hipFree(ordered_psi));

	return;
}