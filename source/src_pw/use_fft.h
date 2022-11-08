#ifndef USE_FFT_H
#define USE_FFT_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "../module_pw/pw_basis.h"

#ifdef __CUDA
#include "cufft.h"
#include "use_fft_kernel.h"
#endif

#ifdef __ROCM
#include "hipfft.h"
#include "use_fft_kernel.h"
#endif 

class Use_FFT
{
	public:

	Use_FFT(){};
	~Use_FFT(){};

#ifdef __CUDA
	double2 *d_porter;
	cufftHandle fft_handle;
	void RoundTrip(const float2 *psi, const float *vr, const int *fft_index, float2 *psic)
	{
		RoundTrip_kernel(psi, vr, fft_index, psic);
	}
	void RoundTrip(const double2 *psi, const double *vr, const int *fft_index, double2 *psic)
	{
		RoundTrip_kernel(psi, vr, fft_index, psic);
	}
#endif

#ifdef __ROCM
	hipfftHandle fft_handle;
	void RoundTrip(const hipblasComplex *psi, const float *vr, const int *fft_index, hipblasComplex *psic)
	{
		RoundTrip_kernel(psi, vr, fft_index, psic);
	}
	void RoundTrip(const hipblasDoubleComplex *psi, const double *vr, const int *fft_index, hipblasDoubleComplex *psic)
	{
		RoundTrip_kernel(psi, vr, fft_index, psic);
	}
#endif

private:

};

#endif
