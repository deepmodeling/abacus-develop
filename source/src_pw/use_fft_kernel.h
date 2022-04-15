#ifndef USE_FFT_KERNEL_H
#define USE_FFT_KERNEL_H

#ifdef __CUDA
#include "cuda_runtime.h"
#include "cufft.h"
#include "device_launch_parameters.h"

void RoundTrip_kernel(const float2 *psi, const float *vr, const int *fft_index, float2 *psic);
void RoundTrip_kernel(const double2 *psi, const double *vr, const int *fft_index, double2 *psic);

#endif

#ifdef __ROCM
#include "hip/hip_runtime.h"
#include "hipblas.h"
#include "hipfft.h"

void RoundTrip_kernel(const hipblasComplex *psi, const float *vr, const int *fft_index, hipblasComplex *psic);
void RoundTrip_kernel(const hipblasDoubleComplex *psi, const double *vr, const int *fft_index, hipblasDoubleComplex *psic);

#endif

#endif