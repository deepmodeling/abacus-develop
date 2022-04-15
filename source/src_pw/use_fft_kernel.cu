#include "use_fft.h"
#include "global.h"
#include "cufft.h"
using namespace CudaCheck;

template<class T2>
__global__ void kernel_set(int size, T2 *dst, const T2 *src, const int *index_list)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int p = index_list[idx];
    if(idx < size)
    {
        dst[p].x = src[idx].x;
        dst[p].y = src[idx].y;
    }
}

template<class T, class T2>
__global__ void kernel_roundtrip(int size, T2 *dst, const T *src)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x *= src[idx];
        dst[idx].y *= src[idx];
    }
}

template<class T, class T2>
__global__ void kernel_normalization(int size, T2 *data, T norm)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        data[idx].x /= norm;
        data[idx].y /= norm;
    }
}

__global__ void kernel_reorder(double2 *dst, double2 *src, int size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < size)
    {
        dst[idx].x = src[idx].x;
        dst[idx].y = src[idx].y;
    }
}

void RoundTrip_kernel(const float2 *psi, const float *vr, const int *fft_index, float2 *psic)
{
    // (1) set value
    int thread = 512;
    int block = (GlobalC::wf.npw + thread - 1) / thread;
    int block2 = (GlobalC::pw.nrxx + thread - 1) / thread;
    kernel_set<float2><<<block, thread>>>(GlobalC::wf.npw, psic, psi, fft_index);

    CHECK_CUFFT(cufftExecC2C(GlobalC::UFFT.fft_handle, psic, psic, CUFFT_INVERSE));
    cudaDeviceSynchronize();

    kernel_roundtrip<float, float2><<<block2, thread>>>(GlobalC::pw.nrxx, psic, vr);

    CHECK_CUFFT(cufftExecC2C(GlobalC::UFFT.fft_handle, psic, psic, CUFFT_FORWARD));
    cudaDeviceSynchronize();

    int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
    kernel_normalization<float, float2><<<block3, thread>>>(GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    return;
}

void RoundTrip_kernel(const double2 *psi, const double *vr, const int *fft_index, double2 *psic)
{
    // (1) set value
    int thread = 512;
    int block = (GlobalC::wf.npw + thread - 1) / thread;
    int block2 = (GlobalC::pw.nrxx + thread - 1) / thread;
    kernel_set<double2><<<block, thread>>>(GlobalC::wf.npw, psic, psi, fft_index);

    CHECK_CUFFT(cufftExecZ2Z(GlobalC::UFFT.fft_handle, psic, psic, CUFFT_INVERSE));
    cudaDeviceSynchronize();

    kernel_roundtrip<double, double2><<<block2, thread>>>(GlobalC::pw.nrxx, psic, vr);

    CHECK_CUFFT(cufftExecZ2Z(GlobalC::UFFT.fft_handle, psic, psic, CUFFT_FORWARD));
    cudaDeviceSynchronize();

    int block3 = (GlobalC::pw.nrxx + thread - 1) / thread;
    kernel_normalization<double, double2><<<block3, thread>>>(GlobalC::pw.nrxx, psic, (double)(GlobalC::pw.nrxx));

    return;
}