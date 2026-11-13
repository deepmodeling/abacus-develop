#ifndef MODULEPW_FUSED_FFT_CHECK_CUH
#define MODULEPW_FUSED_FFT_CHECK_CUH

#include <cuda_runtime.h>

#include <cstdio>
#include <cstdlib>

inline void cudaErrcheck(cudaError_t status)
{
    if (status != cudaSuccess)
    {
        std::fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(status));
        std::exit(EXIT_FAILURE);
    }
}

inline void cudaCheckOnDebug()
{
    cudaErrcheck(cudaGetLastError());
}

#endif
