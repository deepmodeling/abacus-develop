#ifndef GPU_RUNTIME_H_
#define GPU_RUNTIME_H_

/**
 * @file gpu_runtime.h
 * @brief Unified CUDA/ROCm API macros for portable GPU code.
 *
 * This header provides macro abstraction for CUDA/ROCm APIs, allowing
 * a single implementation to work with both CUDA and ROCm backends.
 *
 * Usage:
 *   #include "gpu_runtime.h"
 *   gpuError_t err = gpuGetDeviceCount(&count);
 *   if (err != gpuSuccess) { ... }
 */

#if defined(__CUDA)

#include <cuda_runtime.h>

// Error handling
#define gpuError_t                      cudaError_t
#define gpuSuccess                      cudaSuccess
#define gpuGetErrorString               cudaGetErrorString

// Device management
#define gpuGetDeviceCount               cudaGetDeviceCount
#define gpuGetDevice                    cudaGetDevice
#define gpuSetDevice                    cudaSetDevice
#define gpuGetDeviceProperties          cudaGetDeviceProperties
#define gpuDeviceProp_t                 cudaDeviceProp

// Version info
#define gpuDriverGetVersion             cudaDriverGetVersion
#define gpuRuntimeGetVersion            cudaRuntimeGetVersion

// Peer access
#define gpuDeviceCanAccessPeer          cudaDeviceCanAccessPeer

// Error check macro
#define gpuErrcheck                     CHECK_CUDA

#elif defined(__ROCM)

#include <hip/hip_runtime.h>

// Error handling
#define gpuError_t                      hipError_t
#define gpuSuccess                      hipSuccess
#define gpuGetErrorString               hipGetErrorString

// Device management
#define gpuGetDeviceCount               hipGetDeviceCount
#define gpuGetDevice                    hipGetDevice
#define gpuSetDevice                    hipSetDevice
#define gpuGetDeviceProperties          hipGetDeviceProperties
#define gpuDeviceProp_t                 hipDeviceProp_t

// Version info
#define gpuDriverGetVersion             hipDriverGetVersion
#define gpuRuntimeGetVersion            hipRuntimeGetVersion

// Peer access
#define gpuDeviceCanAccessPeer          hipDeviceCanAccessPeer

// Error check macro
#define gpuErrcheck                     CHECK_CUDA

#endif // __CUDA / __ROCM

#endif // GPU_RUNTIME_H_
