#ifndef DEVICE_CHECK_H
#define DEVICE_CHECK_H

#include <stdio.h>

#ifdef __CUDA
#include "cublas_v2.h"
#include "cufft.h"
#include "source_base/module_device/cuda_compat.h"

static const char* _cublasGetErrorString(cublasStatus_t error)
{
    switch (error)
    {
    case CUBLAS_STATUS_SUCCESS:
        return "CUBLAS_STATUS_SUCCESS";
    case CUBLAS_STATUS_NOT_INITIALIZED:
        return "CUBLAS_STATUS_NOT_INITIALIZED";
    case CUBLAS_STATUS_ALLOC_FAILED:
        return "CUBLAS_STATUS_ALLOC_FAILED";
    case CUBLAS_STATUS_INVALID_VALUE:
        return "CUBLAS_STATUS_INVALID_VALUE";
    case CUBLAS_STATUS_ARCH_MISMATCH:
        return "CUBLAS_STATUS_ARCH_MISMATCH";
    case CUBLAS_STATUS_MAPPING_ERROR:
        return "CUBLAS_STATUS_MAPPING_ERROR";
    case CUBLAS_STATUS_EXECUTION_FAILED:
        return "CUBLAS_STATUS_EXECUTION_FAILED";
    case CUBLAS_STATUS_INTERNAL_ERROR:
        return "CUBLAS_STATUS_INTERNAL_ERROR";
    }
    return "<unknown>";
}

#define CHECK_CUDA(func)                                                                                               \
    {                                                                                                                  \
        cudaError_t status = (func);                                                                                   \
        if (status != cudaSuccess)                                                                                     \
        {                                                                                                              \
            printf("In File %s : CUDA API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,                \
                   cudaGetErrorString(status), status);                                                                \
        }                                                                                                              \
    }

#define CHECK_CUBLAS(func)                                                                                             \
    {                                                                                                                  \
        cublasStatus_t status = (func);                                                                                \
        if (status != CUBLAS_STATUS_SUCCESS)                                                                           \
        {                                                                                                              \
            printf("In File %s : CUBLAS API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,              \
                   _cublasGetErrorString(status), status);                                                             \
        }                                                                                                              \
    }

#define CHECK_CUSOLVER(func)                                                                                           \
    {                                                                                                                  \
        cusolverStatus_t status = (func);                                                                              \
        if (status != CUSOLVER_STATUS_SUCCESS)                                                                         \
        {                                                                                                              \
            printf("In File %s : CUSOLVER API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,            \
                   _cusolverGetErrorString(status), status);                                                           \
        }                                                                                                              \
    }

#define CHECK_CUFFT(func)                                                                                              \
    {                                                                                                                  \
        cufftResult_t status = (func);                                                                                 \
        if (status != CUFFT_SUCCESS)                                                                                   \
        {                                                                                                              \
            printf("In File %s : CUFFT API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,               \
                   ModuleBase::cuda_compat::cufftGetErrorStringCompat(status), status);                                                              \
        }                                                                                                              \
    }
#endif // __CUDA

#ifdef __ROCM
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <hipfft/hipfft.h>

static const char* _hipblasGetErrorString(hipblasStatus_t error)
{
    switch (error)
    {
    case HIPBLAS_STATUS_SUCCESS:
        return "HIPBLAS_STATUS_SUCCESS";
    case HIPBLAS_STATUS_NOT_INITIALIZED:
        return "HIPBLAS_STATUS_NOT_INITIALIZED";
    case HIPBLAS_STATUS_ALLOC_FAILED:
        return "HIPBLAS_STATUS_ALLOC_FAILED";
    case HIPBLAS_STATUS_INVALID_VALUE:
        return "HIPBLAS_STATUS_INVALID_VALUE";
    case HIPBLAS_STATUS_ARCH_MISMATCH:
        return "HIPBLAS_STATUS_ARCH_MISMATCH";
    case HIPBLAS_STATUS_MAPPING_ERROR:
        return "HIPBLAS_STATUS_MAPPING_ERROR";
    case HIPBLAS_STATUS_EXECUTION_FAILED:
        return "HIPBLAS_STATUS_EXECUTION_FAILED";
    case HIPBLAS_STATUS_INTERNAL_ERROR:
        return "HIPBLAS_STATUS_INTERNAL_ERROR";
    case HIPBLAS_STATUS_NOT_SUPPORTED:
        return "HIPBLAS_STATUS_NOT_SUPPORTED";
    case HIPBLAS_STATUS_HANDLE_IS_NULLPTR:
        return "HIPBLAS_STATUS_HANDLE_IS_NULLPTR";
    default:
        return "<unknown>";
    }
    return "<unknown>";
}

static const char* _hipfftGetErrorString(hipfftResult_t error)
{
    switch (error)
    {
    case HIPFFT_SUCCESS:
        return "HIPFFT_SUCCESS";
    case HIPFFT_INVALID_PLAN:
        return "HIPFFT_INVALID_PLAN";
    case HIPFFT_ALLOC_FAILED:
        return "HIPFFT_ALLOC_FAILED";
    case HIPFFT_INVALID_TYPE:
        return "HIPFFT_INVALID_TYPE";
    case HIPFFT_INVALID_VALUE:
        return "HIPFFT_INVALID_VALUE";
    case HIPFFT_INTERNAL_ERROR:
        return "HIPFFT_INTERNAL_ERROR";
    case HIPFFT_EXEC_FAILED:
        return "HIPFFT_EXEC_FAILED";
    case HIPFFT_SETUP_FAILED:
        return "HIPFFT_SETUP_FAILED";
    case HIPFFT_INVALID_SIZE:
        return "HIPFFT_INVALID_SIZE";
    case HIPFFT_UNALIGNED_DATA:
        return "HIPFFT_UNALIGNED_DATA";
    case HIPFFT_INCOMPLETE_PARAMETER_LIST:
        return "HIPFFT_INCOMPLETE_PARAMETER_LIST";
    case HIPFFT_INVALID_DEVICE:
        return "HIPFFT_INVALID_DEVICE";
    case HIPFFT_PARSE_ERROR:
        return "HIPFFT_PARSE_ERROR";
    case HIPFFT_NO_WORKSPACE:
        return "HIPFFT_NO_WORKSPACE";
    case HIPFFT_NOT_IMPLEMENTED:
        return "HIPFFT_NOT_IMPLEMENTED";
    case HIPFFT_NOT_SUPPORTED:
        return "HIPFFT_NOT_SUPPORTED";
    }
    return "<unknown>";
}

#define CHECK_CUDA(func)                                                                                               \
    {                                                                                                                  \
        hipError_t status = (func);                                                                                    \
        if (status != hipSuccess)                                                                                      \
        {                                                                                                              \
            printf("In File %s : HIP API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,                 \
                   hipGetErrorString(status), status);                                                                 \
        }                                                                                                              \
    }

#define CHECK_CUBLAS(func)                                                                                             \
    {                                                                                                                  \
        hipblasStatus_t status = (func);                                                                               \
        if (status != HIPBLAS_STATUS_SUCCESS)                                                                          \
        {                                                                                                              \
            printf("In File %s : HIPBLAS API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,             \
                   _hipblasGetErrorString(status), status);                                                            \
        }                                                                                                              \
    }

#define CHECK_CUFFT(func)                                                                                              \
    {                                                                                                                  \
        hipfftResult_t status = (func);                                                                                \
        if (status != HIPFFT_SUCCESS)                                                                                  \
        {                                                                                                              \
            printf("In File %s : HIPFFT API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,              \
                   _hipfftGetErrorString(status), status);                                                             \
        }                                                                                                              \
    }
#endif // __ROCM

#endif // DEVICE_CHECK_H