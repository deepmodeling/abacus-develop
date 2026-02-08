#ifndef DEVICE_CHECK_H
#define DEVICE_CHECK_H

#include <cstdlib>
#include <cstdio>

#ifdef __CUDA
#include "cublas_v2.h"
#include "cufft.h"
#include "cusolverDn.h"
#include <cuda.h>

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
    case CUBLAS_STATUS_NOT_SUPPORTED:
        return "CUBLAS_STATUS_NOT_SUPPORTED";
    case CUBLAS_STATUS_LICENSE_ERROR:
        return "CUBLAS_STATUS_LICENSE_ERROR";
    default:
        return "<unknown>";
    }
}

static const char* _cusolverGetErrorString(cusolverStatus_t error)
{
    switch (error)
    {
    case CUSOLVER_STATUS_SUCCESS:
        return "CUSOLVER_STATUS_SUCCESS";
    case CUSOLVER_STATUS_NOT_INITIALIZED:
        return "CUSOLVER_STATUS_NOT_INITIALIZED";
    case CUSOLVER_STATUS_ALLOC_FAILED:
        return "CUSOLVER_STATUS_ALLOC_FAILED";
    case CUSOLVER_STATUS_INVALID_VALUE:
        return "CUSOLVER_STATUS_INVALID_VALUE";
    case CUSOLVER_STATUS_ARCH_MISMATCH:
        return "CUSOLVER_STATUS_ARCH_MISMATCH";
    case CUSOLVER_STATUS_MAPPING_ERROR:
        return "CUSOLVER_STATUS_MAPPING_ERROR";
    case CUSOLVER_STATUS_EXECUTION_FAILED:
        return "CUSOLVER_STATUS_EXECUTION_FAILED";
    case CUSOLVER_STATUS_INTERNAL_ERROR:
        return "CUSOLVER_STATUS_INTERNAL_ERROR";
    case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
        return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
    case CUSOLVER_STATUS_NOT_SUPPORTED:
        return "CUSOLVER_STATUS_NOT_SUPPORTED";
    case CUSOLVER_STATUS_ZERO_PIVOT:
        return "CUSOLVER_STATUS_ZERO_PIVOT";
    case CUSOLVER_STATUS_INVALID_LICENSE:
        return "CUSOLVER_STATUS_INVALID_LICENSE";
    case CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED:
        return "CUSOLVER_STATUS_IRS_PARAMS_NOT_INITIALIZED";
    case CUSOLVER_STATUS_IRS_PARAMS_INVALID:
        return "CUSOLVER_STATUS_IRS_PARAMS_INVALID";
    case CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC:
        return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_PREC";
    case CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE:
        return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_REFINE";
    case CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER:
        return "CUSOLVER_STATUS_IRS_PARAMS_INVALID_MAXITER";
    case CUSOLVER_STATUS_IRS_INTERNAL_ERROR:
        return "CUSOLVER_STATUS_IRS_INTERNAL_ERROR";
    case CUSOLVER_STATUS_IRS_NOT_SUPPORTED:
        return "CUSOLVER_STATUS_IRS_NOT_SUPPORTED";
    case CUSOLVER_STATUS_IRS_OUT_OF_RANGE:
        return "CUSOLVER_STATUS_IRS_OUT_OF_RANGE";
    case CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES:
        return "CUSOLVER_STATUS_IRS_NRHS_NOT_SUPPORTED_FOR_REFINE_GMRES";
    case CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED:
        return "CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED";
    case CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED:
        return "CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED";
    case CUSOLVER_STATUS_IRS_MATRIX_SINGULAR:
        return "CUSOLVER_STATUS_IRS_MATRIX_SINGULAR";
    case CUSOLVER_STATUS_INVALID_WORKSPACE:
        return "CUSOLVER_STATUS_INVALID_WORKSPACE";
    default:
        return "<unknown>";
    }
}

static const char* _cufftGetErrorString(cufftResult_t error)
{
    switch (error)
    {
    case CUFFT_SUCCESS:
        return "CUFFT_SUCCESS";
    case CUFFT_INVALID_PLAN:
        return "CUFFT_INVALID_PLAN";
    case CUFFT_ALLOC_FAILED:
        return "CUFFT_ALLOC_FAILED";
    case CUFFT_INVALID_TYPE:
        return "CUFFT_INVALID_TYPE";
    case CUFFT_INVALID_VALUE:
        return "CUFFT_INVALID_VALUE";
    case CUFFT_INTERNAL_ERROR:
        return "CUFFT_INTERNAL_ERROR";
    case CUFFT_EXEC_FAILED:
        return "CUFFT_EXEC_FAILED";
    case CUFFT_SETUP_FAILED:
        return "CUFFT_SETUP_FAILED";
    case CUFFT_INVALID_SIZE:
        return "CUFFT_INVALID_SIZE";
    case CUFFT_UNALIGNED_DATA:
        return "CUFFT_UNALIGNED_DATA";
    case CUFFT_INVALID_DEVICE:
        return "CUFFT_INVALID_DEVICE";
    case CUFFT_NO_WORKSPACE:
        return "CUFFT_NO_WORKSPACE";
    case CUFFT_NOT_IMPLEMENTED:
        return "CUFFT_NOT_IMPLEMENTED";
    case CUFFT_NOT_SUPPORTED:
        return "CUFFT_NOT_SUPPORTED";
#if defined(CUDA_VERSION) && CUDA_VERSION < 13000
    case CUFFT_INCOMPLETE_PARAMETER_LIST:
        return "CUFFT_INCOMPLETE_PARAMETER_LIST";
    case CUFFT_PARSE_ERROR:
        return "CUFFT_PARSE_ERROR";
    case CUFFT_LICENSE_ERROR:
        return "CUFFT_LICENSE_ERROR";
#endif
    default:
        return "<unknown>";
    }
}

#define CHECK_CUDA(func)                                                                                               \
    do                                                                                                                 \
    {                                                                                                                  \
        cudaError_t status = (func);                                                                                   \
        if (status != cudaSuccess)                                                                                     \
        {                                                                                                              \
            fprintf(stderr, "In File %s : CUDA API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,       \
                    cudaGetErrorString(status), status);                                                               \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#define CHECK_CUBLAS(func)                                                                                             \
    do                                                                                                                 \
    {                                                                                                                  \
        cublasStatus_t status = (func);                                                                                \
        if (status != CUBLAS_STATUS_SUCCESS)                                                                           \
        {                                                                                                              \
            fprintf(stderr, "In File %s : CUBLAS API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,     \
                    _cublasGetErrorString(status), status);                                                            \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#define CHECK_CUSOLVER(func)                                                                                           \
    do                                                                                                                 \
    {                                                                                                                  \
        cusolverStatus_t status = (func);                                                                              \
        if (status != CUSOLVER_STATUS_SUCCESS)                                                                         \
        {                                                                                                              \
            fprintf(stderr, "In File %s : CUSOLVER API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,   \
                    _cusolverGetErrorString(status), status);                                                          \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#define CHECK_CUFFT(func)                                                                                              \
    do                                                                                                                 \
    {                                                                                                                  \
        cufftResult_t status = (func);                                                                                 \
        if (status != CUFFT_SUCCESS)                                                                                   \
        {                                                                                                              \
            fprintf(stderr, "In File %s : CUFFT API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,      \
                    _cufftGetErrorString(status), status);                                                             \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#define CHECK_LAST_CUDA_ERROR(msg)                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        cudaError_t status = cudaGetLastError();                                                                       \
        if (status != cudaSuccess)                                                                                     \
        {                                                                                                              \
            fprintf(stderr, "%s(%d) : CUDA error : %s : (%d) %s.\n", __FILE__, __LINE__, msg,                          \
                    static_cast<int>(status), cudaGetErrorString(status));                                             \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#ifdef __DEBUG
#define CHECK_CUDA_SYNC()                                                                                              \
    do                                                                                                                 \
    {                                                                                                                  \
        cudaError_t status = cudaDeviceSynchronize();                                                                  \
        if (status != cudaSuccess)                                                                                     \
        {                                                                                                              \
            fprintf(stderr, "In File %s : CUDA sync failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,      \
                    cudaGetErrorString(status), status);                                                               \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)
#else
#define CHECK_CUDA_SYNC() do {} while (0)
#endif

// cuSOLVER MP support
#ifdef __CUSOLVERMP
#include <cusolverMp.h>

#ifdef __USE_CAL
#include <cal.h>

static const char* _calGetErrorString(calError_t error)
{
    switch (error)
    {
    case CAL_OK:
        return "CAL_OK";
    case CAL_ERROR:
        return "CAL_ERROR";
    case CAL_ERROR_INVALID_PARAMETER:
        return "CAL_ERROR_INVALID_PARAMETER";
    case CAL_ERROR_INTERNAL:
        return "CAL_ERROR_INTERNAL";
    case CAL_ERROR_CUDA:
        return "CAL_ERROR_CUDA";
    case CAL_ERROR_UCC:
        return "CAL_ERROR_UCC";
    case CAL_ERROR_NOT_SUPPORTED:
        return "CAL_ERROR_NOT_SUPPORTED";
    case CAL_ERROR_INPROGRESS:
        return "CAL_ERROR_INPROGRESS";
    default:
        return "<unknown>";
    }
}

#define CHECK_CAL(func)                                                                                                \
    do                                                                                                                 \
    {                                                                                                                  \
        calError_t status = (func);                                                                                    \
        if (status != CAL_OK)                                                                                          \
        {                                                                                                              \
            fprintf(stderr, "In File %s : CAL API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,        \
                    _calGetErrorString(status), status);                                                               \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)
#else // !__USE_CAL (use NCCL)
#include <nccl.h>

#define CHECK_NCCL(func)                                                                                               \
    do                                                                                                                 \
    {                                                                                                                  \
        ncclResult_t status = (func);                                                                                  \
        if (status != ncclSuccess)                                                                                     \
        {                                                                                                              \
            fprintf(stderr, "In File %s : NCCL API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,       \
                    ncclGetErrorString(status), status);                                                               \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)
#endif // __USE_CAL

#endif // __CUSOLVERMP

#endif // __CUDA

#ifdef __ROCM
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <hipfft/hipfft.h>
#include <hipsolver/hipsolver.h>

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
    default:
        return "<unknown>";
    }
}

static const char* _hipsolverGetErrorString(hipsolverStatus_t error)
{
    switch (error)
    {
    case HIPSOLVER_STATUS_SUCCESS:
        return "HIPSOLVER_STATUS_SUCCESS";
    case HIPSOLVER_STATUS_NOT_INITIALIZED:
        return "HIPSOLVER_STATUS_NOT_INITIALIZED";
    case HIPSOLVER_STATUS_ALLOC_FAILED:
        return "HIPSOLVER_STATUS_ALLOC_FAILED";
    case HIPSOLVER_STATUS_INVALID_VALUE:
        return "HIPSOLVER_STATUS_INVALID_VALUE";
    case HIPSOLVER_STATUS_MAPPING_ERROR:
        return "HIPSOLVER_STATUS_MAPPING_ERROR";
    case HIPSOLVER_STATUS_EXECUTION_FAILED:
        return "HIPSOLVER_STATUS_EXECUTION_FAILED";
    case HIPSOLVER_STATUS_INTERNAL_ERROR:
        return "HIPSOLVER_STATUS_INTERNAL_ERROR";
    case HIPSOLVER_STATUS_NOT_SUPPORTED:
        return "HIPSOLVER_STATUS_NOT_SUPPORTED";
    case HIPSOLVER_STATUS_ARCH_MISMATCH:
        return "HIPSOLVER_STATUS_ARCH_MISMATCH";
    case HIPSOLVER_STATUS_HANDLE_IS_NULLPTR:
        return "HIPSOLVER_STATUS_HANDLE_IS_NULLPTR";
    case HIPSOLVER_STATUS_INVALID_ENUM:
        return "HIPSOLVER_STATUS_INVALID_ENUM";
    case HIPSOLVER_STATUS_UNKNOWN:
        return "HIPSOLVER_STATUS_UNKNOWN";
    default:
        return "<unknown>";
    }
}

#define CHECK_CUDA(func)                                                                                               \
    do                                                                                                                 \
    {                                                                                                                  \
        hipError_t status = (func);                                                                                    \
        if (status != hipSuccess)                                                                                      \
        {                                                                                                              \
            fprintf(stderr, "In File %s : HIP API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,        \
                    hipGetErrorString(status), status);                                                                \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#define CHECK_CUBLAS(func)                                                                                             \
    do                                                                                                                 \
    {                                                                                                                  \
        hipblasStatus_t status = (func);                                                                               \
        if (status != HIPBLAS_STATUS_SUCCESS)                                                                          \
        {                                                                                                              \
            fprintf(stderr, "In File %s : HIPBLAS API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,    \
                    _hipblasGetErrorString(status), status);                                                           \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#define CHECK_CUSOLVER(func)                                                                                           \
    do                                                                                                                 \
    {                                                                                                                  \
        hipsolverStatus_t status = (func);                                                                             \
        if (status != HIPSOLVER_STATUS_SUCCESS)                                                                        \
        {                                                                                                              \
            fprintf(stderr, "In File %s : HIPSOLVER API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,  \
                    _hipsolverGetErrorString(status), status);                                                         \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#define CHECK_CUFFT(func)                                                                                              \
    do                                                                                                                 \
    {                                                                                                                  \
        hipfftResult_t status = (func);                                                                                \
        if (status != HIPFFT_SUCCESS)                                                                                  \
        {                                                                                                              \
            fprintf(stderr, "In File %s : HIPFFT API failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,     \
                    _hipfftGetErrorString(status), status);                                                            \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#define CHECK_LAST_CUDA_ERROR(msg)                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        hipError_t status = hipGetLastError();                                                                         \
        if (status != hipSuccess)                                                                                      \
        {                                                                                                              \
            fprintf(stderr, "%s(%d) : HIP error : %s : (%d) %s.\n", __FILE__, __LINE__, msg,                           \
                    static_cast<int>(status), hipGetErrorString(status));                                              \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)

#ifdef __DEBUG
#define CHECK_CUDA_SYNC()                                                                                              \
    do                                                                                                                 \
    {                                                                                                                  \
        hipError_t status = hipDeviceSynchronize();                                                                    \
        if (status != hipSuccess)                                                                                      \
        {                                                                                                              \
            fprintf(stderr, "In File %s : HIP sync failed at line %d with error: %s (%d)\n", __FILE__, __LINE__,       \
                    hipGetErrorString(status), status);                                                                \
            exit(EXIT_FAILURE);                                                                                        \
        }                                                                                                              \
    } while (0)
#else
#define CHECK_CUDA_SYNC() do {} while (0)
#endif

#endif // __ROCM

#endif // DEVICE_CHECK_H
