// by zhanghaochong 20240529
// The reason that I create a new helper file is the header file of helper_cuda.h is tooooo long

#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HSOLVER_KERNELS_CUDA_HELPER_CUSOLVER_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HSOLVER_KERNELS_CUDA_HELPER_CUSOLVER_H
#ifdef __CUSOLVERMP
#include <cusolverMp.h>

const char* calGetErrorString(calError_t status)
{
    switch (status)
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
        return "CAL UNKNOWN ERROR";
    }
}

#define CAL_CHECK(cmd)                                                                                                 \
    do                                                                                                                 \
    {                                                                                                                  \
        calError_t status = cmd;                                                                                       \
        if (status != CAL_OK)                                                                                          \
        {                                                                                                              \
            fprintf(stderr, "ERROR: %s %s %d\n", calGetErrorString(status), __FILE__, __LINE__);                       \
            abort();                                                                                                   \
        }                                                                                                              \
    } while (0)
#endif

const char* cusolverGetErrorString(cusolverStatus_t status)
{
    switch (status)
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
        return "CUSOLVER UNKNOWN ERROR";
    }
}

#define CUSOLVER_CHECK(cmd)                                                                                            \
    do                                                                                                                 \
    {                                                                                                                  \
        cusolverStatus_t status = cmd;                                                                                 \
        if (status != CUSOLVER_STATUS_SUCCESS)                                                                         \
        {                                                                                                              \
            fprintf(stderr, "ERROR: %s %s %d\n", cusolverGetErrorString(status), __FILE__, __LINE__);                  \
            abort();                                                                                                   \
        }                                                                                                              \
    } while (0)

#endif // W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_HSOLVER_KERNELS_CUDA_HELPER_CUSOLVER_H