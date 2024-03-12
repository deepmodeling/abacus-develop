#ifndef BASE_MACROS_ROCM_H_
#define BASE_MACROS_ROCM_H_

#include <complex>
#include <hipblas/hipblas.h>
#include <hip/hip_runtime.h>
#include <hipsolver/hipsolver.h>

#if defined(__HCC__) || defined(__HIP__)
#include <thrust/complex.h>
#endif // defined(__HCC__) || defined(__HIP__)

#define THREADS_PER_BLOCK 256

#if defined(__HCC__) || defined(__HIP__)
template <typename T>
struct GetTypeThrust {
    using type = T;
};

template <>
struct GetTypeThrust<std::complex<float>> {
    using type = thrust::complex<float>; /**< The return type specialization for std::complex<float>. */
};

template <>
struct GetTypeThrust<std::complex<double>> {
    using type = thrust::complex<double>; /**< The return type specialization for std::complex<float>. */
};
#endif // defined(__HCC__) || defined(__HIP__)

static inline
hipblasOperation_t GetHipblasOperation(const char& trans) {
    hipblasOperation_t hip_trans = {};
    if (trans == 'N') {
        hip_trans = HIPBLAS_OP_N;
    } else if (trans == 'T') {
        hip_trans = HIPBLAS_OP_T;
    } else if (trans == 'C') {
        hip_trans = HIPBLAS_OP_C;
    } else {
        // Handle invalid input or provide a default behavior.
        hip_trans = HIPBLAS_OP_N;
    }
    return hip_trans;
}

template <typename T>
struct GetTypeRocm {
    static constexpr hipDataType hip_data_type = HIP_R_32F;
};

// Specializations of GetTypeRocm for supported types.
template <>
struct GetTypeRocm<int> {
    static constexpr hipDataType hip_data_type = HIP_R_32F;
};
template <>
struct GetTypeRocm<float> {
    static constexpr hipDataType hip_data_type = HIP_R_32F;
};
template <>
struct GetTypeRocm<double> {
    static constexpr hipDataType hip_data_type = HIP_R_64F;
};
template <>
struct GetTypeRocm<int64_t> {
    static constexpr hipDataType hip_data_type = HIP_R_64F;
};
template <>
struct GetTypeRocm<std::complex<float>> {
    static constexpr hipDataType hip_data_type = HIP_C_32F;
};
template <>
struct GetTypeRocm<std::complex<double>> {
    static constexpr hipDataType hip_data_type = HIP_C_64F;
};

static inline
hipblasFillMode_t hipblas_fill_mode(const char& uplo) {
    if (uplo == 'U' || uplo == 'u')
        return HIPBLAS_FILL_MODE_UPPER;
    else if (uplo == 'L' || uplo == 'l')
        return HIPBLAS_FILL_MODE_LOWER;
    else
        throw std::runtime_error("hipblas_fill_mode: unknown uplo");
}

static inline
hipblasDiagType_t hipblas_diag_type(const char& diag) {
    if (diag == 'U' || diag == 'u')
        return HIPBLAS_DIAG_UNIT;
    else if (diag == 'N' || diag == 'n')
        return HIPBLAS_DIAG_NON_UNIT;
    else
        throw std::runtime_error("hipblas_diag_type: unknown diag");
}

static inline
hipsolverEigMode_t hipblas_eig_mode(const char& jobz) {
    if (jobz == 'N' || jobz == 'n')
        return HIPSOLVER_EIG_MODE_NOVECTOR;
    else if (jobz == 'V' || jobz == 'v')
        return HIPSOLVER_EIG_MODE_VECTOR;
    else
        throw std::runtime_error("hipblas_eig_mode: unknown diag");
}

static inline
hipsolverEigType_t hipblas_eig_type(const int& itype) {
    if (itype == 1)
        return HIPSOLVER_EIG_TYPE_1;
    else if (itype == 2)
        return HIPSOLVER_EIG_TYPE_2;
    else
        throw std::runtime_error("hipblas_eig_mode: unknown diag");
}

static inline
hipsolverFillMode_t hipsolver_fill_mode(const char& uplo) {
    if (uplo == 'U' || uplo == 'u')
        return HIPSOLVER_FILL_MODE_UPPER;
    else if (uplo == 'L' || uplo == 'l')
        return HIPSOLVER_FILL_MODE_LOWER;
    else
        throw std::runtime_error("hipsolver_fill_mode: unknown uplo");
}

// hipSOLVER API errors
static const char* hipsolverGetErrorEnum(hipsolverStatus_t error) {
    switch (error) {
        case HIPSOLVER_STATUS_SUCCESS:
            return "HIPSOLVER_STATUS_SUCCESS";
        case HIPSOLVER_STATUS_NOT_INITIALIZED:
            return "HIPSOLVER_STATUS_NOT_INITIALIZED";
        case HIPSOLVER_STATUS_ALLOC_FAILED:
            return "HIPSOLVER_STATUS_ALLOC_FAILED";
        case HIPSOLVER_STATUS_INVALID_VALUE:
            return "HIPSOLVER_STATUS_INVALID_VALUE";
        case HIPSOLVER_STATUS_ARCH_MISMATCH:
            return "HIPSOLVER_STATUS_ARCH_MISMATCH";
        case HIPSOLVER_STATUS_MAPPING_ERROR:
            return "HIPSOLVER_STATUS_MAPPING_ERROR";
        case HIPSOLVER_STATUS_EXECUTION_FAILED:
            return "HIPSOLVER_STATUS_EXECUTION_FAILED";
        case HIPSOLVER_STATUS_INTERNAL_ERROR:
            return "HIPSOLVER_STATUS_INTERNAL_ERROR";
        case HIPSOLVER_STATUS_NOT_SUPPORTED:
            return "HIPSOLVER_STATUS_NOT_SUPPORTED ";
        case HIPSOLVER_STATUS_INVALID_ENUM:
            return "HIPSOLVER_STATUS_INVALID_ENUM";
        default:
            return "Unknown hipsolverStatus_t message";
    }
}

inline void hipsolverAssert(hipsolverStatus_t code, const char* file, int line, bool abort = true)
{
    if (code != HIPSOLVER_STATUS_SUCCESS)
    {
        fprintf(stderr, "hipSOLVER Assert: %s %s %d\n", hipsolverGetErrorEnum(code), file, line);
        if (abort)
            exit(code);
    }
}

// hipSOLVER API errors
static const char * hipblasGetErrorEnum(hipblasStatus_t error) {
    switch (error) {
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
        default:
            return "Unknown";
    }
}

inline void hipblasAssert(hipblasStatus_t code, const char *file, int line, bool abort=true) {
    if (code != HIPBLAS_STATUS_SUCCESS) {
        fprintf(stderr, "Unexpected hipBLAS Error: %s %s %d\n", hipblasGetErrorEnum(code), file, line);
        if (abort) exit(code);
    }
}

#define hipsolverErrcheck(res) { hipsolverAssert((res), __FILE__, __LINE__); }

#define hipblasErrcheck(res) { hipblasAssert((res), __FILE__, __LINE__); }

// ROCM API errors
#define hipErrcheck(res) {                                              \
    if (res != hipSuccess) {                                            \
        fprintf(stderr, " Unexpected Device Error %s:%d: %s, %s\n", __FILE__, __LINE__, hipGetErrorName(res), hipGetErrorString(res)); \
        exit(res);                                                      \
    }                                                                   \
}

#endif // BASE_MACROS_ROCM_H_