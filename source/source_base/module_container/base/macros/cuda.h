#ifndef BASE_MACROS_CUDA_H_
#define BASE_MACROS_CUDA_H_

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <thrust/complex.h>

#include "source_base/module_device/device_check.h"

#define THREADS_PER_BLOCK 256

template <typename T>
struct GetTypeThrust
{
    using type = T;
};

template <>
struct GetTypeThrust<std::complex<float>>
{
    using type = thrust::complex<float>; /**< The return type specialization for std::complex<float>. */
};

template <>
struct GetTypeThrust<std::complex<double>>
{
    using type = thrust::complex<double>; /**< The return type specialization for std::complex<double>. */
};

static inline cublasOperation_t GetCublasOperation(const char& trans)
{
    cublasOperation_t cutrans = {};
    if (trans == 'N')
    {
        cutrans = CUBLAS_OP_N;
    }
    else if (trans == 'T')
    {
        cutrans = CUBLAS_OP_T;
    }
    else if (trans == 'C')
    {
        cutrans = CUBLAS_OP_C;
    }
    return cutrans;
}

template <typename T>
struct GetTypeCuda
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_R_32F;
};
// Specializations of DataTypeToEnum for supported types.
template <>
struct GetTypeCuda<int>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_R_32I;
};
template <>
struct GetTypeCuda<float>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_R_32F;
};
template <>
struct GetTypeCuda<double>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_R_64F;
};
template <>
struct GetTypeCuda<int64_t>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_R_64I;
};
template <>
struct GetTypeCuda<std::complex<float>>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_C_32F;
};
template <>
struct GetTypeCuda<std::complex<double>>
{
    static constexpr cudaDataType cuda_data_type = cudaDataType::CUDA_C_64F;
};

static inline cublasFillMode_t cublas_fill_mode(const char& uplo)
{
    if (uplo == 'U' || uplo == 'u')
        return CUBLAS_FILL_MODE_UPPER;
    else if (uplo == 'L' || uplo == 'l')
        return CUBLAS_FILL_MODE_LOWER;
    else
        throw std::runtime_error("cublas_fill_mode: unknown uplo");
}

static inline cublasDiagType_t cublas_diag_type(const char& diag)
{
    if (diag == 'U' || diag == 'u')
        return CUBLAS_DIAG_UNIT;
    else if (diag == 'N' || diag == 'n')
        return CUBLAS_DIAG_NON_UNIT;
    else
        throw std::runtime_error("cublas_diag_type: unknown diag");
}

static inline cusolverEigMode_t cublas_eig_mode(const char& jobz)
{
    if (jobz == 'N' || jobz == 'n')
        return CUSOLVER_EIG_MODE_NOVECTOR;
    else if (jobz == 'V' || jobz == 'v')
        return CUSOLVER_EIG_MODE_VECTOR;
    else
        throw std::runtime_error("cublas_eig_mode: unknown diag");
}

static inline cusolverEigType_t cublas_eig_type(const int& itype)
{
    if (itype == 1)
        return CUSOLVER_EIG_TYPE_1;
    else if (itype == 2)
        return CUSOLVER_EIG_TYPE_2;
    else
        throw std::runtime_error("cublas_eig_mode: unknown diag");
}

/**
 * @brief Converts a character specifying eigenvalue range to cuSOLVER enum.
 *
 *        'A' or 'a' -> CUSOLVER_EIG_RANGE_ALL: all eigenvalues
 *        'V' or 'v' -> CUSOLVER_EIG_RANGE_V:  values in [vl, vu]
 *        'I' or 'i' -> CUSOLVER_EIG_RANGE_I:  indices in [il, iu]
 *
 * @param range Character indicating selection mode ('A', 'V', 'I')
 * @return Corresponding cusolverEigRange_t enum value
 * @throws std::runtime_error if character is invalid
 */
static inline cusolverEigRange_t cublas_eig_range(const char& range)
{
    if (range == 'A' || range == 'a')
        return CUSOLVER_EIG_RANGE_ALL;
    else if (range == 'V' || range == 'v')
        return CUSOLVER_EIG_RANGE_V;
    else if (range == 'I' || range == 'i')
        return CUSOLVER_EIG_RANGE_I;
    else
        throw std::runtime_error("cublas_eig_range: unknown range '" + std::string(1, range) + "'");
}

#endif // BASE_MACROS_CUDA_H_
