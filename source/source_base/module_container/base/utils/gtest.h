#ifndef BASE_UTILS_GTEST_H_
#define BASE_UTILS_GTEST_H_
#include <gtest/gtest.h>
#include <ATen/kernels/blas.h>
#include <ATen/kernels/lapack.h>

namespace base {
namespace utils {

#if __CUDA || __ROCM
using ComplexTypes = ::testing::Types<
        std::tuple<std::complex<float>, ct::DEVICE_CPU>, std::tuple<std::complex<float>, ct::DEVICE_GPU>,
        std::tuple<std::complex<double>, ct::DEVICE_CPU>,  std::tuple<std::complex<double>, ct::DEVICE_GPU>>;
using Types = ::testing::Types<
        std::tuple<float, ct::DEVICE_CPU>, std::tuple<float, ct::DEVICE_GPU>,
        std::tuple<double, ct::DEVICE_CPU>, std::tuple<double, ct::DEVICE_GPU>,
        std::tuple<std::complex<float>, ct::DEVICE_CPU>, std::tuple<std::complex<float>, ct::DEVICE_GPU>,
        std::tuple<std::complex<double>, ct::DEVICE_CPU>,  std::tuple<std::complex<double>, ct::DEVICE_GPU>>;
#else 
using ComplexTypes = ::testing::Types<
        std::tuple<std::complex<float>, ct::DEVICE_CPU>,
        std::tuple<std::complex<double>, ct::DEVICE_CPU>>;
using Types = ::testing::Types<
        std::tuple<float, ct::DEVICE_CPU>, 
        std::tuple<double, ct::DEVICE_CPU>,
        std::tuple<std::complex<float>, ct::DEVICE_CPU>,
        std::tuple<std::complex<double>, ct::DEVICE_CPU>>;
#endif 

static inline void init_blas_handle() {
    #if __CUDA || __ROCM
        ct::kernels::createGpuBlasHandle();
    #endif
}

static inline void delete_blas_handle() {
    #if __CUDA || __ROCM
        ct::kernels::destroyGpuBlasHandle();
    #endif
}

static inline void init_cusolver_handle() {
    #if __CUDA || __ROCM
        ct::kernels::createGpuSolverHandle();
    #endif
}

static inline void delete_cusolver_handle() {
    #if __CUDA || __ROCM
        ct::kernels::destroyGpuSolverHandle();
    #endif
}

} // namespace utils
} // namespace base

#endif // BASE_UTILS_GTEST_H_