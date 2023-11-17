#ifndef ATEN_KERNELS_TEST_OP_TEST_UTILS_H_
#define ATEN_KERNELS_TEST_OP_TEST_UTILS_H_
#include <gtest/gtest.h>
#include "module_base/module_container/ATen/kernels/blas.h"
#include "module_base/module_container/ATen/kernels/lapack.h"

namespace container {
namespace test_utils {

# if __CUDA || __ROCM
using ComplexTypes = ::testing::Types<
        std::tuple<std::complex<float>, DEVICE_CPU>, std::tuple<std::complex<float>, DEVICE_GPU>,
        std::tuple<std::complex<double>, DEVICE_CPU>,  std::tuple<std::complex<double>, DEVICE_GPU>>;
using Types = ::testing::Types<
        std::tuple<float, DEVICE_CPU>, std::tuple<float, DEVICE_GPU>,
        std::tuple<double, DEVICE_CPU>, std::tuple<double, DEVICE_GPU>,
        std::tuple<std::complex<float>, DEVICE_CPU>, std::tuple<std::complex<float>, DEVICE_GPU>,
        std::tuple<std::complex<double>, DEVICE_CPU>,  std::tuple<std::complex<double>, DEVICE_GPU>>;
#else 
using ComplexTypes = ::testing::Types<
        std::tuple<std::complex<float>, DEVICE_CPU>,
        std::tuple<std::complex<double>, DEVICE_CPU>>;
using Types = ::testing::Types<
        std::tuple<float, DEVICE_CPU>, 
        std::tuple<double, DEVICE_CPU>,
        std::tuple<std::complex<float>, DEVICE_CPU>,
        std::tuple<std::complex<double>, DEVICE_CPU>>;
#endif 

static inline void init_blas_handle() {
    #if __CUDA || __ROCM
        kernels::createGpuBlasHandle();
    #endif
}

static inline void delete_blas_handle() {
    #if __CUDA || __ROCM
        kernels::destroyGpuBlasHandle();
    #endif
}

static inline void init_cusolver_handle() {
    #if __CUDA || __ROCM
        kernels::createGpuSolverHandle();
    #endif
}

static inline void delete_cusolver_handle() {
    #if __CUDA || __ROCM
        kernels::destroyGpuSolverHandle();
    #endif
}

} // namespace test_utils
} // namespace container

#endif // ATEN_KERNELS_TEST_OP_TEST_UTILS_H_