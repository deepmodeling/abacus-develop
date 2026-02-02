/**
 * @file cuda_compat.h
 * @brief Compatibility layer for CUDA and NVTX headers across different CUDA Toolkit versions.
 *
 * This header abstracts the differences in NVTX (NVIDIA Tools Extension) header locations
 * between CUDA Toolkit versions.
 *
 * @note Depends on the CUDA_VERSION macro defined in <cuda.h>.
 *
 */

#ifndef CUDA_COMPAT_H_
#define CUDA_COMPAT_H_

#include <iostream> // For std::ostream
#include <stdexcept> // For std::invalid_argument
#include <cuda.h> // defines CUDA_VERSION
#include <cuda_runtime.h>
#include <cufft.h>


// NVTX header for CUDA versions prior to 12.9 vs. 12.9+
// This block ensures the correct NVTX header path is used based on CUDA_VERSION.
// - For CUDA Toolkit < 12.9, the legacy header "nvToolsExt.h" is included.
// - For CUDA Toolkit >= 12.9, the modern header "nvtx3/nvToolsExt.h" is included,
// and NVTX v2 is removed from 12.9.
// This allows NVTX profiling APIs (e.g. nvtxRangePush) to be used consistently
// across different CUDA versions.
// See:
// https://docs.nvidia.com/cuda/archive/12.9.0/cuda-toolkit-release-notes/index.html#id4
#if defined(__CUDA) && defined(__USE_NVTX)
#if CUDA_VERSION < 12090
    #include "nvToolsExt.h"
#else
    #include "nvtx3/nvToolsExt.h"
#endif
#endif

//-------------------------------------------------------------------------------------------------
// Compatibility Layer Declarations
//-------------------------------------------------------------------------------------------------
namespace ModuleBase {
namespace cuda_compat {

/**
 * @brief Prints device information that was deprecated or removed in CUDA 13.0.
 *
 * This function handles properties like clockRate, memoryClockRate, memoryBusWidth,
 * and concurrency flags, which are not available in newer CUDA toolkits.
 *
 * @param os The output stream (e.g., std::cout, std::ofstream).
 * @param prop The cudaDeviceProp structure containing device properties.
 */
void printDeprecatedDeviceInfo(std::ostream& os, const cudaDeviceProp& prop);

/**
 * @brief Prints the device's compute mode using a legacy string mapping.
 *
 * The compute mode display logic is encapsulated here as it relies on aspects
 * of the driver model that have changed.
 *
 * @param os The output stream (e.g., std::cout, std::ofstream).
 * @param prop The cudaDeviceProp structure containing device properties.
 */
void printComputeModeInfo(std::ostream& os, const cudaDeviceProp& prop);

} // namespace cuda_compat
} // namespace ModuleBase

#endif // CUDA_COMPAT_H_
