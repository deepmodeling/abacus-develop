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

#include <cuda.h> // defines CUDA_VERSION

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

#endif // CUDA_COMPAT_H_
