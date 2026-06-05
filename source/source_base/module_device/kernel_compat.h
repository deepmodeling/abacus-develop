#ifndef KERNEL_COMPAT_H_
#define KERNEL_COMPAT_H_

/**
 * @file kernel_compat.h
 * @brief Device-side kernel compatibility polyfills for older GPU architectures.
 *
 * This is a lightweight header (no heavy includes) for GPU kernel device-side
 * compatibility code. Include this header in .cu/.hip.cu files that need
 * legacy GPU support.
 *
 * Note: The existing cuda_compat.h is for host-side CUDA compatibility
 * (NVTX, deprecated APIs, cuFFT) and includes heavy headers, so we keep
 * this separate.
 */

// atomicAdd for double precision - required for CUDA architectures < 600 (pre-Pascal)
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600 && !defined(__CUDA_ON_DCU)
static __inline__ __device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do
    {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif

#endif // KERNEL_COMPAT_H_
