#pragma once
#include <cstdio>

// if exponent is an integer between 0 and 5 (the most common cases in gint) and
// and exp is a variable that cannot be determined at compile time (which means the compiler cannot optimize the code),
// pow_int is much faster than std::pow
template<typename T>
__forceinline__ __device__ T pow_int(const T base, const int exp)
{
    switch (exp)
    {
    case 0:
        return 1.0;
    case 1:
        return base;
    case 2:
        return base * base;
    case 3:
        return base * base * base;
    case 4:
        return base * base * base * base;
    case 5:
        return base * base * base * base * base;
    default:
        double result = std::pow(base, exp);
        return result;
    }
}

template<typename T>
__forceinline__ __device__ T warpReduceSum(T val)
{
    val += __shfl_xor_sync(0xffffffff, val, 16, 32);
    val += __shfl_xor_sync(0xffffffff, val, 8, 32);
    val += __shfl_xor_sync(0xffffffff, val, 4, 32);
    val += __shfl_xor_sync(0xffffffff, val, 2, 32);
    val += __shfl_xor_sync(0xffffffff, val, 1, 32);
    return val;
}

inline int ceil_div(const int a, const int b)
{
    return a / b + (a % b != 0 && (a ^ b) > 0);
}

// ---------------------------------------------------------------------------
// gemm_vec_traits<T> -- the wide-load primitive used by the GEMM inner loop,
// which reads VK consecutive K elements per shared-memory load instead of one.
//
//   VK    = number of T elements in one 16-byte load (4 for FP32, 2 for FP64)
//   vec_t = the 16-byte vector type used for that load (float4 / double2)
//   PAD   = padding added to the shared-memory K-stride so that (BLK_K + PAD)
//           elements span a whole number of 16-byte words, keeping the
//           vectorized shared-memory loads aligned and spreading the warp's
//           strided reads across banks.
//
// The load is one *reinterpret_cast<vec_t*>(&sA(m, k)); unpack() then fans the
// vector out into the per-thread registers. The explicit per-component copy is
// deliberate: nvcc reliably vectorizes float4 but not double2, and writing
// .x/.y(/.z/.w) by hand guarantees the wide load instruction is emitted.
// ---------------------------------------------------------------------------
template <typename T> struct gemm_vec_traits;

template <> struct gemm_vec_traits<float>
{
    using vec_t = float4;
    static constexpr int VK = 4;
    static constexpr int PAD = 4;
    __forceinline__ __device__
    static void unpack(const vec_t& v, float* d)
    {
        d[0] = v.x; d[1] = v.y; d[2] = v.z; d[3] = v.w;
    }
};

template <> struct gemm_vec_traits<double>
{
    using vec_t = double2;
    static constexpr int VK = 2;
    static constexpr int PAD = 2;
    __forceinline__ __device__
    static void unpack(const vec_t& v, double* d)
    {
        d[0] = v.x; d[1] = v.y;
    }
};

