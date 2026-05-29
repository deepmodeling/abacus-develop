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
// gemm_vec_traits<T> -- wide-LDS primitive for the V1 K-inner inner loop.
//
// VK  = how many T elements pack into one 16-byte LDS (4 for FP32, 2 for FP64)
// vec_t = the 16-byte CUDA vector type used for the LDS
// PAD = K-stride padding that makes (BLK_K + PAD) * sizeof(T) a multiple of
//       16 *and* keeps the warp's idx-strided shmem access bank-conflict-free
//       (gcd((BLK_K+PAD) % 32, 32) == VK).
//
// The load is issued as one *reinterpret_cast<vec_t*>(&sA(m,k)); the
// component fan-out is done by unpack(). FP64 needs the explicit cast --
// the compiler's auto-vectorizer is reliable for float4 but not for
// double2; per-component .x/.y/.z/.w writes guarantee the LDS.{64,128}
// SASS forms emit.
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

