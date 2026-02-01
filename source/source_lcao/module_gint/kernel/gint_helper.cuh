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

