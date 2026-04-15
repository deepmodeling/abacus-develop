#ifndef MODULE_BASE_TRUNCATED_FUNC_H
#define MODULE_BASE_TRUNCATED_FUNC_H

#include "source_base/libm/libm.h"
#include <cstdint>
#include <cstring>
#include <complex>

namespace ModuleBase
{

/**
 * @brief Truncated exponential function to avoid underflow.
 *
 * This function returns 0 if the real part of the input is less than -230.0,
 * otherwise it calls ModuleBase::libm::exp(x).
 *
 * @tparam FPTYPE The floating point type (float, double, or complex).
 * @param x The input value.
 * @return FPTYPE The result of the exponential function.
 */
template <typename FPTYPE>
inline FPTYPE truncated_exp(FPTYPE x)
{
    if (std::real(x) < -230.0)
    {
        return static_cast<FPTYPE>(0.0);
    }
    return ModuleBase::libm::exp(x);
}

/**
 * @brief Truncated complementary error function to avoid underflow for large arguments.
 *
 * This function returns 0 if the real part of the input is greater than 20.0,
 * otherwise it calls std::erfc(x).
 *
 * @tparam FPTYPE The floating point type (float, double, or complex).
 * @param x The input value.
 * @return FPTYPE The result of the erfc function.
 */
template <typename FPTYPE>
inline FPTYPE truncated_erfc(FPTYPE x)
{
    if (std::real(x) > 20.0)
    {
        return static_cast<FPTYPE>(0.0);
    }
    return std::erfc(x);
}

/**
 * @brief Truncated value function to avoid underflow.
 *
 * This function returns 0 if the absolute value of the input is less than 1.0e-30,
 * otherwise it returns the input x.
 *
 * @tparam FPTYPE The floating point type (float, double, or complex).
 * @param x The input value.
 * @return FPTYPE The result of the truncation.
 */
/**
 * @brief Truncated value function to avoid underflow.
 *
 * This function modifies the input to 0 if its absolute value is less than 1.0e-30.
 *
 * @tparam FPTYPE The floating point type (float, double, or complex).
 * @param x The input value to be checked and possibly truncated.
 */
template <typename FPTYPE>
inline void truncated_underflow(FPTYPE& x)
{
    if (std::abs(x) < 1.0e-30)
    {
        x = static_cast<FPTYPE>(0.0);
    }
}

template <>
inline void truncated_underflow(double& x)
{
    const uint64_t u = *reinterpret_cast<const uint64_t*>(&x);
    // The exponent bits are 52-62 (11 bits). The bias is 1023.
    // 1e-30 corresponds to -100 in base-2 exponent roughly.
    // 923 = 1023 - 100.
    if (((u >> 52) & 0x7FF) <= 923)
    {
        x = 0.0;
    }
}

template <>
inline void truncated_underflow(float& x)
{
    const uint32_t u = *reinterpret_cast<const uint32_t*>(&x);
    // The exponent bits are 23-30 (8 bits). The bias is 127.
    // 1e-30 corresponds to -100 in base-2 exponent roughly.
    // 27 = 127 - 100.
    if (((u >> 23) & 0xFF) <= 27)
    {
        x = 0.0f;
    }
}

template <typename T>
inline void truncated_underflow(std::complex<T>& x)
{
    T* ptr = reinterpret_cast<T*>(&x);
    truncated_underflow(ptr[0]);
    truncated_underflow(ptr[1]);
}


} // namespace ModuleBase

#endif // MODULE_BASE_TRUNCATED_FUNC_H