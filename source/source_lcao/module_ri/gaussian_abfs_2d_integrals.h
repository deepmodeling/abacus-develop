#ifndef GAUSSIAN_ABFS_2D_INTEGRALS_H
#define GAUSSIAN_ABFS_2D_INTEGRALS_H

#include <cmath>
#include <complex>

namespace GaussianAbfs2D
{
constexpr double pi = 3.141592653589793238462643383279502884;

inline double i2(const double q, const double tau_z, const double beta)
{
    if (q < 1.0e-12)
    {
        return 0.0;
    }

    const double sqrt_beta = std::sqrt(beta);
    const double x = sqrt_beta * q;
    if (std::abs(tau_z) < 1.0e-14)
    {
        return pi / q * std::exp(x * x) * std::erfc(x);
    }

    const double a = std::abs(tau_z);
    const double term_plus = std::exp(q * a) * std::erfc(x + a / (2.0 * sqrt_beta));
    const double term_minus = std::exp(-q * a) * std::erfc(x - a / (2.0 * sqrt_beta));
    return pi / (2.0 * q) * std::exp(x * x) * (term_plus + term_minus);
}

inline double di2_dtau(const double q, const double tau_z, const double beta)
{
    if (q < 1.0e-12 || std::abs(tau_z) < 1.0e-14)
    {
        return 0.0;
    }

    const double a = std::abs(tau_z);
    const double sqrt_beta = std::sqrt(beta);
    const double x = sqrt_beta * q;
    const double term_plus = std::exp(q * a) * std::erfc(x + a / (2.0 * sqrt_beta));
    const double term_minus = std::exp(-q * a) * std::erfc(x - a / (2.0 * sqrt_beta));
    const double derivative_a = 0.5 * pi * std::exp(x * x) * (term_plus - term_minus);
    return tau_z < 0.0 ? -derivative_a : derivative_a;
}

inline std::complex<double> coulomb_l1(const int m,
                                       const double qx,
                                       const double qy,
                                       const double tau_z,
                                       const double beta)
{
    const double q = std::hypot(qx, qy);
    const double coefficient = std::sqrt(3.0 / (4.0 * pi));
    if (m == 0)
    {
        return std::complex<double>(0.0, -coefficient * di2_dtau(q, tau_z, beta));
    }
    if (m == 1)
    {
        return -coefficient * qx * i2(q, tau_z, beta);
    }
    if (m == 2)
    {
        return -coefficient * qy * i2(q, tau_z, beta);
    }
    return 0.0;
}
} // namespace GaussianAbfs2D

#endif
