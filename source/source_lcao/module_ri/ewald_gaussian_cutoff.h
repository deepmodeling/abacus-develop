#ifndef EWALD_GAUSSIAN_CUTOFF_H
#define EWALD_GAUSSIAN_CUTOFF_H

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace EwaldGaussianCutoff
{
constexpr double historical_real_space_exponent = 17.5;
constexpr double real_space_tail_tolerance = 1.0e-14;
constexpr double historical_reciprocal_exponent = 35.0;
constexpr double reciprocal_tail_tolerance = 1.0e-14;

inline double real_space_multipole_tail(const double exponent, const int angular_momentum)
{
    if (exponent < 0.0)
    {
        throw std::invalid_argument("Ewald Gaussian decay exponent must be non-negative");
    }
    if (angular_momentum < 0)
    {
        throw std::invalid_argument("Ewald Gaussian angular momentum must be non-negative");
    }

    double shape = 0.5;
    double tail = std::erfc(std::sqrt(exponent));
    for (int step = 0; step <= angular_momentum; ++step)
    {
        tail += std::pow(exponent, shape) * std::exp(-exponent) / std::tgamma(shape + 1.0);
        shape += 1.0;
    }
    return tail;
}

inline double real_space_decay_exponent(const int angular_momentum)
{
    if (angular_momentum < 0)
    {
        throw std::invalid_argument("Ewald Gaussian angular momentum must be non-negative");
    }

    double lower = historical_real_space_exponent;
    double upper = lower;
    while (real_space_multipole_tail(upper, angular_momentum) > real_space_tail_tolerance)
    {
        upper += 1.0;
    }
    for (int iteration = 0; iteration != 80; ++iteration)
    {
        const double midpoint = 0.5 * (lower + upper);
        if (real_space_multipole_tail(midpoint, angular_momentum) > real_space_tail_tolerance)
        {
            lower = midpoint;
        }
        else
        {
            upper = midpoint;
        }
    }
    return upper;
}

inline double real_space_radius(const double lambda, const int angular_momentum)
{
    if (lambda <= 0.0)
    {
        throw std::invalid_argument("Ewald Gaussian lambda must be positive");
    }
    return std::sqrt(2.0 * real_space_decay_exponent(angular_momentum) / lambda);
}

inline double reciprocal_tail_scale(const double exponent, const int lmax)
{
    if (exponent < 0.0)
    {
        throw std::invalid_argument("Ewald Gaussian decay exponent must be non-negative");
    }
    if (lmax < 0)
    {
        throw std::invalid_argument("Ewald Gaussian angular momentum must be non-negative");
    }
    const int polynomial_power = std::max(0, lmax - 2);
    return std::pow(exponent, 0.5 * polynomial_power) * std::exp(-exponent);
}

inline double reciprocal_decay_exponent(const int lmax)
{
    if (lmax < 0)
    {
        throw std::invalid_argument("Ewald Gaussian angular momentum must be non-negative");
    }

    double lower = historical_reciprocal_exponent;
    double upper = lower;
    while (reciprocal_tail_scale(upper, lmax) > reciprocal_tail_tolerance)
    {
        upper += 1.0;
    }
    for (int iteration = 0; iteration != 80; ++iteration)
    {
        const double midpoint = 0.5 * (lower + upper);
        if (reciprocal_tail_scale(midpoint, lmax) > reciprocal_tail_tolerance)
        {
            lower = midpoint;
        }
        else
        {
            upper = midpoint;
        }
    }
    return upper;
}

inline double reciprocal_radius(const double lambda, const int lmax)
{
    if (lambda <= 0.0)
    {
        throw std::invalid_argument("Ewald Gaussian lambda must be positive");
    }
    return std::sqrt(reciprocal_decay_exponent(lmax) * lambda);
}
} // namespace EwaldGaussianCutoff

#endif
