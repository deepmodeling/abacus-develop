#ifndef GAUSSIAN_ABFS_2D_INTEGRALS_H
#define GAUSSIAN_ABFS_2D_INTEGRALS_H

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

// Closed-form in-plane momentum integrals of the 2D Ewald Gaussian
// compensation for auxiliary functions.
//
// In the mixed representation (in-plane momentum q, real z) every needed
// integral reduces to z moments of a shifted Gaussian, which are analytic
// (physicists' Hermite polynomials), times in-plane factors that are
// polynomials in q arranged by real solid harmonics. The base Coulomb
// integrals i2/di2_dtau carry the erfc form of the z integration. Only even
// in-plane kernel powers q^p with p >= -2 occur; p = -2 (the Coulomb
// singularity) is handled by synthetic division onto the base integrals.
namespace GaussianAbfs2D
{
constexpr double pi = 3.141592653589793238462643383279502884;

// Scaled complementary error function erfcx(x) = e^{x^2} erfc(x) for
// non-negative arguments. The direct product overflows once
// x^2 > log(DBL_MAX) while erfc simultaneously underflows, although the
// scaled value stays finite; large arguments use the asymptotic expansion.
inline double erfcx_nonneg(const double x)
{
    if (x * x < 709.0)
    {
        return std::exp(x * x) * std::erfc(x);
    }
    const double x2 = x * x;
    return (1.0 / (std::sqrt(pi) * x))
           * (1.0 - 0.5 / x2 + 0.75 / (x2 * x2) - 1.875 / (x2 * x2 * x2)
              + 6.5625 / (x2 * x2 * x2 * x2));
}

// The z-shifted erfc products of the base integral, each premultiplied by
// e^{beta q^2 +- q a} so the cancellation between that exponential and the
// erfc argument is carried analytically:
//     plus  = e^{-a^2/(4 beta)} erfcx(sqrt(beta) q + a/(2 sqrt(beta))),
//     minus = e^{-a^2/(4 beta)} erfcx(sqrt(beta) q - a/(2 sqrt(beta))).
// The minus argument can be negative; erfcx(-v) = 2 e^{v^2} - erfcx(v) then
// keeps the a^2/(4 beta) cancellation inside the exponent as well.
inline void scaled_erfc_pair(const double q,
                             const double a,
                             const double beta,
                             double& plus,
                             double& minus)
{
    const double sqrt_beta = std::sqrt(beta);
    const double x = sqrt_beta * q;
    const double half_shift = a / (2.0 * sqrt_beta);
    const double decay = std::exp(-half_shift * half_shift);
    plus = decay * erfcx_nonneg(x + half_shift);
    if (x >= half_shift)
    {
        minus = decay * erfcx_nonneg(x - half_shift);
    }
    else
    {
        const double v = half_shift - x;
        minus = 2.0 * std::exp(beta * q * q - q * a) - decay * erfcx_nonneg(v);
    }
}

// Base integral: the z integration of the Gaussian-compensated Coulomb at
// in-plane momentum q and center offset tau_z,
//     i2 = (pi/2q) e^{beta q^2} [ e^{q|tau_z|} erfc(sqrt(beta) q + |tau_z|/2 sqrt(beta))
//                                + e^{-q|tau_z|} erfc(sqrt(beta) q - |tau_z|/2 sqrt(beta)) ],
// which reduces to (pi/q) e^{beta q^2} erfc(sqrt(beta) q) at tau_z = 0.
inline double i2(const double q, const double tau_z, const double beta)
{
    if (q < 1.0e-12)
    {
        // Diverges like pi/q as q->0; the Gamma sample of this channel is
        // dropped by the lattice sum and restored analytically by the
        // callers together with the chi correction.
        return 0.0;
    }

    const double sqrt_beta = std::sqrt(beta);
    const double x = sqrt_beta * q;
    if (std::abs(tau_z) < 1.0e-14)
    {
        return pi / q * erfcx_nonneg(x);
    }

    const double a = std::abs(tau_z);
    double plus = 0.0;
    double minus = 0.0;
    scaled_erfc_pair(q, a, beta, plus, minus);
    const double result = pi / (2.0 * q) * (plus + minus);
    if (!std::isfinite(result))
    {
        throw std::runtime_error("GaussianAbfs2D::i2 produced a nonfinite base integral");
    }
    return result;
}

inline double di2_dtau(const double q, const double tau_z, const double beta)
{
    if (q < 1.0e-12)
    {
        // Unique q->0 limit, independent of the in-plane direction: the
        // normal-dipole (s-pz) channel of the Gamma term. At tau_z = 0 this
        // reduces to zero, consistent with the odd symmetry below.
        return -pi * std::erf(tau_z / (2.0 * std::sqrt(beta)));
    }

    if (std::abs(tau_z) < 1.0e-14)
    {
        return 0.0;
    }

    const double a = std::abs(tau_z);
    double plus = 0.0;
    double minus = 0.0;
    scaled_erfc_pair(q, a, beta, plus, minus);
    const double derivative_a = 0.5 * pi * (plus - minus);
    const double result = tau_z < 0.0 ? -derivative_a : derivative_a;
    if (!std::isfinite(result))
    {
        throw std::runtime_error("GaussianAbfs2D::di2_dtau produced a nonfinite derivative");
    }
    return result;
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

namespace detail
{
inline double factorial(const int n)
{
    if (n < 0) throw std::invalid_argument("factorial argument must be non-negative");
    return std::tgamma(static_cast<double>(n + 1));
}

inline double physicists_hermite(const int degree, const double x)
{
    if (degree < 0) throw std::invalid_argument("Hermite degree must be non-negative");
    if (degree == 0) return 1.0;
    if (degree == 1) return 2.0 * x;
    double hm2 = 1.0;
    double hm1 = 2.0 * x;
    for (int n = 2; n <= degree; ++n)
    {
        const double h = 2.0 * x * hm1 - 2.0 * static_cast<double>(n - 1) * hm2;
        hm2 = hm1;
        hm1 = h;
    }
    return hm1;
}

inline std::complex<double> i_power(const int degree)
{
    switch ((degree % 4 + 4) % 4)
    {
        case 0: return {1.0, 0.0};
        case 1: return {0.0, 1.0};
        case 2: return {-1.0, 0.0};
        default: return {0.0, -1.0};
    }
}

// z moment of the shifted Gaussian, <(z)^degree exp(-beta (z-tau_z)^2)>,
// expressed through physicists' Hermite polynomials; i^degree carries the
// Fourier convention of the mixed representation.
inline std::complex<double> gaussian_fourier_moment(const int degree,
                                                    const double tau_z,
                                                    const double beta)
{
    if (beta <= 0.0) throw std::invalid_argument("Gaussian exponent must be positive");
    const double scaled_tau = tau_z / std::sqrt(beta);
    const double radial_prefactor = std::pow(beta, -0.5 * static_cast<double>(degree + 1));
    const double fourier_prefactor = std::sqrt(pi) / std::pow(2.0, degree);
    return radial_prefactor * fourier_prefactor
           * physicists_hermite(degree, 0.5 * scaled_tau)
           * std::exp(-0.25 * scaled_tau * scaled_tau) * i_power(degree);
}

inline std::vector<std::complex<double>> complex_solid_harmonic_coefficients(const int L,
                                                                             const int m,
                                                                             const double qx,
                                                                             const double qy)
{
    if (L < 0 || m < 0 || m > L) throw std::invalid_argument("invalid solid-harmonic indices");
    std::vector<std::complex<double>> coefficients(static_cast<std::size_t>(L + 1),
                                                   std::complex<double>{});
    const double normalization
        = std::sqrt(static_cast<double>(2 * L + 1) / (4.0 * pi)
                    * factorial(L - m) / factorial(L + m));
    const double condon_shortley = (m % 2 == 0) ? 1.0 : -1.0;
    // std::pow of the zero in-plane vector with exponent zero is not
    // portable (exp(0*log(0)) is NaN on some standard libraries but exactly
    // one on others); spell out the m=0 factor. The q=0 evaluation point of
    // the Gamma term hits this on every call.
    const std::complex<double> in_plane_power
        = m == 0 ? std::complex<double>(1.0, 0.0)
                 : std::pow(std::complex<double>(qx, qy), static_cast<double>(m));
    const std::complex<double> angular_prefactor = normalization * condon_shortley * in_plane_power;
    const double q_sq = qx * qx + qy * qy;
    for (int k = 0; k <= L / 2; ++k)
    {
        const int z_degree = L - 2 * k - m;
        if (z_degree < 0) continue;
        const double legendre_derivative
            = ((k % 2 == 0) ? 1.0 : -1.0) * factorial(2 * L - 2 * k)
              / (std::pow(2.0, L) * factorial(k) * factorial(L - k) * factorial(z_degree));
        for (int j = 0; j <= k; ++j)
        {
            const int degree = z_degree + 2 * j;
            const double binomial = factorial(k) / (factorial(j) * factorial(k - j));
            coefficients[static_cast<std::size_t>(degree)]
                += angular_prefactor * legendre_derivative * binomial * std::pow(q_sq, k - j);
        }
    }
    return coefficients;
}

inline std::vector<double> real_solid_harmonic_coefficients(const int L,
                                                            const int M,
                                                            const double qx,
                                                            const double qy)
{
    if (M < 0 || M > 2 * L) throw std::invalid_argument("invalid real solid-harmonic index");
    const int m = M == 0 ? 0 : (M + 1) / 2;
    const auto complex_coefficients = complex_solid_harmonic_coefficients(L, m, qx, qy);
    std::vector<double> coefficients(complex_coefficients.size(), 0.0);
    if (M == 0)
    {
        for (std::size_t i = 0; i < coefficients.size(); ++i)
            coefficients[i] = complex_coefficients[i].real();
    }
    else if (M % 2 == 1)
    {
        for (std::size_t i = 0; i < coefficients.size(); ++i)
            coefficients[i] = std::sqrt(2.0) * complex_coefficients[i].real();
    }
    else
    {
        for (std::size_t i = 0; i < coefficients.size(); ++i)
            coefficients[i] = std::sqrt(2.0) * complex_coefficients[i].imag();
    }
    return coefficients;
}
} // namespace detail

// Master integral: in-plane kernel q^power times the (L,M) real solid
// harmonic of the auxiliary function, with the z dependence of the Gaussian
// integrated analytically. The solid harmonic is first expanded into
// coefficients indexed by z degree; then either
//   power = -2: synthetic division of the polynomial by q^2 maps everything
//               onto the Coulomb base integrals i2 / di2_dtau, or
//   power >= 0 (even): binomial expansion of (q_parallel^2 + qz^2)^{p/2}
//               maps onto the Gaussian z moments.
inline std::complex<double> integral_lm(const int L,
                                        const int M,
                                        const double power,
                                        const double qx,
                                        const double qy,
                                        const double tau_z,
                                        const double beta)
{
    if (beta <= 0.0) throw std::invalid_argument("Gaussian exponent must be positive");
    const int integer_power = static_cast<int>(std::lround(power));
    if (std::abs(power - static_cast<double>(integer_power)) > 1.0e-12
        || integer_power < -2 || integer_power % 2 != 0)
    {
        throw std::invalid_argument("2D Gaussian kernel power must be an even integer >= -2");
    }

    const double q_sq = qx * qx + qy * qy;
    const double q = std::sqrt(q_sq);
    const auto solid = detail::real_solid_harmonic_coefficients(L, M, qx, qy);
    std::vector<std::complex<double>> polynomial(solid.begin(), solid.end());

    if (integer_power == -2)
    {
        std::vector<std::complex<double>> quotient(polynomial.size(), std::complex<double>{});
        for (int degree = static_cast<int>(polynomial.size()) - 1; degree >= 2; --degree)
        {
            const auto leading = polynomial[static_cast<std::size_t>(degree)];
            quotient[static_cast<std::size_t>(degree - 2)] += leading;
            polynomial[static_cast<std::size_t>(degree)] = 0.0;
            polynomial[static_cast<std::size_t>(degree - 2)] -= q_sq * leading;
        }
        std::complex<double> result{};
        for (std::size_t degree = 0; degree < quotient.size(); ++degree)
            result += quotient[degree]
                      * detail::gaussian_fourier_moment(static_cast<int>(degree), tau_z, beta);
        result += polynomial[0] * i2(q, tau_z, beta);
        if (polynomial.size() > 1)
            result += polynomial[1] * std::complex<double>(0.0, -di2_dtau(q, tau_z, beta));
        return result;
    }

    const int exponent = integer_power / 2;
    std::vector<std::complex<double>> expanded(
        polynomial.size() + static_cast<std::size_t>(2 * exponent), std::complex<double>{});
    for (std::size_t degree = 0; degree < polynomial.size(); ++degree)
    {
        for (int j = 0; j <= exponent; ++j)
        {
            const double binomial
                = detail::factorial(exponent) / (detail::factorial(j) * detail::factorial(exponent - j));
            expanded[degree + static_cast<std::size_t>(2 * j)]
                += polynomial[degree] * binomial * std::pow(q_sq, exponent - j);
        }
    }
    std::complex<double> result{};
    for (std::size_t degree = 0; degree < expanded.size(); ++degree)
        result += expanded[degree]
                  * detail::gaussian_fourier_moment(static_cast<int>(degree), tau_z, beta);
    return result;
}
} // namespace GaussianAbfs2D

#endif
