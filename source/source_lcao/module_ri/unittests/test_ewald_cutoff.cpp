#include "../ewald_gaussian_cutoff.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition)
    {
        std::cerr << "FAILED: " << message << '\n';
        std::exit(1);
    }
}
} // namespace

int main()
{
    const double historical_exponent = 17.5;
    const double l8_tail = EwaldGaussianCutoff::real_space_multipole_tail(historical_exponent, 8);
    require(std::abs(l8_tail - 1.396683e-2) < 1.0e-8,
            "the historical real-space cutoff must expose the l=8 multipole loss");

    const double exponent = EwaldGaussianCutoff::real_space_decay_exponent(8);
    require(exponent > historical_exponent, "high angular momentum must enlarge the real-space cutoff");
    require(EwaldGaussianCutoff::real_space_multipole_tail(exponent, 8)
                <= EwaldGaussianCutoff::real_space_tail_tolerance,
            "the l=8 real-space multipole tail must satisfy the configured tolerance");

    const double lambda = 1.7;
    const double rcut = EwaldGaussianCutoff::real_space_radius(lambda, 8);
    require(std::abs(0.5 * lambda * rcut * rcut - exponent) < 1.0e-12,
            "real-space radius must use exp(-lambda*r^2/2)");

    bool rejected = false;
    try
    {
        (void)EwaldGaussianCutoff::real_space_radius(0.0, 8);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "non-positive lambda must be rejected");

    const int lmax = 16;
    const double reciprocal_exponent = EwaldGaussianCutoff::reciprocal_decay_exponent(lmax);
    require(reciprocal_exponent > 35.0, "L=16 must enlarge the historical reciprocal cutoff");
    require(EwaldGaussianCutoff::reciprocal_tail_scale(reciprocal_exponent, lmax)
                <= EwaldGaussianCutoff::reciprocal_tail_tolerance,
            "the L=16 reciprocal Gaussian tail must satisfy the configured tolerance");

    const double gmax = EwaldGaussianCutoff::reciprocal_radius(lambda, lmax);
    require(std::abs(gmax * gmax / lambda - reciprocal_exponent) < 1.0e-12,
            "reciprocal radius must use exp(-G^2/lambda)");

    return 0;
}
