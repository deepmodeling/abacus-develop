#include "../gaussian_abfs.h"

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
}

int main()
{
    for (const double lambda: {0.5, 1.0, 2.0})
    {
        // The strict-2D kernel controls the complete multipole tail at 1e-14.
        const double l0_cutoff = Gaussian_Abfs::radial_cutoff(lambda, 0, 2);
        const double l8_cutoff = Gaussian_Abfs::radial_cutoff(lambda, 8, 2);
        require(l0_cutoff > std::sqrt(35.0 / lambda),
                "strict-2D Gaussian radial cutoff must control the complete l=0 multipole tail");
        require(l8_cutoff > l0_cutoff,
                "strict-2D high angular momentum must enlarge the Gaussian radial cutoff");

        // 3D keeps the historical 17.5 real-space and 35.0 reciprocal
        // exponents at every angular momentum.
        const double historical_cutoff = std::sqrt(35.0 / lambda);
        for (const int angular_momentum: {0, 4, 8})
        {
            require(std::abs(Gaussian_Abfs::radial_cutoff(lambda, angular_momentum, 3) - historical_cutoff)
                        < 1.0e-14 * historical_cutoff,
                    "3D Gaussian radial cutoff must keep the historical exponent");
        }
        const double historical_gmax = std::sqrt(35.0 * lambda);
        require(std::abs(EwaldGaussianCutoff::reciprocal_radius(lambda, 0, 3) - historical_gmax)
                    < 1.0e-14 * historical_gmax,
                "3D reciprocal cutoff must keep the historical exponent");
        require(std::abs(EwaldGaussianCutoff::reciprocal_radius(lambda, 8, 3) - historical_gmax)
                    < 1.0e-14 * historical_gmax,
                "3D reciprocal cutoff must not depend on the angular momentum");
        require(EwaldGaussianCutoff::reciprocal_radius(lambda, 8, 2) > historical_gmax,
                "strict-2D reciprocal cutoff must account for the higher-l tail");
    }

    bool rejected = false;
    try
    {
        (void)Gaussian_Abfs::radial_cutoff(0.0, 8, 2);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "non-positive lambda must be rejected");
    return 0;
}
