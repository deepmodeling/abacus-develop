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
        const double expected = std::sqrt(35.0 / lambda);
        require(std::abs(Gaussian_Abfs::radial_cutoff(lambda) - expected) < 1e-14,
                "Gaussian radial cutoff must match sqrt(35/lambda)");
    }

    bool rejected = false;
    try
    {
        (void)Gaussian_Abfs::radial_cutoff(0.0);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "non-positive lambda must be rejected");
    return 0;
}
