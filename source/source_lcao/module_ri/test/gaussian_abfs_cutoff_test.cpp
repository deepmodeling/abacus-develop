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
        const double l0_cutoff = Gaussian_Abfs::radial_cutoff(lambda, 0);
        const double l8_cutoff = Gaussian_Abfs::radial_cutoff(lambda, 8);
        require(l0_cutoff > std::sqrt(35.0 / lambda),
                "Gaussian radial cutoff must control the complete l=0 multipole tail");
        require(l8_cutoff > l0_cutoff, "high angular momentum must enlarge the Gaussian radial cutoff");
    }

    bool rejected = false;
    try
    {
        (void)Gaussian_Abfs::radial_cutoff(0.0, 8);
    }
    catch (const std::invalid_argument&)
    {
        rejected = true;
    }
    require(rejected, "non-positive lambda must be rejected");
    return 0;
}
