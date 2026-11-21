#include "window_response_math.h"
#include <cmath>
#include <stdexcept>
namespace LR
{
std::vector<double> window_multipliers(const int nocc, const int nvirt, const int nactive,
    const double* virtual_energies, const double* X, const double* AX)
{
    if (nocc <= 0 || nactive <= 0 || nactive > nvirt)
    {
        throw std::invalid_argument("Invalid orbital dimensions for window response");
    }
    const int matrix_size = nvirt * nvirt;
    std::vector<double> W(matrix_size, 0.0);
    for (int a = 0; a < nactive; ++a)
    {
        for (int r = nactive; r < nvirt; ++r)
        {
            const double gap = virtual_energies[a] - virtual_energies[r];
            // A spectral projector is not differentiable across a degenerate
            // boundary. Include the entire multiplet instead of clipping gaps.
            if (std::abs(gap) < 1e-8)
            {
                throw std::runtime_error("LR window cuts a degenerate virtual multiplet; increase nvirt");
            }
            double coupling = 0.0;
            for (int i = 0; i < nocc; ++i)
            {
                coupling += X[i * nvirt + a] * AX[i * nvirt + r];
            }
            W[a * nvirt + r] = 2.0 * coupling / gap;
        }
    }
    return W;
}
}
