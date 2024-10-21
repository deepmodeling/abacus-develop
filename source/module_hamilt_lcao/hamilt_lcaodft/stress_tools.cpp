#include "stress_tools.h"
namespace StressTools
{
void stress_fill(const double& lat0_, const double& omega_, ModuleBase::matrix& stress_matrix)
{
    assert(omega_ > 0.0);
    double weight = lat0_ / omega_;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (j > i)
            {
                stress_matrix(j, i) = stress_matrix(i, j);
            }
            stress_matrix(i, j) *= weight;
        }
    }
}
} // namespace StressTools