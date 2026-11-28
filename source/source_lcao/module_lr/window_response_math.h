#ifndef ABACUS_LR_WINDOW_RESPONSE_MATH_H
#define ABACUS_LR_WINDOW_RESPONSE_MATH_H
#include <vector>
namespace LR
{
// Column-major matrices; X and AX span all available virtuals. Only the
// retained columns of the virtual-virtual multiplier can be nonzero.
std::vector<double> window_multipliers(int nocc, int nvirt, int nactive,
    const double* virtual_energies, const double* X, const double* AX);
}
#endif
