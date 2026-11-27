#ifndef HSOLVER_OVERLAP_SPECTRUM_H
#define HSOLVER_OVERLAP_SPECTRUM_H

#include <vector>

namespace hsolver
{
/** @brief Compute inverse square root weights for an ascending overlap spectrum.
 * At zero cutoff every value must be positive. At positive cutoff, weights for
 * values at or below relative_cutoff * max(values) are zero.
 */
bool overlap_inverse_sqrt(const std::vector<double>& values,
                          double relative_cutoff,
                          std::vector<double>& weights);
}

#endif
