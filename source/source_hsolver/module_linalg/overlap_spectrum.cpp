#include "overlap_spectrum.h"

#include <algorithm>
#include <cmath>

namespace hsolver
{
bool overlap_inverse_sqrt(const std::vector<double>& values,
                          const double relative_cutoff,
                          std::vector<double>& weights)
{
    weights.clear();
    const auto is_finite = [](double value) { return std::isfinite(value); };
    if (values.empty() || !(relative_cutoff >= 0.0 && relative_cutoff < 1.0)
        || !std::all_of(values.begin(), values.end(), is_finite)
        || !std::is_sorted(values.begin(), values.end()) || values.back() <= 0.0
        || (relative_cutoff == 0.0 && values.front() <= 0.0))
    {
        return false;
    }

    weights.resize(values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (values[i] > 0.0 && (relative_cutoff == 0.0 || values[i] / values.back() > relative_cutoff))
        {
            weights[i] = 1.0 / std::sqrt(values[i]);
        }
    }
    return true;
}
} // namespace hsolver
