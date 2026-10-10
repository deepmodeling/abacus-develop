#ifndef ABACUS_LR_STATE_GROUP_H
#define ABACUS_LR_STATE_GROUP_H
#include <complex>
#include <vector>
namespace LR
{
struct GroupMatch
{
    int group = 0;
    double score = 0.0;
    double subspace_score = 0.0;
    std::vector<double> reference_weights;
    double old_coverage = 0.0;
    double new_coverage = 0.0;
    double runner_up = 0.0;
    std::vector<double> singular_values;
};
// Row-major old-group by current-root overlap. JT ranks by the actual reference
// projection weight; otherwise total subspace overlap ranks groups. Ties use the other score.
// No dimension cap: losing old directions is diagnosed, rather than rejected.
GroupMatch match_state_group(const std::vector<std::complex<double>>& overlap,
    int old_dimension, int nstates, const std::vector<std::vector<int>>& groups,
    const std::vector<double>& reference_overlap, bool reference_first);
}
#endif
