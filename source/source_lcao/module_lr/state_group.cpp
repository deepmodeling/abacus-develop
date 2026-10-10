#include "state_group.h"
#include "source_base/module_external/lapack_connector.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace LR
{
GroupMatch match_state_group(const std::vector<std::complex<double>>& overlap,
    const int old_dimension, const int nstates, const std::vector<std::vector<int>>& groups,
    const std::vector<double>& reference_overlap, const bool reference_first)
{
    if (old_dimension < 1 || groups.empty()
        || overlap.size() != static_cast<size_t>(old_dimension) * nstates
        || reference_overlap.size() != static_cast<size_t>(nstates))
    {
        throw std::invalid_argument("LR group overlap: invalid dimensions");
    }
    GroupMatch result;
    double best = -1.0;
    double best_secondary = -1.0;
    result.reference_weights.resize(groups.size());
    for (size_t ig = 0; ig < groups.size(); ++ig)
    {
        const auto& group = groups[ig];
        const int dimension = group.size();
        if (dimension == 0) { throw std::invalid_argument("LR group overlap: empty group"); }
        double weight = 0.0;
        double reference = 0.0;
        for (const int root : group)
        {
            if (root < 0 || root >= nstates) { throw std::invalid_argument("LR group overlap: invalid root"); }
            reference += reference_overlap[root] * reference_overlap[root];
            for (int a = 0; a < old_dimension; ++a) { weight += std::norm(overlap[a * nstates + root]); }
        }
        // Symmetric normalization penalizes unrelated extra roots as well as lost directions.
        const double subspace_score = weight / std::max(old_dimension, dimension);
        result.reference_weights[ig] = reference;
        const double score = reference_first ? reference : subspace_score;
        const double secondary = reference_first ? subspace_score : reference;
        const bool tied = std::abs(score - best) <= 1e-10;
        if (score < best - 1e-10 || (tied && secondary <= best_secondary))
        {
            result.runner_up = std::max(result.runner_up, score);
            continue;
        }
        result.runner_up = std::max(result.runner_up, best);
        best = score;
        best_secondary = secondary;
        result.group = ig;
        result.score = score;
        result.subspace_score = subspace_score;
        result.old_coverage = weight / old_dimension;
        result.new_coverage = weight / dimension;
    }
    const auto& selected = groups[result.group];
    const int rank = std::min(old_dimension, static_cast<int>(selected.size()));
    std::vector<std::complex<double>> gram(rank * rank, 0.0);
    for (int col = 0; col < rank; ++col)
    {
        for (int row = 0; row < rank; ++row)
        {
            if (old_dimension <= static_cast<int>(selected.size()))
            {
                for (const int root : selected)
                {
                    gram[col * rank + row] += overlap[row * nstates + root]
                        * std::conj(overlap[col * nstates + root]);
                }
            }
            else
            {
                for (int a = 0; a < old_dimension; ++a)
                {
                    gram[col * rank + row] += std::conj(overlap[a * nstates + selected[row]])
                        * overlap[a * nstates + selected[col]];
                }
            }
        }
    }
    // Eigenvalues of the smaller Hermitian Gram matrix are squared singular values.
    result.singular_values.resize(rank);
    const int lwork = std::max(1, 2 * rank);
    const int real_work_size = std::max(1, 3 * rank);
    std::vector<std::complex<double>> work(lwork);
    std::vector<double> real_work(real_work_size);
    const char job = 'N';
    const char triangle = 'U';
    int info = 0;
    zheev_(&job, &triangle, &rank, gram.data(), &rank, result.singular_values.data(),
        work.data(), &lwork, real_work.data(), &info);
    if (info != 0) { throw std::runtime_error("LR group overlap: singular spectrum failed"); }
    for (double& value : result.singular_values) { value = std::sqrt(std::max(0.0, value)); }
    std::reverse(result.singular_values.begin(), result.singular_values.end());
    return result;
}
}
