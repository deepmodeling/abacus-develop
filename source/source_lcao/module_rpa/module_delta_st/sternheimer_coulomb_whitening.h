#ifndef STERNHEIMER_COULOMB_WHITENING_H
#define STERNHEIMER_COULOMB_WHITENING_H

#include <vector>

namespace ModuleRI
{

struct SternheimerCoulombWhitening
{
    int raw_dimension = 0;
    int retained_rank = 0;
    int discarded_rank = 0;
    double relative_threshold = 0.0;
    double largest_eigenvalue = 0.0;
    double smallest_retained_eigenvalue = 0.0;
    double max_orthonormality_error = 0.0;
    std::vector<double> eigenvalues;
    std::vector<double> transform; ///< row-major raw_dimension x retained_rank
};

SternheimerCoulombWhitening make_sternheimer_coulomb_whitening(const std::vector<double>& metric,
                                                               int dimension,
                                                               double relative_threshold);

std::vector<double> apply_sternheimer_channel_transform(const std::vector<double>& rows_by_raw_channel,
                                                        int row_count,
                                                        int raw_dimension,
                                                        const SternheimerCoulombWhitening& whitening);

} // namespace ModuleRI

#endif
