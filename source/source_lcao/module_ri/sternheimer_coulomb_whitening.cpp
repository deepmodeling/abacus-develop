#include "source_lcao/module_ri/sternheimer_coulomb_whitening.h"

#include "source_base/module_external/lapack_connector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace ModuleRI
{
namespace
{

std::vector<double> diagonalize_symmetric(std::vector<double>& matrix, const int dimension)
{
    const char jobz = 'V';
    const char uplo = 'U';
    const int lda = dimension;
    int info = 0;
    std::vector<double> eigenvalues(static_cast<std::size_t>(dimension), 0.0);
    double work_query = 0.0;
    const int query = -1;

    dsyev_(&jobz,
           &uplo,
           &dimension,
           matrix.data(),
           &lda,
           eigenvalues.data(),
           &work_query,
           &query,
           &info);
    if (info != 0 || !std::isfinite(work_query))
    {
        throw std::runtime_error("Sternheimer Coulomb whitening eigensolver workspace query failed.");
    }
    const int lwork = std::max(1, static_cast<int>(std::ceil(work_query)));
    std::vector<double> work(static_cast<std::size_t>(lwork), 0.0);
    dsyev_(&jobz,
           &uplo,
           &dimension,
           matrix.data(),
           &lda,
           eigenvalues.data(),
           work.data(),
           &lwork,
           &info);
    if (info != 0)
    {
        throw std::runtime_error("Sternheimer Coulomb whitening eigensolver failed.");
    }
    return eigenvalues;
}

void validate_shape(const std::vector<double>& values, const int rows, const int columns, const char* label)
{
    if (rows < 0 || columns < 0
        || values.size() != static_cast<std::size_t>(rows) * static_cast<std::size_t>(columns))
    {
        throw std::invalid_argument(std::string(label) + " has inconsistent dimensions.");
    }
}

} // namespace

SternheimerCoulombWhitening make_sternheimer_coulomb_whitening(const std::vector<double>& metric,
                                                               const int dimension,
                                                               const double relative_threshold)
{
    if (dimension <= 0)
    {
        throw std::invalid_argument("Sternheimer Coulomb whitening dimension must be positive.");
    }
    validate_shape(metric, dimension, dimension, "Sternheimer Coulomb metric");
    if (!std::isfinite(relative_threshold) || relative_threshold <= 0.0 || relative_threshold >= 1.0)
    {
        throw std::invalid_argument(
            "Sternheimer Coulomb whitening relative threshold must be finite and between zero and one.");
    }

    std::vector<double> symmetric(metric.size(), 0.0);
    for (int row = 0; row != dimension; ++row)
    {
        for (int column = 0; column != dimension; ++column)
        {
            const double left = metric[static_cast<std::size_t>(row * dimension + column)];
            const double right = metric[static_cast<std::size_t>(column * dimension + row)];
            if (!std::isfinite(left) || !std::isfinite(right))
            {
                throw std::invalid_argument("Sternheimer Coulomb metric contains a non-finite value.");
            }
            symmetric[static_cast<std::size_t>(row * dimension + column)] = 0.5 * (left + right);
        }
    }
    const std::vector<double> original_symmetric = symmetric;
    std::vector<double> eigenvalues = diagonalize_symmetric(symmetric, dimension);
    const double largest = eigenvalues.back();
    if (!std::isfinite(largest) || largest <= 0.0)
    {
        throw std::runtime_error("Sternheimer Coulomb metric has no positive eigenspace.");
    }
    const double cutoff = relative_threshold * largest;
    const double roundoff_floor
        = 64.0 * std::numeric_limits<double>::epsilon() * static_cast<double>(dimension) * largest;
    const double negative_tolerance = std::max(cutoff, roundoff_floor);
    if (eigenvalues.front() < -negative_tolerance)
    {
        throw std::runtime_error("Sternheimer Coulomb metric has a materially negative eigenvalue.");
    }

    int first_retained = 0;
    while (first_retained != dimension && eigenvalues[static_cast<std::size_t>(first_retained)] <= cutoff)
    {
        ++first_retained;
    }
    if (first_retained == dimension)
    {
        throw std::runtime_error("Sternheimer Coulomb whitening threshold discarded every eigenvector.");
    }

    SternheimerCoulombWhitening result;
    result.raw_dimension = dimension;
    result.retained_rank = dimension - first_retained;
    result.discarded_rank = first_retained;
    result.relative_threshold = relative_threshold;
    result.largest_eigenvalue = largest;
    result.smallest_retained_eigenvalue = eigenvalues[static_cast<std::size_t>(first_retained)];
    result.eigenvalues = eigenvalues;
    result.transform.assign(static_cast<std::size_t>(dimension)
                                * static_cast<std::size_t>(result.retained_rank),
                            0.0);

    // DSYEV returns eigenvectors as columns in column-major storage.
    for (int retained = 0; retained != result.retained_rank; ++retained)
    {
        const int eigenvector = first_retained + retained;
        const double inverse_sqrt
            = 1.0 / std::sqrt(eigenvalues[static_cast<std::size_t>(eigenvector)]);
        for (int raw = 0; raw != dimension; ++raw)
        {
            result.transform[static_cast<std::size_t>(raw * result.retained_rank + retained)]
                = symmetric[static_cast<std::size_t>(eigenvector * dimension + raw)] * inverse_sqrt;
        }
    }

    std::vector<double> metric_times_transform(
        static_cast<std::size_t>(dimension) * static_cast<std::size_t>(result.retained_rank),
        0.0);
    for (int row = 0; row != dimension; ++row)
    {
        for (int retained = 0; retained != result.retained_rank; ++retained)
        {
            for (int column = 0; column != dimension; ++column)
            {
                metric_times_transform[static_cast<std::size_t>(row * result.retained_rank + retained)]
                    += original_symmetric[static_cast<std::size_t>(row * dimension + column)]
                       * result.transform[static_cast<std::size_t>(column * result.retained_rank + retained)];
            }
        }
    }

    double max_error = 0.0;
    for (int left = 0; left != result.retained_rank; ++left)
    {
        for (int right = 0; right != result.retained_rank; ++right)
        {
            double value = 0.0;
            for (int row = 0; row != dimension; ++row)
            {
                value += result.transform[static_cast<std::size_t>(row * result.retained_rank + left)]
                         * metric_times_transform[static_cast<std::size_t>(row * result.retained_rank + right)];
            }
            const double expected = left == right ? 1.0 : 0.0;
            max_error = std::max(max_error, std::abs(value - expected));
        }
    }
    result.max_orthonormality_error = max_error;
    if (!std::isfinite(max_error) || max_error > 1.0e-8)
    {
        throw std::runtime_error("Sternheimer Coulomb whitening failed its W^T V W identity check.");
    }
    return result;
}

std::vector<double> apply_sternheimer_channel_transform(const std::vector<double>& rows_by_raw_channel,
                                                        const int row_count,
                                                        const int raw_dimension,
                                                        const SternheimerCoulombWhitening& whitening)
{
    validate_shape(rows_by_raw_channel,
                   row_count,
                   raw_dimension,
                   "Sternheimer raw-channel matrix");
    if (raw_dimension != whitening.raw_dimension || whitening.retained_rank <= 0)
    {
        throw std::invalid_argument("Sternheimer channel transform dimension does not match its whitening metric.");
    }
    validate_shape(whitening.transform,
                   whitening.raw_dimension,
                   whitening.retained_rank,
                   "Sternheimer Coulomb whitening transform");

    std::vector<double> transformed(
        static_cast<std::size_t>(row_count) * static_cast<std::size_t>(whitening.retained_rank),
        0.0);
    for (int row = 0; row != row_count; ++row)
    {
        for (int retained = 0; retained != whitening.retained_rank; ++retained)
        {
            double value = 0.0;
            for (int raw = 0; raw != raw_dimension; ++raw)
            {
                value += rows_by_raw_channel[static_cast<std::size_t>(row * raw_dimension + raw)]
                         * whitening.transform[static_cast<std::size_t>(raw * whitening.retained_rank + retained)];
            }
            transformed[static_cast<std::size_t>(row * whitening.retained_rank + retained)] = value;
        }
    }
    return transformed;
}

} // namespace ModuleRI
