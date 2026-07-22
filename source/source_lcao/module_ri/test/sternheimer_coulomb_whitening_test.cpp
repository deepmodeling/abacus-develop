#include "source_lcao/module_ri/sternheimer_coulomb_whitening.h"

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{

std::vector<double> gram(const std::vector<double>& columns, const int rows, const int count)
{
    std::vector<double> result(static_cast<std::size_t>(count * count), 0.0);
    for (int left = 0; left != count; ++left)
    {
        for (int right = 0; right != count; ++right)
        {
            for (int row = 0; row != rows; ++row)
            {
                result[static_cast<std::size_t>(left * count + right)]
                    += columns[static_cast<std::size_t>(row * count + left)]
                       * columns[static_cast<std::size_t>(row * count + right)];
            }
        }
    }
    return result;
}

std::vector<double> right_multiply(const std::vector<double>& left,
                                   const int rows,
                                   const int inner,
                                   const std::vector<double>& right,
                                   const int columns)
{
    std::vector<double> result(static_cast<std::size_t>(rows * columns), 0.0);
    for (int row = 0; row != rows; ++row)
    {
        for (int column = 0; column != columns; ++column)
        {
            for (int index = 0; index != inner; ++index)
            {
                result[static_cast<std::size_t>(row * columns + column)]
                    += left[static_cast<std::size_t>(row * inner + index)]
                       * right[static_cast<std::size_t>(index * columns + column)];
            }
        }
    }
    return result;
}

std::vector<double> row_projector(const std::vector<double>& columns, const int rows, const int count)
{
    std::vector<double> result(static_cast<std::size_t>(rows * rows), 0.0);
    for (int left = 0; left != rows; ++left)
    {
        for (int right = 0; right != rows; ++right)
        {
            for (int column = 0; column != count; ++column)
            {
                result[static_cast<std::size_t>(left * rows + right)]
                    += columns[static_cast<std::size_t>(left * count + column)]
                       * columns[static_cast<std::size_t>(right * count + column)];
            }
        }
    }
    return result;
}

std::vector<double> left_transpose_multiply(const std::vector<double>& left,
                                            const int rows,
                                            const int inner,
                                            const std::vector<double>& right,
                                            const int columns)
{
    std::vector<double> result(static_cast<std::size_t>(inner * columns), 0.0);
    for (int left_column = 0; left_column != inner; ++left_column)
    {
        for (int right_column = 0; right_column != columns; ++right_column)
        {
            for (int row = 0; row != rows; ++row)
            {
                result[static_cast<std::size_t>(left_column * columns + right_column)]
                    += left[static_cast<std::size_t>(row * inner + left_column)]
                       * right[static_cast<std::size_t>(row * columns + right_column)];
            }
        }
    }
    return result;
}

} // namespace

TEST(SternheimerCoulombWhitening, ProducesIdentityInTheRetainedCoulombMetric)
{
    const std::vector<double> metric = {4.0, 1.0, 1.0, 2.0};
    const ModuleRI::SternheimerCoulombWhitening whitening
        = ModuleRI::make_sternheimer_coulomb_whitening(metric, 2, 1.0e-12);

    EXPECT_EQ(whitening.raw_dimension, 2);
    EXPECT_EQ(whitening.retained_rank, 2);
    EXPECT_LT(whitening.max_orthonormality_error, 1.0e-12);

    const std::vector<double> metric_times_w
        = right_multiply(metric, 2, 2, whitening.transform, whitening.retained_rank);
    const std::vector<double> identity = left_transpose_multiply(whitening.transform,
                                                                 2,
                                                                 whitening.retained_rank,
                                                                 metric_times_w,
                                                                 whitening.retained_rank);
    EXPECT_NEAR(identity[0], 1.0, 1.0e-12);
    EXPECT_NEAR(identity[1], 0.0, 1.0e-12);
    EXPECT_NEAR(identity[2], 0.0, 1.0e-12);
    EXPECT_NEAR(identity[3], 1.0, 1.0e-12);
}

TEST(SternheimerCoulombWhitening, IsInvariantToInvertibleRawAuxiliaryChanges)
{
    const std::vector<double> raw_potentials = {
        1.0, 0.0,
        0.5, 1.0,
        0.0, 2.0,
    };
    const std::vector<double> change = {2.0, 0.5, -0.25, 1.5};
    const std::vector<double> changed_potentials = right_multiply(raw_potentials, 3, 2, change, 2);
    const std::vector<double> raw_metric = gram(raw_potentials, 3, 2);
    const std::vector<double> changed_metric = gram(changed_potentials, 3, 2);

    const auto raw_whitening = ModuleRI::make_sternheimer_coulomb_whitening(raw_metric, 2, 1.0e-12);
    const auto changed_whitening
        = ModuleRI::make_sternheimer_coulomb_whitening(changed_metric, 2, 1.0e-12);
    const auto raw_orthonormal
        = ModuleRI::apply_sternheimer_channel_transform(raw_potentials, 3, 2, raw_whitening);
    const auto changed_orthonormal
        = ModuleRI::apply_sternheimer_channel_transform(changed_potentials, 3, 2, changed_whitening);

    const auto raw_projector = row_projector(raw_orthonormal, 3, raw_whitening.retained_rank);
    const auto changed_projector = row_projector(changed_orthonormal, 3, changed_whitening.retained_rank);
    ASSERT_EQ(raw_projector.size(), changed_projector.size());
    for (std::size_t index = 0; index != raw_projector.size(); ++index)
    {
        EXPECT_NEAR(raw_projector[index], changed_projector[index], 1.0e-12);
    }
}

TEST(SternheimerCoulombWhitening, TruncatesNearNullAndRejectsMateriallyNegativeDirections)
{
    const auto truncated
        = ModuleRI::make_sternheimer_coulomb_whitening({4.0, 0.0, 0.0, 1.0e-14}, 2, 1.0e-12);
    EXPECT_EQ(truncated.retained_rank, 1);
    EXPECT_EQ(truncated.discarded_rank, 1);

    EXPECT_THROW(ModuleRI::make_sternheimer_coulomb_whitening({1.0, 0.0, 0.0, -1.0e-3}, 2, 1.0e-12),
                 std::runtime_error);
}
