#include "../nacv_factorization.h"
#include "source_base/parallel_reduce.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include <algorithm>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

TEST(NacvFactorization, ReusesPivotedLuForDifferentRhsCounts)
{
    const std::vector<double> matrix{0.0,
        1.0,
        2.0,
        3.0};
    LR::NacvFactorization factors(matrix, 2);
    std::vector<double> first{4.0, 7.0};
    factors.solve(first, 1);
    EXPECT_NEAR(first[0], 1.0, 1e-14);
    EXPECT_NEAR(first[1], 2.0, 1e-14);
    std::vector<double> second{1.0,
        0.5,
        4.0,
        7.0};
    factors.solve(second, 2);
    EXPECT_NEAR(second[0], -1.0, 1e-14);
    EXPECT_NEAR(second[1], 0.5, 1e-14);
    EXPECT_NEAR(second[2], 1.0, 1e-14);
    EXPECT_NEAR(second[3], 2.0, 1e-14);
}

TEST(NacvFactorization, DifferentGeometryOwnsIndependentFactors)
{
    const std::vector<double> first_matrix{2.0};
    const std::vector<double> second_matrix{4.0};
    LR::NacvFactorization first(first_matrix, 1);
    LR::NacvFactorization second(second_matrix, 1);
    std::vector<double> rhs{1.0};
    second.solve(rhs, 1);
    EXPECT_DOUBLE_EQ(rhs[0], 0.25);
    rhs[0] = 1.0;
    first.solve(rhs, 1);
    EXPECT_DOUBLE_EQ(rhs[0], 0.5);
}

TEST(NacvFactorization, RejectsSingularLu)
{
    const std::vector<double> singular{1.0,
        2.0,
        2.0,
        4.0};
    EXPECT_THROW(LR::NacvFactorization factors(singular, 2), std::runtime_error);
}

#ifdef __MPI
namespace
{
std::vector<double> distribute(const std::vector<double> &full, const Parallel_Blkcyclic &layout)
{
    const std::size_t storage = std::max<std::size_t>(1, layout.size());
    std::vector<double> local(storage, 0.0);
    for (int col = 0; col < layout.ncloc(); ++col)
    {
        for (int row = 0; row < layout.nrloc(); ++row)
        {
            local[col * layout.nrloc() + row] = full[layout.icglb(col) * layout.nrglb() + layout.irglb(row)];
        }
    }
    return local;
}

void check_distributed_reuse(bool cholesky, int dimension)
{
    std::vector<double> matrix(dimension * dimension, 0.0);
    std::vector<double> symmetric(matrix.size(), 0.0);
    for (int col = 0; col < dimension; ++col)
    {
        for (int row = 0; row < dimension; ++row)
        {
            const int index = col * dimension + row;
            symmetric[index] = row == col ? dimension + 2.0 : 0.125;
            matrix[index] = symmetric[index];
            if (cholesky)
            {
                matrix[index] += 0.03125 * (row - col);
            }
        }
    }
    Parallel_Blkcyclic layout;
    LR_Util::setup_2d_division(
        layout,
        1,
        dimension,
        dimension);
    auto local = distribute(matrix, layout);
    LR::NacvFactorization factors(local, layout, cholesky);
    for (const int count : {1, 3, 1})
    {
        Parallel_Blkcyclic rhs_layout;
        LR_Util::setup_2d_division(
            rhs_layout,
            1,
            dimension,
            count,
            layout.ictxt());
        std::vector<double> expected(dimension * count);
        std::vector<double> rhs(expected.size(), 0.0);
        for (int col = 0; col < count; ++col)
        {
            for (int row = 0; row < dimension; ++row)
            {
                expected[col * dimension + row] = 0.25 * (col + 1) - row;
            }
            for (int row = 0; row < dimension; ++row)
            {
                for (int inner = 0; inner < dimension; ++inner)
                {
                    rhs[col * dimension + row] +=
                        symmetric[inner * dimension + row] * expected[col * dimension + inner];
                }
            }
        }
        auto solution = distribute(rhs, rhs_layout);
        factors.solve(solution, rhs_layout);
        std::vector<double> full(expected.size(), 0.0);
        LR_Util::gather_2d_to_full(
            rhs_layout,
            solution.data(),
            full.data(),
            false,
            dimension,
            count);
        for (std::size_t i = 0; i < full.size(); ++i)
        {
            EXPECT_NEAR(full[i], expected[i], 1e-12);
        }
    }
}
} // namespace
TEST(NacvFactorization, DistributedLuReusesFactorsAndPivotLayout) { check_distributed_reuse(false, 5); }
TEST(NacvFactorization, DistributedCholeskyPreservesSymmetrization) { check_distributed_reuse(true, 5); }
TEST(NacvFactorization, DistributedZeroStorageRanksStillParticipate)
{
    check_distributed_reuse(false, 1);
    check_distributed_reuse(true, 1);
}
TEST(NacvFactorization, DistributedFailuresAreReportedOnAllRanks)
{
    Parallel_Blkcyclic layout;
    LR_Util::setup_2d_division(
        layout,
        1,
        2,
        2);
    std::vector<double> matrix{1.0,
        0.0,
        0.0,
        0.0};
    auto local = distribute(matrix, layout);
    EXPECT_THROW(LR::NacvFactorization factors(local, layout, false), std::runtime_error);
    matrix[3] = -1.0;
    local = distribute(matrix, layout);
    EXPECT_THROW(LR::NacvFactorization factors(local, layout, true), std::runtime_error);
}
#endif
