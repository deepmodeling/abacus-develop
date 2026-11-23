#include "../nacv_solver.h"
#include "source_base/opt_cg.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_reduce.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include <algorithm>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
bool owns_rhs() { return ModuleBase::world_comm_domain().rank() == 0; }

void diagonal_operator(const double *input, double *output)
{
    output[0] = 2.0 * input[0];
    output[1] = 5.0 * input[1];
}
} // namespace

TEST(ZeqCgPrecision, TightToleranceMatchesIndependentDiagonalInverse)
{
    double rhs[] = {0.0, 0.0};
    if (owns_rhs())
    {
        rhs[0] = 1.0;
        rhs[1] = -3.0;
    }
    double solution[] = {9.0, 7.0};
    LR::solve_nac_cg(
        solution,
        rhs,
        2,
        1,
        diagonal_operator,
        false,
        1e-12,
        10);
    EXPECT_NEAR(solution[0], rhs[0] / 2.0, 1e-12);
    EXPECT_NEAR(solution[1], rhs[1] / 5.0, 1e-12);
}

TEST(ZeqCgPrecision, RequestedThresholdChangesSmallResponseResolution)
{
    double rhs[] = {0.0, 0.0};
    if (owns_rhs())
    {
        rhs[0] = 1e-7;
        rhs[1] = -2e-7;
    }
    double loose[] = {7.0, 7.0};
    double tight[] = {7.0, 7.0};
    LR::solve_nac_cg(
        loose,
        rhs,
        2,
        1,
        diagonal_operator,
        false,
        1e-6,
        10);
    LR::solve_nac_cg(
        tight,
        rhs,
        2,
        1,
        diagonal_operator,
        false,
        1e-12,
        10);
    EXPECT_DOUBLE_EQ(loose[0], 0.0);
    EXPECT_DOUBLE_EQ(loose[1], 0.0);
    EXPECT_NEAR(tight[0], rhs[0] / 2.0, 1e-14);
    EXPECT_NEAR(tight[1], rhs[1] / 5.0, 1e-14);
}

TEST(ZeqCgPrecision, ConditionedSystemMatchesIndependentInverseAndTrueResidual)
{
    const std::vector<double> eigenvalues = {1.0,
        10.0,
        100.0,
        1000.0,
        10000.0,
        100000.0,
        1000000.0,
        10000000.0};
    const int size = static_cast<int>(eigenvalues.size());
    std::vector<double> rhs(size, 0.0);
    if (owns_rhs())
    {
        for (int index = 0; index < size; ++index)
        {
            rhs[index] = index + 1.0;
        }
    }
    std::vector<double> solution(size, 0.0);
    const std::vector<double> original_rhs = rhs;
    const std::function<void(const double *, double *)> diagonal = [&eigenvalues, size](
                                                                       const double *input, double *output)
    {
        for (int index = 0; index < size; ++index)
        {
            output[index] = eigenvalues[index] * input[index];
        }
    };
    LR::solve_nac_cg(
        solution.data(),
        rhs.data(),
        size,
        1,
        diagonal,
        false,
        1e-12,
        400);
    double residual_squared = 0.0;
    for (int index = 0; index < size; ++index)
    {
        const double difference = eigenvalues[index] * solution[index] - rhs[index];
        residual_squared += difference * difference;
        EXPECT_NEAR(solution[index], rhs[index] / eigenvalues[index], 1e-12);
    }
    Parallel_Reduce::reduce_all(residual_squared);
    EXPECT_LT(std::sqrt(residual_squared), 1e-12);
    EXPECT_EQ(rhs, original_rhs);
}

/** @brief Require bounded update progress when valid norm reductions round differently. */
TEST(ZeqCgPrecision, ResidualReplacementMakesProgressAtTheThreshold)
{
    const int size = 32;
    for (std::uint32_t seed = 1; seed <= 64; ++seed)
    {
        std::vector<double> rhs(size, 0.0);
        std::uint32_t state = seed;
        if (owns_rhs())
        {
            for (int index = 0; index < size; ++index)
            {
                state = 1664525u * state + 1013904223u;
                const double mantissa = 1.0 + static_cast<double>(state & 65535u) / 65536.0;
                const int exponent = static_cast<int>((state >> 16) % 24) - 12;
                rhs[index] = std::ldexp(mantissa, exponent);
            }
        }
        double scalar_squared = 0.0;
        for (int index = 0; index < size; ++index)
        {
            scalar_squared += rhs[index] * rhs[index];
        }
        Parallel_Reduce::reduce_all(scalar_squared);
        double *const rhs_data = rhs.data();
        double blas_squared = BlasConnector::dot(
            size,
            rhs_data,
            1,
            rhs_data,
            1);
        Parallel_Reduce::reduce_all(blas_squared);
        const double scalar_norm = std::sqrt(scalar_squared);
        const double blas_norm = std::sqrt(blas_squared);
        if (!(blas_norm < scalar_norm))
        {
            continue;
        }
        std::vector<double> solution(size, 0.0);
        double *const solution_data = solution.data();
        int operator_calls = 0;
        const std::function<void(const double *, double *)> identity = [&](const double *input, double *output)
        {
            ++operator_calls;
            if (operator_calls > 4)
            {
                // Bound a regression rather than letting the unit test hang.
                throw std::runtime_error("A residual replacement repeated without a CG update");
            }
            for (int index = 0; index < size; ++index)
            {
                output[index] = input[index];
            }
        };
        const bool test_force = false;
        const int maximum_updates = 1;
        const int states = 1;
        EXPECT_NO_THROW(
            LR::solve_nac_cg(
                solution_data,
                rhs_data,
                size,
                states,
                identity,
                test_force,
                scalar_norm,
                maximum_updates));
        EXPECT_LE(operator_calls, 4);
        return;
    }
    GTEST_SKIP() << "This BLAS/compiler did not produce a residual-norm rounding boundary";
}

TEST(ZeqCgPrecision, IterationLimitRejectsAnUncertifiedVector)
{
    double rhs[] = {0.0, 0.0};
    if (owns_rhs())
    {
        rhs[0] = 1.0;
        rhs[1] = -3.0;
    }
    double solution[] = {0.0, 0.0};
    EXPECT_THROW(
        LR::solve_nac_cg(
            solution,
            rhs,
            2,
            1,
            diagonal_operator,
            false,
            1e-12,
            1),
        std::runtime_error);
}

TEST(ZeqCgPrecision, CertifiesTheActualOperatorAtTheReturnedVector)
{
    double rhs[] = {0.0};
    if (owns_rhs())
    {
        rhs[0] = 1.0;
    }
    double solution[] = {0.0};
    const std::function<void(const double *, double *)> nonlinear = [](const double *input, double *output)
    {
        const double scale = std::abs(input[0]) > .75 ? 2.0 : 1.0;
        output[0] = scale * input[0];
    };
    EXPECT_THROW(
        LR::solve_nac_cg(
            solution,
            rhs,
            1,
            1,
            nonlinear,
            false,
            1e-10,
            10),
        std::runtime_error);
}

TEST(ZeqCgPrecision, EmptyLocalRanksParticipateInTheGlobalResidual)
{
    const int local_size = owns_rhs() ? 2 : 0;
    std::vector<double> rhs(local_size, 1.0);
    std::vector<double> solution(local_size, 7.0);
    const std::function<void(const double *, double *)> diagonal = [local_size](const double *input, double *output)
    {
        for (int index = 0; index < local_size; ++index)
        {
            output[index] = (index + 2.0) * input[index];
        }
    };
    LR::solve_nac_cg(
        solution.data(),
        rhs.data(),
        local_size,
        1,
        diagonal,
        false,
        1e-12,
        10);
    for (int index = 0; index < local_size; ++index)
    {
        EXPECT_NEAR(solution[index], 1.0 / (index + 2.0), 1e-12);
    }
}

TEST(ZeqCgPrecision, ReusesTheSuccessfulTrueResidualCertification)
{
    double rhs[] = {0.0};
    if (owns_rhs())
    {
        rhs[0] = 1.0;
    }
    double solution[] = {9.0};
    int applications = 0;
    const auto diagonal = [&applications](const double *input, double *output)
    {
        ++applications;
        output[0] = 2.0 * input[0];
    };
    LR::solve_nac_cg(
        solution,
        rhs,
        1,
        1,
        diagonal,
        false,
        1e-12,
        4);
    EXPECT_EQ(applications, 2);
    EXPECT_DOUBLE_EQ(solution[0], rhs[0] * 0.5);
}

TEST(ZeqCgPrecision, ZeroRhsIsCertifiedExactlyOnce)
{
    double rhs[] = {0.0};
    double solution[] = {9.0};
    int applications = 0;
    const auto diagonal = [&applications](const double *input, double *output)
    {
        ++applications;
        output[0] = 2.0 * input[0];
    };
    LR::solve_nac_cg(
        solution,
        rhs,
        1,
        1,
        diagonal,
        false,
        1e-12,
        4);
    EXPECT_EQ(applications, 1);
    EXPECT_DOUBLE_EQ(solution[0], 0.0);
}

TEST(ZeqCgPrecision, CoupledMultipleRightHandSidesMatchKnownSolutions)
{
    const int local_size = owns_rhs() ? 3 : 0;
    const int states = 2;
    const int local_rhs_size = local_size * states;
    std::vector<double> rhs(local_rhs_size, 0.0);
    std::vector<double> expected;
    if (owns_rhs())
    {
        // The same positive-definite Hessian acts on both state-major RHS rows.
        // L = {{4,1,0},{1,3,1},{0,1,2}}.
        rhs = {2.0,
            -2.0,
            4.0,
            -3.5,
            2.5,
            4.5};
        expected = {1.0,
            -2.0,
            3.0,
            -1.0,
            0.5,
            2.0};
    }
    const std::vector<double> original_rhs = rhs;
    const std::size_t rhs_size = rhs.size();
    std::vector<double> solution(rhs_size, 7.0);
    const auto hessian = [local_size, states](const double *input, double *output)
    {
        if (local_size == 0)
        {
            return;
        }
        for (int state = 0; state < states; ++state)
        {
            const int offset = state * local_size;
            output[offset] = 4.0 * input[offset] + input[offset + 1];
            output[offset + 1] = input[offset] + 3.0 * input[offset + 1] + input[offset + 2];
            output[offset + 2] = input[offset + 1] + 2.0 * input[offset + 2];
        }
    };
    double *const solution_data = solution.data();
    double *const rhs_data = rhs.data();
    LR::solve_nac_cg(
        solution_data,
        rhs_data,
        local_size,
        states,
        hessian,
        false,
        1e-12,
        10);
    EXPECT_EQ(rhs, original_rhs);
    for (std::size_t index = 0; index < expected.size(); ++index)
    {
        EXPECT_NEAR(solution[index], expected[index], 1e-12);
    }
    std::vector<double> applied(rhs_size, 0.0);
    double *const applied_data = applied.data();
    hessian(solution_data, applied_data);
    double residual_squared = 0.0;
    for (std::size_t index = 0; index < rhs.size(); ++index)
    {
        const double residual = applied[index] - rhs[index];
        residual_squared += residual * residual;
    }
    Parallel_Reduce::reduce_all(residual_squared);
    EXPECT_LT(std::sqrt(residual_squared), 1e-12);
}

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
    const std::size_t local_size = layout.size();
    const std::size_t storage = std::max<std::size_t>(1, local_size);
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
    const int matrix_size = dimension * dimension;
    std::vector<double> matrix(matrix_size, 0.0);
    std::vector<double> symmetric(matrix_size, 0.0);
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
        const int rhs_size = dimension * count;
        std::vector<double> expected(rhs_size);
        std::vector<double> rhs(rhs_size, 0.0);
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
        std::vector<double> full(rhs_size, 0.0);
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
