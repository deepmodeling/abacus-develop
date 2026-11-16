#include "source_base/opt_cg.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_reduce.h"
#include "source_lcao/module_lr/nacv_cg.h"

#include <cmath>
#include <cstdint>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <utility>
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
