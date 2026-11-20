/** @file test_nacv_tools.cpp
 * @brief Exercise NAC solvers, result conventions and exact state serialization.
 */
#include "../nacv_tools.h"
#include "source_base/opt_cg.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_reduce.h"
#include "source_lcao/module_lr/utils/lr_util.h"
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <gtest/gtest.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifndef _WIN32
#include <unistd.h>
#endif

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

namespace
{
LR::NacvResult example_result(double gap)
{
    LR::NacvResult result;
    result.gap_rydberg = gap;
    result.ao_connection.create(1, 3);
    result.ao_connection(0, 0) = 0.25;
    result.ao_connection(0, 1) = -0.5;
    return result;
}

ModuleBase::matrix example_force()
{
    ModuleBase::matrix force(1, 3);
    force(0, 0) = -2.0;
    force(0, 1) = 4.0;
    return force;
}
} // namespace

TEST(NacvResult, RequiresCompleteOrbitalWindowWithoutIntegerOverflow)
{
    EXPECT_NO_THROW(LR::validate_nac_window(8, 20, 28));
    EXPECT_THROW(LR::validate_nac_window(7, 20, 28), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_window(8, 19, 28), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_window(0, 28, 28), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_window(28, 0, 28), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_window(-1, 29, 28), std::invalid_argument);
    const int largest = std::numeric_limits<int>::max();
    EXPECT_THROW(LR::validate_nac_window(largest, largest, 28), std::invalid_argument);
}

TEST(NacvResult, RejectsUnavailableResolvedFunctionalNames)
{
    const std::vector<std::string> unsupported = {"PBESOL",
        "REVPBE",
        "WC",
        "BLYP",
        "BP",
        "PW91",
        "HCTH",
        "OLYP",
        "SCAN",
        "SCAN0",
        "SCANL",
        "OPT_ORB",
        "MULLER",
        "POWER",
        "WP22",
        "CWP22",
        "BLYP_LR"};
    for (const std::string &name : unsupported)
    {
        EXPECT_THROW(LR::validate_nac_functionals(name, "pbe"), std::invalid_argument);
        EXPECT_THROW(LR::validate_nac_functionals("PBE", name), std::invalid_argument);
        EXPECT_THROW(LR::validate_nac_functionals("default", name), std::invalid_argument);
        EXPECT_THROW(LR::validate_nac_functionals(name, "default"), std::invalid_argument);
    }
}

TEST(NacvResult, AcceptsIntegerGammaWeightsForBothClosedShellConventions)
{
    for (int spins : {1, 2})
    {
        ModuleBase::matrix weights(spins, 4);
        for (int spin = 0; spin < spins; ++spin)
        {
            weights(spin, 0) = 2.0 / spins;
            weights(spin, 1) = 2.0 / spins;
        }
        EXPECT_NO_THROW(
            LR::validate_nac_occupations(
                weights,
                2,
                4,
                spins));
        weights(0, 2) = std::numeric_limits<double>::epsilon();
        EXPECT_NO_THROW(
            LR::validate_nac_occupations(
                weights,
                2,
                4,
                spins));
    }
}

TEST(NacvResult, RejectsFractionalNonfiniteAndIncompleteOccupations)
{
    ModuleBase::matrix weights(1, 4);
    weights(0, 0) = 2.0;
    weights(0, 1) = 1.9;
    weights(0, 2) = 0.1;
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            2,
            4,
            1),
        std::invalid_argument);
    weights(0, 1) = 2.0;
    weights(0, 2) = 1e-8;
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            2,
            4,
            1),
        std::invalid_argument);
    weights(0, 2) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            2,
            4,
            1),
        std::invalid_argument);
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            2,
            5,
            1),
        std::invalid_argument);
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            0,
            4,
            1),
        std::invalid_argument);
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            2,
            4,
            4),
        std::invalid_argument);
}

TEST(NacvResult, RequiresImplementedGroundAndExcitedResponseKernels)
{
    EXPECT_NO_THROW(LR::validate_nac_functionals("PBE", "pbe"));
    EXPECT_NO_THROW(LR::validate_nac_functionals("PZ", "lda"));
    EXPECT_NO_THROW(LR::validate_nac_functionals("SLAPWPBXPBC", "pbe"));
    EXPECT_NO_THROW(LR::validate_nac_functionals("PBE0", "pbe"));
    EXPECT_NO_THROW(LR::validate_nac_functionals("HF", "hf"));
    EXPECT_NO_THROW(LR::validate_nac_functionals("NOX+NOC", "rpa"));
    EXPECT_THROW(LR::validate_nac_functionals("SCAN", "pbe"), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_functionals("pbe", "scan"), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_functionals("pbe", "bse"), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_functionals("default", "pbe"), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_functionals("MULLER", "hf"), std::invalid_argument);
}

TEST(NacvResult, GapResolutionUsesOperatorResidualsAndRoundoff)
{
    const double roundoff = 16.0 * std::numeric_limits<double>::epsilon() * 2.0;
    EXPECT_DOUBLE_EQ(
        LR::nacv_gap_resolution(
            1.0,
            2.0,
            1e-10,
            2e-10),
        3e-10 + roundoff);
    EXPECT_DOUBLE_EQ(
        LR::nacv_gap_resolution(
            0.0,
            0.0,
            0.0,
            0.0),
        0.0);
    EXPECT_DOUBLE_EQ(
        LR::nacv_gap_resolution(
            -2.0,
            1.0,
            0.0,
            0.0),
        roundoff);
}

TEST(NacvResult, ResolvesSignedGapAndPreservesUnits)
{
    const ModuleBase::matrix force = example_force();
    for (double gap : {2.0, -2.0})
    {
        LR::NacvResult result = example_result(gap);
        LR::assemble_nacv_result(result, force, 1e-8);
        ASSERT_TRUE(result.gap_resolved);
        EXPECT_DOUBLE_EQ(result.numerator_etf(0, 0), 2.0);
        EXPECT_DOUBLE_EQ(result.numerator_raw(0, 0), 2.0 + gap * 0.25);
        EXPECT_DOUBLE_EQ(result.etf(0, 0), 2.0 / gap);
        EXPECT_DOUBLE_EQ(result.raw(0, 0), 2.0 / gap + 0.25);
        EXPECT_DOUBLE_EQ(result.raw(0, 1), -4.0 / gap - 0.5);
        EXPECT_DOUBLE_EQ(result.raw(0, 2), 0.0);
    }
}

TEST(NacvResult, UnresolvedGapRetainsNumeratorsAndClearsDividedVectors)
{
    LR::NacvResult result = example_result(2.0);
    const ModuleBase::matrix force = example_force();
    LR::assemble_nacv_result(result, force, 1e-8);
    for (double gap : {0.0, 1e-8, -1e-8})
    {
        result.gap_rydberg = gap;
        LR::assemble_nacv_result(result, force, 1e-8);
        EXPECT_FALSE(result.gap_resolved);
        EXPECT_DOUBLE_EQ(result.numerator_etf(0, 0), 2.0);
        EXPECT_DOUBLE_EQ(result.numerator_raw(0, 0), 2.0 + gap * 0.25);
        EXPECT_DOUBLE_EQ(result.raw(0, 0), 0.0);
        EXPECT_DOUBLE_EQ(result.etf(0, 1), 0.0);
        EXPECT_DOUBLE_EQ(result.ao_connection(0, 0), 0.25);
    }
}

TEST(NacvResult, ReversePreservesGapWeightedNumerators)
{
    LR::NacvResult result = example_result(2.0);
    const ModuleBase::matrix force = example_force();
    LR::assemble_nacv_result(result, force, 1e-8);
    const LR::NacvResult original = result;
    LR::reverse_nac_pair(result);
    EXPECT_TRUE(result.gap_resolved);
    EXPECT_DOUBLE_EQ(result.gap_rydberg, -original.gap_rydberg);
    EXPECT_DOUBLE_EQ(result.raw(0, 0), -original.raw(0, 0));
    EXPECT_DOUBLE_EQ(result.etf(0, 1), -original.etf(0, 1));
    EXPECT_DOUBLE_EQ(result.ao_connection(0, 0), -original.ao_connection(0, 0));
    EXPECT_DOUBLE_EQ(result.numerator_raw(0, 0), original.numerator_raw(0, 0));
    EXPECT_DOUBLE_EQ(result.numerator_etf(0, 1), original.numerator_etf(0, 1));
    EXPECT_DOUBLE_EQ(result.raw(0, 0) * result.gap_rydberg, result.numerator_raw(0, 0));
    LR::reverse_nac_pair(result);
    EXPECT_DOUBLE_EQ(result.raw(0, 0), original.raw(0, 0));
}

TEST(NacvResult, RejectsNonfiniteConnectionAndForce)
{
    LR::NacvResult result = example_result(2.0);
    ModuleBase::matrix force = example_force();
    result.ao_connection(0, 0) = std::numeric_limits<double>::infinity();
    EXPECT_THROW(LR::assemble_nacv_result(result, force, 1e-8), std::runtime_error);
    result.ao_connection(0, 0) = 0.0;
    force(0, 0) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(LR::assemble_nacv_result(result, force, 1e-8), std::runtime_error);
}

TEST(NacvResult, OutputConvertsNumeratorsToHartreeAndRestoresStream)
{
    LR::NacvResult result = example_result(2.0);
    const ModuleBase::matrix force = example_force();
    LR::assemble_nacv_result(result, force, 1e-8);
    std::ostringstream output;
    output << std::fixed << std::setprecision(3);
    const std::ios::fmtflags flags = output.flags();
    LR::print_nac_result(
        output,
        result,
        1,
        0,
        1,
        "singlet");
    EXPECT_EQ(output.flags(), flags);
    EXPECT_EQ(output.precision(), 3);
    std::istringstream input(output.str());
    std::string tag;
    std::string spin;
    int step = 0;
    int bra = 0;
    int ket = 0;
    bool resolved = false;
    double gap = 0.0;
    input >> tag >> step >> bra >> ket >> spin >> resolved >> gap;
    EXPECT_EQ(tag, "NACV_BEGIN");
    EXPECT_EQ(step, 1);
    EXPECT_EQ(ket, 1);
    EXPECT_TRUE(resolved);
    EXPECT_DOUBLE_EQ(gap, 1.0);
    std::string ignored;
    std::getline(input, ignored);
    std::getline(input, ignored);
    int atom = -1;
    int axis = -1;
    double raw = 0.0;
    double etf = 0.0;
    double connection = 0.0;
    double numerator_raw = 0.0;
    double numerator_etf = 0.0;
    input >> tag >> atom >> axis >> raw >> etf >> connection >> numerator_raw >> numerator_etf;
    ASSERT_TRUE(input.good());
    EXPECT_EQ(tag, "NACV_COMPONENT");
    EXPECT_EQ(atom, 0);
    EXPECT_EQ(axis, 0);
    EXPECT_DOUBLE_EQ(raw, 1.25);
    EXPECT_DOUBLE_EQ(etf, 1.0);
    EXPECT_DOUBLE_EQ(connection, 0.25);
    EXPECT_DOUBLE_EQ(numerator_raw, 1.25);
    EXPECT_DOUBLE_EQ(numerator_etf, 1.0);
}

TEST(NacvResult, OutputReportsStreamFailure)
{
    LR::NacvResult result = example_result(2.0);
    const ModuleBase::matrix force = example_force();
    LR::assemble_nacv_result(result, force, 1e-8);
    std::ostringstream output;
    output.setstate(std::ios::badbit);
    EXPECT_THROW(
        LR::print_nac_result(
            output,
            result,
            1,
            0,
            1,
            "singlet"),
        std::runtime_error);
}

TEST(NacvResult, EigenpairDiagnosticsKeepPhysicalUnitsAndStreamFormat)
{
    const double energies[] = {2.0, 4.0};
    const std::vector<double> residuals{2e-10, 4e-10};
    std::ostringstream output;
    output << std::fixed << std::setprecision(4);
    const std::ios::fmtflags flags = output.flags();
    LR::print_nac_eigenpairs(
        output,
        3,
        "singlet",
        energies,
        residuals);
    EXPECT_EQ(output.flags(), flags);
    EXPECT_EQ(output.precision(), 4);
    std::istringstream input(output.str());
    for (int state = 1; state <= 2; ++state)
    {
        std::string tag;
        std::string spin;
        int step = 0;
        int index = 0;
        double energy = 0.0;
        double residual = 0.0;
        input >> tag >> step >> spin >> index >> energy >> residual;
        ASSERT_FALSE(input.fail());
        EXPECT_EQ(tag, "LR_EIGENPAIR");
        EXPECT_EQ(step, 3);
        EXPECT_EQ(index, state);
        EXPECT_DOUBLE_EQ(energy, state);
        EXPECT_DOUBLE_EQ(residual, state * 1e-10);
    }
}

TEST(NacvResult, AcceptsUnequalIntegerUnrestrictedSpinSpaces)
{
    ModuleBase::matrix weights(2, 8);
    const std::vector<int> occupied{4, 2};
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int band = 0; band < occupied[spin]; ++band)
        {
            weights(spin, band) = 1.0;
        }
    }
    EXPECT_NO_THROW(LR::validate_nac_occupations(weights, occupied, 8));
    EXPECT_THROW(
        LR::validate_nac_occupations(
            weights,
            4,
            8,
            2),
        std::invalid_argument);
    weights(1, 1) = 0.5;
    EXPECT_THROW(LR::validate_nac_occupations(weights, occupied, 8), std::invalid_argument);
    weights(1, 1) = 1.0;
    weights(0, 7) = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(LR::validate_nac_occupations(weights, occupied, 8), std::invalid_argument);
}

TEST(NacvResult, RejectsMissingUnrestrictedSpinWindows)
{
    ModuleBase::matrix weights(2, 8);
    EXPECT_THROW(LR::validate_nac_occupations(weights, std::vector<int>{4}, 8), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_occupations(weights, std::vector<int>{4, 0}, 8), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_occupations(weights, std::vector<int>{4, 8}, 8), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_occupations(weights, std::vector<int>{4, 2}, 7), std::invalid_argument);
}

TEST(NacvResult, AvoidedCrossingMatchesTheExactTwoLevelDerivativeCoupling)
{
    const double slope = 0.5;
    const double coupling = 0.2;
    for (const double position : {-0.8, 0.0, 0.8})
    {
        // H(R) = {{slope*R,coupling},{coupling,-slope*R}} in Ry.
        // Its normalized eigenvectors are (-sin(theta),cos(theta)) and
        // (cos(theta),sin(theta)); contract dH/dR independently of NAC assembly.
        const double diagonal = slope * position;
        const double energy_squared = diagonal * diagonal + coupling * coupling;
        const double energy = std::sqrt(energy_squared);
        const double theta = 0.5 * std::atan2(coupling, diagonal);
        const double bra_first = -std::sin(theta);
        const double bra_second = std::cos(theta);
        const double ket_first = std::cos(theta);
        const double ket_second = std::sin(theta);
        const double hamiltonian_derivative = slope * (bra_first * ket_first - bra_second * ket_second);
        ModuleBase::matrix force(1, 3);
        force(0, 0) = -hamiltonian_derivative;
        LR::NacvResult result;
        result.gap_rydberg = 2.0 * energy;
        result.ao_connection.create(1, 3);
        LR::assemble_nacv_result(result, force, 1e-12);
        const double exact_coupling = -slope * coupling / (2.0 * energy_squared);
        ASSERT_TRUE(result.gap_resolved);
        EXPECT_NEAR(result.raw(0, 0), exact_coupling, 1e-14);
        EXPECT_NEAR(result.etf(0, 0), exact_coupling, 1e-14);
        EXPECT_DOUBLE_EQ(result.raw(0, 1), 0.0);
        EXPECT_DOUBLE_EQ(result.raw(0, 2), 0.0);
        EXPECT_DOUBLE_EQ(result.numerator_raw(0, 0), hamiltonian_derivative);
    }
}

TEST(NacvResult, ElectronicPhaseFlipChangesVectorsAndNumeratorsTogether)
{
    const ModuleBase::matrix force = example_force();
    LR::NacvResult original = example_result(2.0);
    LR::assemble_nacv_result(original, force, 1e-8);
    LR::NacvResult flipped = example_result(2.0);
    flipped.ao_connection *= -1.0;
    const ModuleBase::matrix flipped_force = force * -1.0;
    LR::assemble_nacv_result(flipped, flipped_force, 1e-8);
    EXPECT_DOUBLE_EQ(flipped.gap_rydberg, original.gap_rydberg);
    EXPECT_EQ(flipped.gap_resolved, original.gap_resolved);
    for (int axis = 0; axis < 3; ++axis)
    {
        EXPECT_DOUBLE_EQ(flipped.raw(0, axis), -original.raw(0, axis));
        EXPECT_DOUBLE_EQ(flipped.etf(0, axis), -original.etf(0, axis));
        EXPECT_DOUBLE_EQ(flipped.numerator_raw(0, axis), -original.numerator_raw(0, axis));
        EXPECT_DOUBLE_EQ(flipped.numerator_etf(0, axis), -original.numerator_etf(0, axis));
    }
}

namespace
{
template <typename T> void expect_file_values(const std::string &filename, const std::vector<T> &expected)
{
    std::ifstream input(filename);
    ASSERT_TRUE(input.is_open()) << filename;
    for (const T &reference : expected)
    {
        T actual{};
        input >> actual;
        ASSERT_FALSE(input.fail()) << filename;
        EXPECT_EQ(actual, reference);
    }
    input >> std::ws;
    EXPECT_TRUE(input.eof()) << filename;
}

/** @brief Give each writer test its own filenames in the gtest scratch directory. */
class NacvIoTest : public ::testing::Test
{
  protected:
    std::ostringstream output;
    std::string directory;
    std::string label;
    std::string suffix;
    std::string energy_filename;
    std::vector<std::string> amplitude_filenames;

    void SetUp() override
    {
        directory = ::testing::TempDir();
        label = "nacv_";
        label += ::testing::UnitTest::GetInstance()->current_test_info()->name();
        suffix = "_g2";
        energy_filename = directory + "Excitation_Energy_" + label + suffix + ".dat";
        for (int rank = 0; rank < 3; ++rank)
        {
            const std::string filename =
                directory + "Excitation_Amplitude_" + label + suffix + "_" + std::to_string(rank + 1) + ".dat";
            amplitude_filenames.push_back(filename);
        }
    }

    void TearDown() override
    {
        std::remove(energy_filename.c_str());
        for (const std::string &filename : amplitude_filenames)
        {
            std::remove(filename.c_str());
        }
    }
};
} // namespace

TEST_F(NacvIoTest, RealStatesRoundTripAtDoublePrecision)
{
    const std::vector<double> energies{1.1234567890123457, 2.9876543210987654};
    const std::vector<double> amplitudes{0.12345678901234568,
        -0.9876543210987654,
        1.2345678901234568e-9,
        -0.3333333333333333,
        0.7777777777777778,
        0.0};
    const double *const energy_data = energies.data();
    const double *const amplitude_data = amplitudes.data();
    LR::write_nac_states(
        directory,
        label,
        suffix,
        0,
        energy_data,
        amplitude_data,
        3,
        2,
        output);
    expect_file_values(energy_filename, energies);
    expect_file_values(amplitude_filenames[0], amplitudes);
}

TEST_F(NacvIoTest, ComplexStatesKeepRealAndImaginaryComponents)
{
    const std::vector<double> energies{0.12345678901234568, 0.9876543210987654};
    const std::vector<std::complex<double>> amplitudes{{0.12345678901234568, -0.9876543210987654},
        {-0.3333333333333333, 0.7777777777777778}, {1.2345678901234568e-9, -2.345678901234568e-9}, {0.0, 0.0}};
    const double *const energy_data = energies.data();
    const std::complex<double> *const amplitude_data = amplitudes.data();
    LR::write_nac_states(
        directory,
        label,
        suffix,
        0,
        energy_data,
        amplitude_data,
        2,
        2,
        output);
    expect_file_values(energy_filename, energies);
    expect_file_values(amplitude_filenames[0], amplitudes);
}

TEST_F(NacvIoTest, NonwriterCreatesOnlyItsLocalAmplitudeFile)
{
    const std::vector<double> energies{1.5};
    const std::vector<double> amplitudes{0.25, -0.75};
    const double *const energy_data = energies.data();
    const double *const amplitude_data = amplitudes.data();
    LR::write_nac_states(
        directory,
        label,
        suffix,
        1,
        energy_data,
        amplitude_data,
        2,
        1,
        output);
    std::ifstream energy_file(energy_filename);
    EXPECT_FALSE(energy_file.is_open());
    expect_file_values(amplitude_filenames[1], amplitudes);
    std::ifstream other_rank(amplitude_filenames[0]);
    EXPECT_FALSE(other_rank.is_open());
}

TEST_F(NacvIoTest, NonwriterPreservesTheWriterEnergyFile)
{
    {
        std::ofstream output(energy_filename);
        ASSERT_TRUE(output.is_open());
        output << "0.125\n";
    }
    const std::vector<double> energies{1.5};
    const std::vector<double> amplitudes{-0.5};
    const double *const energy_data = energies.data();
    const double *const amplitude_data = amplitudes.data();
    LR::write_nac_states(
        directory,
        label,
        suffix,
        2,
        energy_data,
        amplitude_data,
        1,
        1,
        output);
    const std::vector<double> original_energies{0.125};
    expect_file_values(energy_filename, original_energies);
    expect_file_values(amplitude_filenames[2], amplitudes);
}

TEST_F(NacvIoTest, EmptyLocalAmplitudeWindowProducesAnEmptyRecord)
{
    const double energies[]{0.5, 1.0};
    const double unused_amplitude = 7.0;
    LR::write_nac_states(
        directory,
        label,
        suffix,
        1,
        energies,
        &unused_amplitude,
        0,
        2,
        output);
    const std::vector<double> empty;
    expect_file_values(amplitude_filenames[1], empty);
    std::ifstream energy_file(energy_filename);
    EXPECT_FALSE(energy_file.is_open());
}

TEST(XcKernelNames, CanonicalizesOnlyActualBuiltinAliases)
{
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("PZ"), "lda");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("SLAPZNOGXNOGC"), "lda");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("SLAPWPBXPBC"), "pbe");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("NOX+NOC"), "none");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("RPA"), "rpa");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("OPT_ORB"), "opt_orb");
    EXPECT_EQ(LR::NacvXC::canonical_xc_name("Unknown"), "unknown");
}

TEST(XcKernelNames, PreservesTheExistingCanonicalKernelCatalogue)
{
    EXPECT_TRUE(LR_Util::has_local_xc("lda"));
    EXPECT_TRUE(LR_Util::has_local_xc("pwlda"));
    EXPECT_TRUE(LR_Util::has_local_xc("pbe"));
    EXPECT_EQ(LR_Util::hybrid_xc_list().size(), 9u);
    for (const std::string &name : LR_Util::hybrid_xc_list())
    {
        EXPECT_EQ(LR_Util::has_local_xc(name), name != "hf");
    }
    for (const std::string &name : {"rpa",
             "none",
             "pz",
             "pbesol",
             "scan",
             "blyp_lr",
             "muller"})
    {
        EXPECT_FALSE(LR_Util::has_local_xc(name));
    }
}

TEST(NacvSurfaces, ChecksTheResolvedRootCountRatherThanTheInputDefault)
{
    const std::vector<int> boundary{0,
        4,
        1,
        4};
    EXPECT_NO_THROW(LR::validate_nac_surfaces(boundary, 4));
    EXPECT_THROW(LR::validate_nac_surfaces(boundary, 3), std::invalid_argument);
    EXPECT_THROW(LR::validate_nac_surfaces(boundary, 0), std::invalid_argument);
    const std::vector<int> disabled;
    EXPECT_NO_THROW(LR::validate_nac_surfaces(disabled, 0));
}

TEST(NacvIoMpi, ReportsEveryCompletedRankFileToTheWriterLog)
{
    const auto domain = ModuleBase::world_comm_domain();
    const int rank = domain.rank();
    const int processes = domain.size();
    const std::string directory = ::testing::TempDir();
    const std::string label = "collective_nacv_io";
    const std::string suffix = "_g1";
    const double energies[]{0.125};
    const double amplitudes[]{5e-11};
    std::ostringstream output;
    LR::write_nac_states(
        directory,
        label,
        suffix,
        rank,
        energies,
        amplitudes,
        1,
        1,
        output);
    const std::string energy = directory + "Excitation_Energy_" + label + suffix + ".dat";
    if (rank == 0)
    {
        EXPECT_NE(output.str().find("Write file: " + energy), std::string::npos);
        for (int process = 0; process < processes; ++process)
        {
            const int output_rank = process + 1;
            const std::string amplitude =
                directory + "Excitation_Amplitude_" + label + suffix + "_" + std::to_string(output_rank) + ".dat";
            EXPECT_NE(output.str().find("Write file: " + amplitude), std::string::npos);
        }
        std::remove(energy.c_str());
    }
    else
    {
        EXPECT_TRUE(output.str().empty());
    }
    const int output_rank = rank + 1;
    const std::string amplitude =
        directory + "Excitation_Amplitude_" + label + suffix + "_" + std::to_string(output_rank) + ".dat";
    std::remove(amplitude.c_str());
}

TEST_F(NacvIoTest, RetainsValuesBelowAndAtTheLegacyFilterThreshold)
{
    const std::vector<double> energies{5e-11, 1e-10};
    const std::vector<double> amplitudes{5e-11,
        -5e-11,
        1e-10,
        -1e-10,
        1e-30,
        0.0};
    const double *energy_data = energies.data();
    const double *amplitude_data = amplitudes.data();
    LR::write_nac_states(
        directory,
        label,
        suffix,
        0,
        energy_data,
        amplitude_data,
        3,
        2,
        output);
    expect_file_values(energy_filename, energies);
    expect_file_values(amplitude_filenames[0], amplitudes);
    EXPECT_NE(output.str().find("Write file: " + energy_filename), std::string::npos);
    EXPECT_NE(output.str().find("Write file: " + amplitude_filenames[0]), std::string::npos);
}

TEST_F(NacvIoTest, RetainsSmallComplexComponentsIndependently)
{
    const double energies[]{0.125};
    const std::vector<std::complex<double>> amplitudes{{5e-11, -1e-10}, {1e-30, -5e-11}};
    const std::complex<double> *data = amplitudes.data();
    LR::write_nac_states(
        directory,
        label,
        suffix,
        0,
        energies,
        data,
        2,
        1,
        output);
    expect_file_values(amplitude_filenames[0], amplitudes);
}

TEST_F(NacvIoTest, RejectsAnOpenFailureWithoutReportingSuccess)
{
    {
        std::ofstream existing(energy_filename);
        existing << "sentinel\n";
    }
    const std::string invalid_directory = energy_filename + "/";
    const double energies[]{0.125};
    const double amplitudes[]{0.25};
    EXPECT_THROW(
        LR::write_nac_states(
            invalid_directory,
            label,
            suffix,
            0,
            energies,
            amplitudes,
            1,
            1,
            output),
        std::runtime_error);
    EXPECT_TRUE(output.str().empty());
}

#ifndef _WIN32
TEST(NacvIoMpi, ReportsOneRankWriteFailureOnEveryRank)
{
    const auto domain = ModuleBase::world_comm_domain();
    const int processes = domain.size();
    if (processes < 2)
    {
        GTEST_SKIP() << "Requires a nonwriter rank";
    }
    const int rank = domain.rank();
    const std::string directory = ::testing::TempDir();
    const std::string label = "collective_nacv_failure";
    const std::string suffix = "_g1";
    const std::string energy = directory + "Excitation_Energy_" + label + suffix + ".dat";
    const int output_rank = rank + 1;
    const std::string amplitude =
        directory + "Excitation_Amplitude_" + label + suffix + "_" + std::to_string(output_rank) + ".dat";
    const bool failing_rank = output_rank == processes;
    int setup_failure = 0;
    if (failing_rank)
    {
        setup_failure = ::symlink("/dev/full", amplitude.c_str()) != 0;
    }
    ASSERT_EQ(domain.max(setup_failure), 0);

    const double energies[]{0.125};
    const double amplitudes[]{0.25};
    std::ostringstream output;
    EXPECT_THROW(
        LR::write_nac_states(
            directory,
            label,
            suffix,
            rank,
            energies,
            amplitudes,
            1,
            1,
            output),
        std::runtime_error);
    EXPECT_TRUE(output.str().empty());
    std::remove(amplitude.c_str());
    if (rank == 0)
    {
        std::remove(energy.c_str());
    }
}

TEST_F(NacvIoTest, DetectsBufferedWriteFailureBeforeReportingSuccess)
{
    ASSERT_EQ(::symlink("/dev/full", energy_filename.c_str()), 0);
    const double energies[]{0.125};
    const double amplitudes[]{0.25};
    EXPECT_THROW(
        LR::write_nac_states(
            directory,
            label,
            suffix,
            0,
            energies,
            amplitudes,
            1,
            1,
            output),
        std::runtime_error);
    EXPECT_TRUE(output.str().empty());
}
#endif
