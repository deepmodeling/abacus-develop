#include "diago_ppcg_test_support.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace hsolver
{

class DiagoPPCGRefinementTest : public ::testing::Test
{
  protected:
    void check_eigenvalue_only_refinement(const double overlap)
    {
        using ppcg_test::T;
        const int dimension = 2;
        const double exact_value = 16.0;
        const double inverse_norm = 1.0 / std::sqrt(overlap);
        std::vector<T> vectors = {T(inverse_norm), T(0)};
        std::vector<double> values(1, 0.0);
        const std::vector<double> thresholds(1, 1e-14);
        const std::vector<double> preconditioner(dimension, 1.0);
        const auto hop = [overlap](const T *input, T *output, int ld, int count)
        {
            for (int band = 0; band < count; ++band)
            {
                output[band * ld] = 16.0 * overlap * input[band * ld];
                output[band * ld + 1] = 32.0 * overlap * input[band * ld + 1];
            }
        };
        const auto sop = [overlap](const T *input, T *output, int ld, int count)
        {
            for (int entry = 0; entry < ld * count; ++entry)
            {
                output[entry] = overlap * input[entry];
            }
        };
        const ppcg_test::CallbackHSOperator<T> op(hop, sop);
        DiagoPPCG<T, base_device::DEVICE_CPU> solver(1e-14, 0, 1, 1);
        ppcg_test::diagonalize(
            solver, hop, sop, dimension, 1, dimension, vectors.data(), values.data(), thresholds,
            preconditioner.data());
        values[0] = exact_value;
        const double infinity = std::numeric_limits<double>::infinity();
        for (int step = 0; step < 4; ++step)
        {
            values[0] = std::nextafter(values[0], infinity);
        }
        solver.compute_eigenpair_residual(values.data());
        const std::vector<T> original = vectors;
        const T original_residual_value = (exact_value - values[0]) * overlap * vectors[0];
        const double original_residual = std::abs(original_residual_value);
        ASSERT_GT(original_residual, thresholds[0]);
        const std::vector<int> columns(1, 0);
        const bool accepted =
            solver.refine_stalled_residual(op, vectors.data(), values.data(), thresholds, columns, -1);
        EXPECT_TRUE(accepted);
        EXPECT_EQ(vectors, original);
        EXPECT_DOUBLE_EQ(values[0], exact_value);
        const T residual_value = (exact_value - values[0]) * overlap * vectors[0];
        const double residual = std::abs(residual_value);
        EXPECT_LT(residual, thresholds[0]);
    }

    void check_refinement(
        const bool generalized,
        const bool exact,
        const int invalid_application,
        const double requested_threshold,
        const int band_hint)
    {
        using ppcg_test::T;
        const int dimension = 4;
        const int leading_dimension = 6;
        const int bands = 2;
        const int size = leading_dimension * bands;
        const std::vector<double> spectrum = {0.5, 2.0, 5.0, 8.0};
        std::vector<T> h(dimension * dimension, T(0));
        std::vector<T> s(dimension * dimension, T(0));
        std::vector<double> overlap(dimension, 1.0);
        for (int row = 0; row < dimension; ++row)
        {
            if (generalized)
            {
                overlap[row] = 0.5 + row;
            }
            h[row + row * dimension] = spectrum[row] * overlap[row];
            s[row + row * dimension] = overlap[row];
        }
        std::vector<T> vectors(size, T(0));
        vectors[0] = 1.0 / std::sqrt(overlap[0]);
        const double contamination = exact ? 0.0 : 0.1;
        vectors[leading_dimension + 1] = std::sqrt(1.0 - contamination * contamination) / std::sqrt(overlap[1]);
        vectors[leading_dimension + 2] = T(0, contamination / std::sqrt(overlap[2]));
        int applications_until_failure = 0;
        const auto hop = [&](const T *input, T *output, int ld, int count)
        {
            if (applications_until_failure > 0 && --applications_until_failure == 0)
            {
                const double invalid = std::numeric_limits<double>::quiet_NaN();
                std::fill(output, output + ld * count, T(invalid));
                return;
            }
            ppcg_test::dense_h_multiply(h.data(), dimension, input, output, ld, count);
        };
        const auto sop = [&](const T *input, T *output, int ld, int count)
        { ppcg_test::dense_h_multiply(s.data(), dimension, input, output, ld, count); };
        const ppcg_test::CallbackHSOperator<T> op(hop, sop);
        const std::vector<double> thresholds(bands, requested_threshold);
        const std::vector<double> preconditioner(dimension, 1.0);
        std::vector<double> values(bands);
        DiagoPPCG<T, base_device::DEVICE_CPU> solver(1e-14, 0, 2, 1);
        ppcg_test::diagonalize(
            solver, hop, sop, leading_dimension, bands, dimension, vectors.data(), values.data(), thresholds,
            preconditioner.data());
        const std::vector<T> original = vectors;
        const std::vector<double> original_values = values;
        const double original_residual = solver.max_residual();
        const std::vector<int> columns = {0, 1};
        applications_until_failure = invalid_application;
        const bool accepted =
            solver.refine_stalled_residual(op, vectors.data(), values.data(), thresholds, columns, band_hint);
        const bool expected_acceptance = !exact && invalid_application == 0;
        EXPECT_EQ(accepted, expected_acceptance);
        if (!expected_acceptance)
        {
            EXPECT_EQ(vectors, original);
            EXPECT_EQ(values, original_values);
        }
        for (int row = 0; row < leading_dimension; ++row)
        {
            EXPECT_EQ(vectors[row], original[row]);
        }
        double maximum_residual = 0.0;
        for (int band = 0; band < bands; ++band)
        {
            double residual_squared = 0.0;
            for (int row = 0; row < dimension; ++row)
            {
                const T x = vectors[row + band * leading_dimension];
                const T residual = (spectrum[row] - values[band]) * overlap[row] * x;
                residual_squared += std::norm(residual);
            }
            maximum_residual = std::max(maximum_residual, std::sqrt(residual_squared));
            for (int other = 0; other < bands; ++other)
            {
                T inner_product(0);
                for (int row = 0; row < dimension; ++row)
                {
                    inner_product += std::conj(vectors[row + band * leading_dimension]) * overlap[row] *
                                     vectors[row + other * leading_dimension];
                }
                const T expected = band == other ? T(1) : T(0);
                const double error = std::abs(inner_product - expected);
                EXPECT_LT(error, 2e-14);
            }
            for (int row = dimension; row < leading_dimension; ++row)
            {
                EXPECT_EQ(vectors[row + band * leading_dimension], T(0));
            }
        }
        if (accepted)
        {
            EXPECT_LT(maximum_residual, original_residual);
            EXPECT_LE(values[1], original_values[1]);
            EXPECT_GE(values[1], spectrum[1]);
        }
    }
};

TEST_F(DiagoPPCGRefinementTest, CorrectsEigenvalueRoundoffWithoutMovingExactStates)
{
    check_eigenvalue_only_refinement(1.0);
    check_eigenvalue_only_refinement(2.0);
}

TEST_F(DiagoPPCGRefinementTest, PreservesStandardGeometry) { check_refinement(false, false, 0, 1e-14, -1); }

TEST_F(DiagoPPCGRefinementTest, PreservesGeneralizedGeometry) { check_refinement(true, false, 0, 1e-14, -1); }

TEST_F(DiagoPPCGRefinementTest, LeavesConvergedPairsUntouched) { check_refinement(true, true, 0, 1e-14, -1); }

TEST_F(DiagoPPCGRefinementTest, RejectsNonfiniteDirection) { check_refinement(true, false, 1, 1e-14, -1); }

TEST_F(DiagoPPCGRefinementTest, RestoresAfterNonfiniteCandidate) { check_refinement(true, false, 2, 1e-14, -1); }

TEST_F(DiagoPPCGRefinementTest, RestoresAfterNonfiniteNormalizedCandidate)
{
    check_refinement(true, false, 3, 1e-14, -1);
}

TEST_F(DiagoPPCGRefinementTest, RefinesConvergedNeighborWithoutLosingOrthogonality)
{
    check_refinement(true, false, 0, 0.5, 1);
}

} // namespace hsolver
