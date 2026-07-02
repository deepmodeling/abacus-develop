#include "source_lcao/module_ri/sternheimer_dft_zero_order.h"

#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

TEST(SternheimerDFTZeroOrder, ExactEigenvalueReferencePasses)
{
    ModuleRI::SternheimerFDZeroOrderStates fd_states;
    fd_states.eigenvalues = {-0.5, 0.125, 0.5};
    fd_states.residual_norms = {1.0e-13, 2.0e-13, 3.0e-13};

    const std::vector<double> dft_eigenvalues = {-0.5, 0.125, 0.5};
    const auto comparison
        = ModuleRI::compare_sternheimer_fd_zero_order_to_dft(fd_states, dft_eigenvalues, 1.0e-10);

    ASSERT_EQ(comparison.bands.size(), 3);
    EXPECT_TRUE(comparison.all_eigenvalues_within_tolerance);
    EXPECT_DOUBLE_EQ(comparison.max_abs_fd_minus_dft, 0.0);
    EXPECT_EQ(comparison.bands[1].band_index, 1);
    EXPECT_DOUBLE_EQ(comparison.bands[1].fd_eigenvalue, 0.125);
    EXPECT_DOUBLE_EQ(comparison.bands[1].dft_eigenvalue, 0.125);
    EXPECT_DOUBLE_EQ(comparison.bands[1].fd_minus_dft, 0.0);
    EXPECT_DOUBLE_EQ(comparison.bands[1].fd_residual_norm, 2.0e-13);
    EXPECT_TRUE(comparison.bands[1].eigenvalue_within_tolerance);
}

TEST(SternheimerDFTZeroOrder, EigenvalueMismatchIsReportedPerBand)
{
    ModuleRI::SternheimerFDZeroOrderStates fd_states;
    fd_states.eigenvalues = {-0.5, 0.125, 0.5};
    fd_states.residual_norms = {1.0e-13, 2.0e-13, 3.0e-13};

    const std::vector<double> dft_eigenvalues = {-0.5, 0.1, 0.5};
    const auto comparison
        = ModuleRI::compare_sternheimer_fd_zero_order_to_dft(fd_states, dft_eigenvalues, 1.0e-3);

    ASSERT_EQ(comparison.bands.size(), 3);
    EXPECT_FALSE(comparison.all_eigenvalues_within_tolerance);
    EXPECT_NEAR(comparison.max_abs_fd_minus_dft, 0.025, 1.0e-15);
    EXPECT_EQ(comparison.bands[1].band_index, 1);
    EXPECT_DOUBLE_EQ(comparison.bands[1].fd_eigenvalue, 0.125);
    EXPECT_DOUBLE_EQ(comparison.bands[1].dft_eigenvalue, 0.1);
    EXPECT_NEAR(comparison.bands[1].fd_minus_dft, 0.025, 1.0e-15);
    EXPECT_FALSE(comparison.bands[1].eigenvalue_within_tolerance);
}

TEST(SternheimerDFTZeroOrder, RejectsInconsistentInputs)
{
    ModuleRI::SternheimerFDZeroOrderStates fd_states;
    fd_states.eigenvalues = {-0.5, 0.125};
    fd_states.residual_norms = {1.0e-13};

    EXPECT_THROW(ModuleRI::compare_sternheimer_fd_zero_order_to_dft(fd_states, {-0.5, 0.125}, 1.0e-10),
                 std::invalid_argument);

    fd_states.residual_norms = {1.0e-13, 2.0e-13};
    EXPECT_THROW(ModuleRI::compare_sternheimer_fd_zero_order_to_dft(fd_states, {-0.5}, 1.0e-10),
                 std::invalid_argument);
    EXPECT_THROW(ModuleRI::compare_sternheimer_fd_zero_order_to_dft(fd_states, {-0.5, 0.125}, -1.0),
                 std::invalid_argument);
}
