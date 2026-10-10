#include <gtest/gtest.h>
#include "../lr_amp.h"
#include "source_base/parallel_global.h"
#include <complex>
#include <sstream>

namespace { int test_rank = 0; }

TEST(GradientAmplitudes, InitialRootSeedsTheReference)
{
    const double roots[] = {1.0, 0.0, 0.0, 1.0};
    int target = -1;
    std::vector<double> previous;
    std::ostringstream log;
    LR::follow_root(roots, 2, 2, 1, target, previous, log);
    EXPECT_EQ(target, 1);
    EXPECT_EQ(previous, (std::vector<double>{0.0, 1.0}));
}

TEST(GradientAmplitudes, FollowsACrossingIndependentOfSign)
{
    const double roots[] = {0.0, 1.0, -1.0, 0.0};
    int target = 0;
    std::vector<double> previous{1.0, 0.0};
    std::ostringstream log;
    LR::follow_root(roots, 2, 2, 0, target, previous, log);
    EXPECT_EQ(target, 1);
    EXPECT_EQ(previous, (std::vector<double>{-1.0, 0.0}));
    EXPECT_NE(log.str().find("root moved"), std::string::npos);
}

TEST(GradientAmplitudes, ComplexPhaseDoesNotChangeTheFollowedRoot)
{
    using Complex = std::complex<double>;
    const Complex roots[] = {Complex(0.0), Complex(1.0), Complex(0.0, 1.0), Complex(0.0)};
    int target = 0;
    std::vector<Complex> previous{Complex(1.0), Complex(0.0)};
    std::ostringstream log;
    LR::follow_root(roots, 2, 2, 0, target, previous, log);
    EXPECT_EQ(target, 1);
    EXPECT_EQ(previous[0], Complex(0.0, 1.0));
}

TEST(GradientAmplitudes, OpenShellDensitiesReadTheEntireDownSpinKBlock)
{
    // Two roots, two k points; up has two pairs per k and down has one.
    const double roots[] = {10, 11, 12, 13, 20, 21, 30, 31, 32, 33, 40, 41};
    const int first_down = LR::electron_hole_offset(0, 2, 2, 1, 1, true);
    const int second_down = LR::electron_hole_offset(1, 2, 2, 1, 1, true);
    const int second_up = LR::electron_hole_offset(1, 2, 2, 1, 0, true);
    EXPECT_DOUBLE_EQ(roots[first_down], 20);
    EXPECT_DOUBLE_EQ(roots[first_down + 1], 21);
    EXPECT_DOUBLE_EQ(roots[second_down], 40);
    EXPECT_DOUBLE_EQ(roots[second_down + 1], 41);
    EXPECT_DOUBLE_EQ(roots[second_up], 30);
    const int closed = LR::electron_hole_offset(1, 2, 2, 1, 1, false);
    EXPECT_EQ(closed, 2);
    const int gamma = LR::electron_hole_offset(1, 1, 2, 1, 1, true);
    EXPECT_EQ(gamma, 5);
}

TEST(GradientAmplitudes, FollowsTheSelectedJTBranchAfterSplitting)
{
    const double degenerate_roots[] = {1.0, 0.0, 0.0, 1.0};
    const std::vector<int> group{0, 1};
    const std::vector<double> mixing{0.0, 1.0};
    std::vector<double> previous{1.0, 0.0};
    LR::save_mixed_root(degenerate_roots, 2, group, mixing, previous);
    // The selected second component becomes the first root at the next geometry.
    const double split_roots[] = {0.0, 1.0, 1.0, 0.0};
    int target = 0;
    std::ostringstream log;
    LR::follow_root(split_roots, 2, 2, 0, target, previous, log);
    EXPECT_EQ(target, 0);
    EXPECT_EQ(previous, (std::vector<double>{0.0, 1.0}));
}

TEST(GradientAmplitudes, MixedReferenceHasGlobalUnitNorm)
{
    using Complex = std::complex<double>;
    const Complex roots[] = {Complex(0, 1), Complex(0), Complex(0), Complex(1)};
    const std::vector<int> group{0, 1};
    const std::vector<double> mixing{0.6, 0.8};
    std::vector<Complex> previous;
    LR::save_mixed_root(roots, 2, group, mixing, previous);
    ASSERT_EQ(previous.size(), 2);
    double norm_squared = std::norm(previous[0]) + std::norm(previous[1]);
    Parallel_Reduce::reduce_all(norm_squared);
    EXPECT_NEAR(norm_squared, 1.0, 1e-14);
    EXPECT_NEAR(previous[0].imag() / previous[1].real(), 0.75, 1e-14);
}

TEST(GradientAmplitudes, KeepsJTMixtureAfterIndividualOrbitalSignChange)
{
    const double old_roots[] = {1, 0, 0, 1};
    const std::vector<int> group{0, 1};
    const std::vector<double> mixing{0.6, 0.8};
    std::vector<double> previous;
    LR::save_mixed_root(old_roots, 2, group, mixing, previous);
    // The sign-flipped second virtual orbital changes that reference component.
    // The full projection is tested against its formula in test_state_ovlp.cpp.
    previous[1] = -previous[1];
    const double current_roots[] = {0.6, -0.8, 0.8, 0.6};
    int target = 0;
    std::ostringstream log;
    LR::follow_root(current_roots, 2, 2, 0, target, previous, log);
    EXPECT_EQ(target, 0);
}

TEST(GradientAmplitudes, EmptyLocalRanksParticipateInRootSelection)
{
    const int local_size = test_rank == 0 ? 2 : 0;
    const double first[] = {1, 0, 0, 1};
    const double second[] = {0, 1, 1, 0};
    int target = 0;
    std::vector<double> previous;
    std::ostringstream log;
    LR::follow_root(first, local_size, 2, 0, target, previous, log);
    LR::follow_root(second, local_size, 2, 0, target, previous, log);
    EXPECT_EQ(target, 1);
    EXPECT_EQ(previous.size(), static_cast<std::size_t>(local_size));
}

int main(int argc, char** argv)
{
    int processes = 1;
    int threads = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, processes, threads, rank);
    test_rank = rank;
#ifdef __MPI
    // This focused test has no pool/grid setup; the cleanup wrapper needs null handles.
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
#endif
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    Parallel_Global::finalize_mpi();
#endif
    return result;
}
