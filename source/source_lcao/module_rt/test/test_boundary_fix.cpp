#include "boundary_test.h"

namespace
{
class BoundaryFixTest : public rt_test::BoundaryTest
{
};

TEST_P(BoundaryFixTest, SavedMatricesAndCoefficientsFollowBoundaryPhase)
{
    types[0].boundary_shift[0].set(-1, 0, 0);
    types[0].boundary_shift[1].set(1, 0, 0);
    const rt_test::Vector no_field;
    const rt_test::Vector saved_A(0.07, -0.11, 0.05);
    ASSERT_NO_FATAL_FAILURE(expect_phase_correction(no_field));
    ASSERT_NO_FATAL_FAILURE(expect_phase_correction(saved_A));
}

INSTANTIATE_TEST_SUITE_P(MatrixLayouts, BoundaryFixTest, testing::Bool());
} // namespace
