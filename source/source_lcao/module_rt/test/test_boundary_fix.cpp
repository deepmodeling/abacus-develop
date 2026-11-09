#include "boundary_test.h"

namespace
{
class BoundaryFixTest : public rt_test::BoundaryTest
{
    void SetUp() override
    {
        BoundaryTest::SetUp();
        types[0].boundary_shift[0].set(-1, 1, -1);
        types[0].boundary_shift[1].set(1, 0, 1);
    }
};

TEST_P(BoundaryFixTest, BlochPhaseCovariance)
{
    const rt_test::Vector no_field;
    check_correction(no_field, rt_test::analytic_matrix);
}

TEST_P(BoundaryFixTest, HybridUsesSavedFieldInNonorthogonalCell)
{
    const rt_test::Vector saved_field(0.07, -0.11, 0.05);
    check_correction(saved_field, rt_test::analytic_matrix);
}

TEST_P(BoundaryFixTest, NoBoundaryShiftPreservesSavedState)
{
    for (auto &shift : types[0].boundary_shift)
    {
        shift.set(0, 0, 0);
    }
    const rt_test::Vector no_field;
    check_correction(no_field, rt_test::analytic_matrix);
}

INSTANTIATE_TEST_SUITE_P(MatrixLayouts, BoundaryFixTest, testing::Bool());
} // namespace
