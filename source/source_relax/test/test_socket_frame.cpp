#include "../socket_frame.h"

#include "gtest/gtest.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace
{
using SocketFrame::CellValidation;
using SocketFrame::Matrix9;
using SocketFrame::VirialConversion;
using SocketFrame::checked_position_count;
using SocketFrame::make_ipi_virial;
using SocketFrame::transpose_matrix9;
using SocketFrame::validate_ipi_cell;
using SocketFrame::validate_positions;

const double EPSILON = std::numeric_limits<double>::epsilon();

CellValidation validate_with_driver_thresholds(const Matrix9& cell, const Matrix9& inverse)
{
    return validate_ipi_cell(cell, inverse, 1.0e12, 64.0 * EPSILON, 64.0);
}

void expect_matrix_near(const Matrix9& expected, const Matrix9& actual, double tolerance)
{
    for (std::size_t index = 0; index < expected.size(); ++index)
    {
        EXPECT_NEAR(expected[index], actual[index], tolerance) << "matrix index " << index;
    }
}
} // namespace

TEST(SocketFrameTest, TransposeKeepsAllNineUniqueEntries)
{
    Matrix9 in = {{1, 2, 3, 4, 5, 6, 7, 8, 9}};
    Matrix9 expected = {{1, 4, 7, 2, 5, 8, 3, 6, 9}};
    EXPECT_EQ(expected, transpose_matrix9(in));
}

TEST(SocketFrameTest, VirialUsesPositiveHalfVolumeAndWireTranspose)
{
    Matrix9 stress = {{1, 2, 3, 2, 5, 6, 3, 6, 9}};
    VirialConversion out = make_ipi_virial(stress, 4.0, 1e-12, 1e-12);
    Matrix9 expected = {{2, 4, 6, 4, 10, 12, 6, 12, 18}};
    ASSERT_TRUE(out.ok) << out.message;
    EXPECT_EQ(expected, out.wire_virial_hartree);
}

TEST(SocketFrameTest, RightHandedTriclinicCellReturnsKnownInverse)
{
    const Matrix9 cell = {{2.0, 1.0, 0.0, 0.0, 3.0, 1.0, 0.0, 0.0, 4.0}};
    const Matrix9 inverse = {{0.5, -1.0 / 6.0, 1.0 / 24.0,
                              0.0, 1.0 / 3.0, -1.0 / 12.0,
                              0.0, 0.0, 0.25}};

    const CellValidation out = validate_with_driver_thresholds(cell, inverse);

    ASSERT_TRUE(out.ok) << out.message;
    EXPECT_DOUBLE_EQ(24.0, out.determinant_bohr3);
    EXPECT_NEAR(0.0, out.inverse_residual, 16.0 * EPSILON);
    expect_matrix_near(inverse, out.computed_inverse_wire_bohr_inv, 16.0 * EPSILON);
}

TEST(SocketFrameTest, AseTriclinicInverseWireLayoutIsAcceptedAndRecomputedFromCell)
{
    // ASE stores row lattice vectors A in Angstrom, sends H = A^T / Bohr,
    // and sends pinv(A) * Bohr as the inverse field. For this nonsingular
    // cell that received field is inv(H)^T, not inv(H).
    const double bohr_angstrom = 0.5291772105638411;
    const Matrix9 cell_wire = {{5.0 / bohr_angstrom, 0.5 / bohr_angstrom, 0.25 / bohr_angstrom,
                                0.0, 4.0 / bohr_angstrom, 0.75 / bohr_angstrom,
                                0.0, 0.0, 3.0 / bohr_angstrom}};
    const Matrix9 ase_inverse_wire = {{bohr_angstrom / 5.0, 0.0, 0.0,
                                       -bohr_angstrom / 40.0, bohr_angstrom / 4.0, 0.0,
                                       -bohr_angstrom / 96.0, -bohr_angstrom / 16.0,
                                       bohr_angstrom / 3.0}};
    const Matrix9 inverse_computed_from_cell = {{bohr_angstrom / 5.0,
                                                  -bohr_angstrom / 40.0,
                                                  -bohr_angstrom / 96.0,
                                                  0.0,
                                                  bohr_angstrom / 4.0,
                                                  -bohr_angstrom / 16.0,
                                                  0.0,
                                                  0.0,
                                                  bohr_angstrom / 3.0}};

    const CellValidation out = validate_with_driver_thresholds(cell_wire, ase_inverse_wire);

    ASSERT_TRUE(out.ok) << out.message;
    EXPECT_NEAR(0.0, out.inverse_residual, 16.0 * EPSILON);
    expect_matrix_near(inverse_computed_from_cell,
                       out.computed_inverse_wire_bohr_inv,
                       16.0 * EPSILON);
}

TEST(SocketFrameTest, RotatedDiagonalTracksRightSingularVectorsAndInverseOrder)
{
    // Hand-multiplied U diag(5, 2, 0.5) V^T, with rational plane rotations.
    const Matrix9 cell = {{4.0, -0.72, -0.96,
                           3.0, 0.96, 1.28,
                           0.0, -0.4, 0.3}};
    const Matrix9 inverse = {{0.16, 0.12, 0.0,
                              -0.18, 0.24, -1.6,
                              -0.24, 0.32, 1.2}};

    const CellValidation out = validate_with_driver_thresholds(cell, inverse);

    ASSERT_TRUE(out.ok) << out.message;
    EXPECT_NEAR(5.0, out.determinant_bohr3, 64.0 * EPSILON);
    EXPECT_NEAR(10.0, out.condition_number_2, 256.0 * EPSILON);
    expect_matrix_near(inverse, out.computed_inverse_wire_bohr_inv, 64.0 * EPSILON);
}

TEST(SocketFrameTest, InconsistentReceivedInverseIsRejected)
{
    const Matrix9 cell = {{2.0, 1.0, 0.0, 0.0, 3.0, 1.0, 0.0, 0.0, 4.0}};
    // This is neither inv(cell) nor inv(cell)^T, so both supported wire
    // layouts must reject it.
    const Matrix9 wrong_inverse = {{0.6, -1.0 / 6.0, 1.0 / 24.0,
                                     0.0, 1.0 / 3.0, -1.0 / 12.0,
                                     0.0, 0.0, 0.25}};

    const CellValidation out = validate_with_driver_thresholds(cell, wrong_inverse);

    EXPECT_FALSE(out.ok);
    EXPECT_NE(std::string::npos, out.message.find("inverse"));
    EXPECT_GT(out.inverse_residual, 0.1);
}

TEST(SocketFrameTest, ReceivedInverseResidualUsesConditionScaledRelativeTolerance)
{
    const Matrix9 cell = {{1.0, 0.0, 0.0,
                           0.0, 1.0, 0.0,
                           0.0, 0.0, 1.0e-6}};
    Matrix9 accepted_inverse = {{1.0 + 1.0e-8, 0.0, 0.0,
                                 0.0, 1.0, 0.0,
                                 0.0, 0.0, 1.0e6}};
    Matrix9 rejected_inverse = accepted_inverse;
    rejected_inverse[0] = 1.0 + 2.0e-8;

    const CellValidation accepted
        = validate_ipi_cell(cell, accepted_inverse, 1.0e12, 0.0, 64.0);
    const CellValidation rejected
        = validate_ipi_cell(cell, rejected_inverse, 1.0e12, 0.0, 64.0);

    ASSERT_TRUE(accepted.ok) << accepted.message;
    EXPECT_DOUBLE_EQ(1.0e6, accepted.condition_number_2);
    EXPECT_NEAR(1.0e-8, accepted.inverse_residual, EPSILON);
    EXPECT_FALSE(rejected.ok);
    EXPECT_NE(std::string::npos, rejected.message.find("inverse"));
    EXPECT_NEAR(2.0e-8, rejected.inverse_residual, EPSILON);
}

TEST(SocketFrameTest, NegativeAndZeroDeterminantsAreRejected)
{
    const Matrix9 identity = {{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}};
    const Matrix9 left_handed = {{-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}};
    const Matrix9 singular = {{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0}};

    EXPECT_FALSE(validate_with_driver_thresholds(left_handed, identity).ok);
    EXPECT_FALSE(validate_with_driver_thresholds(singular, identity).ok);
}

TEST(SocketFrameTest, NonrepresentablePositiveCellVolumeIsRejected)
{
    const Matrix9 huge_cell = {{1.0e200, 0.0, 0.0,
                                0.0, 1.0e200, 0.0,
                                0.0, 0.0, 1.0e200}};
    const Matrix9 tiny_inverse = {{1.0e-200, 0.0, 0.0,
                                   0.0, 1.0e-200, 0.0,
                                   0.0, 0.0, 1.0e-200}};

    const CellValidation out = validate_with_driver_thresholds(huge_cell, tiny_inverse);

    EXPECT_FALSE(out.ok);
    EXPECT_NE(std::string::npos, out.message.find("determinant"));
}

TEST(SocketFrameTest, UnderflowedCellVolumeIsRejectedAsZeroDeterminant)
{
    const Matrix9 tiny_cell = {{1.0e-200, 0.0, 0.0,
                                0.0, 1.0e-200, 0.0,
                                0.0, 0.0, 1.0e-200}};
    const Matrix9 huge_inverse = {{1.0e200, 0.0, 0.0,
                                   0.0, 1.0e200, 0.0,
                                   0.0, 0.0, 1.0e200}};

    const CellValidation out = validate_with_driver_thresholds(tiny_cell, huge_inverse);

    EXPECT_FALSE(out.ok);
    EXPECT_NE(std::string::npos, out.message.find("determinant"));
}

TEST(SocketFrameTest, ConditionNumberMustBeStrictlyBelowMaximum)
{
    const Matrix9 below = {{1.0, 0.0, 0.0, 0.0, 1.0e-6, 0.0, 0.0, 0.0, 2.0e-12}};
    const Matrix9 below_inverse = {{1.0, 0.0, 0.0, 0.0, 1.0e6, 0.0, 0.0, 0.0, 5.0e11}};
    const Matrix9 at = {{1.0, 0.0, 0.0, 0.0, 1.0e-6, 0.0, 0.0, 0.0, 1.0e-12}};
    const Matrix9 at_inverse = {{1.0, 0.0, 0.0, 0.0, 1.0e6, 0.0, 0.0, 0.0, 1.0e12}};
    const Matrix9 above = {{1.0, 0.0, 0.0, 0.0, 1.0e-6, 0.0, 0.0, 0.0, 5.0e-13}};
    const Matrix9 above_inverse = {{1.0, 0.0, 0.0, 0.0, 1.0e6, 0.0, 0.0, 0.0, 2.0e12}};

    EXPECT_TRUE(validate_with_driver_thresholds(below, below_inverse).ok);
    const CellValidation boundary = validate_with_driver_thresholds(at, at_inverse);
    EXPECT_FALSE(boundary.ok);
    EXPECT_NE(std::string::npos, boundary.message.find("condition"));
    EXPECT_DOUBLE_EQ(1.0e12, boundary.condition_number_2);
    EXPECT_FALSE(validate_with_driver_thresholds(above, above_inverse).ok);
}

TEST(SocketFrameTest, NonfiniteCellOrReceivedInverseIsRejected)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();
    const Matrix9 identity = {{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}};
    Matrix9 bad_cell = identity;
    Matrix9 bad_inverse = identity;
    bad_cell[4] = nan;
    EXPECT_FALSE(validate_with_driver_thresholds(bad_cell, identity).ok);
    bad_cell = identity;
    bad_cell[7] = infinity;
    EXPECT_FALSE(validate_with_driver_thresholds(bad_cell, identity).ok);
    bad_inverse[1] = nan;
    EXPECT_FALSE(validate_with_driver_thresholds(identity, bad_inverse).ok);
    bad_inverse = identity;
    bad_inverse[8] = -infinity;
    EXPECT_FALSE(validate_with_driver_thresholds(identity, bad_inverse).ok);
}

TEST(SocketFrameTest, PositionCountRequiresMatchingNonnegativeAtomCount)
{
    std::size_t coordinate_count = 77;
    std::string message;

    EXPECT_FALSE(checked_position_count(-1, 2, coordinate_count, message));
    EXPECT_EQ(77u, coordinate_count);
    EXPECT_NE(std::string::npos, message.find("match"));

    message.clear();
    EXPECT_FALSE(checked_position_count(-1, -1, coordinate_count, message));
    EXPECT_EQ(77u, coordinate_count);
    EXPECT_NE(std::string::npos, message.find("negative"));

    message.clear();
    EXPECT_TRUE(checked_position_count(3, 3, coordinate_count, message)) << message;
    EXPECT_EQ(9u, coordinate_count);
    EXPECT_TRUE(message.empty());
}

TEST(SocketFrameTest, PositionCountRejectsMismatchBeforeDerivingAllocationSize)
{
    std::size_t coordinate_count = 123;
    std::string message;

    EXPECT_FALSE(checked_position_count(std::numeric_limits<std::int32_t>::max(),
                                        1,
                                        coordinate_count,
                                        message));
    EXPECT_EQ(123u, coordinate_count);
    EXPECT_NE(std::string::npos, message.find("match"));
}

TEST(SocketFrameTest, PositionsRequireExactSizeAndFiniteCoordinates)
{
    std::string message;
    const std::vector<double> valid = {1.0, -2.0, 3.0};
    EXPECT_TRUE(validate_positions(valid, 3, message)) << message;

    message.clear();
    EXPECT_FALSE(validate_positions(valid, 6, message));
    EXPECT_NE(std::string::npos, message.find("count"));

    std::vector<double> nonfinite = valid;
    nonfinite[1] = std::numeric_limits<double>::quiet_NaN();
    message.clear();
    EXPECT_FALSE(validate_positions(nonfinite, 3, message));
    EXPECT_NE(std::string::npos, message.find("finite"));

    nonfinite[1] = std::numeric_limits<double>::infinity();
    message.clear();
    EXPECT_FALSE(validate_positions(nonfinite, 3, message));
    EXPECT_NE(std::string::npos, message.find("finite"));
}

TEST(SocketFrameTest, SmallStressAsymmetryIsAveragedBeforeConversion)
{
    const Matrix9 stress = {{1.0, 2.1, 3.2,
                             1.9, 5.0, 6.3,
                             2.8, 5.7, 9.0}};
    const Matrix9 expected = {{1.0, 2.0, 3.0,
                               2.0, 5.0, 6.0,
                               3.0, 6.0, 9.0}};

    const VirialConversion out = make_ipi_virial(stress, 2.0, 0.61, 0.0);

    ASSERT_TRUE(out.ok) << out.message;
    expect_matrix_near(expected, out.wire_virial_hartree, 4.0 * EPSILON);
    EXPECT_NEAR(0.6, out.max_antisymmetric_component, 4.0 * EPSILON);
}

TEST(SocketFrameTest, ExcessiveStressAsymmetryIsRejected)
{
    const Matrix9 stress = {{1.0, 2.1, 3.2,
                             1.9, 5.0, 6.3,
                             2.8, 5.7, 9.0}};

    const VirialConversion out = make_ipi_virial(stress, 2.0, 0.59, 0.0);

    EXPECT_FALSE(out.ok);
    EXPECT_NE(std::string::npos, out.message.find("symmetric"));
    EXPECT_NEAR(0.6, out.max_antisymmetric_component, 4.0 * EPSILON);
}

TEST(SocketFrameTest, StressAsymmetryUsesAbsolutePlusRelativeTolerance)
{
    const Matrix9 accepted_stress = {{10.0, 2.0 + 4.0e-8, 3.0,
                                      2.0 - 4.0e-8, 5.0, 6.0,
                                      3.0, 6.0, 9.0}};
    Matrix9 rejected_stress = accepted_stress;
    rejected_stress[1] = 2.0 + 6.0e-8;
    rejected_stress[3] = 2.0 - 6.0e-8;
    const Matrix9 expected = {{10.0, 2.0, 3.0,
                               2.0, 5.0, 6.0,
                               3.0, 6.0, 9.0}};

    const VirialConversion accepted
        = make_ipi_virial(accepted_stress, 2.0, 1.0e-10, 1.0e-8);
    const VirialConversion rejected
        = make_ipi_virial(rejected_stress, 2.0, 1.0e-10, 1.0e-8);

    ASSERT_TRUE(accepted.ok) << accepted.message;
    expect_matrix_near(expected, accepted.wire_virial_hartree, 4.0 * EPSILON);
    EXPECT_NEAR(8.0e-8, accepted.max_antisymmetric_component, EPSILON);
    EXPECT_FALSE(rejected.ok);
    EXPECT_NE(std::string::npos, rejected.message.find("symmetric"));
    EXPECT_NEAR(1.2e-7, rejected.max_antisymmetric_component, EPSILON);
}

TEST(SocketFrameTest, NonpositiveOrNonfiniteVolumeIsRejected)
{
    const Matrix9 zero_stress = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

    EXPECT_FALSE(make_ipi_virial(zero_stress, 0.0, 1.0e-10, 1.0e-8).ok);
    EXPECT_FALSE(make_ipi_virial(zero_stress, -1.0, 1.0e-10, 1.0e-8).ok);
    EXPECT_FALSE(make_ipi_virial(zero_stress,
                                 std::numeric_limits<double>::infinity(),
                                 1.0e-10,
                                 1.0e-8)
                     .ok);
}

TEST(SocketFrameTest, FiniteStressAndVolumeRejectConvertedVirialOverflow)
{
    const double largest_finite = std::numeric_limits<double>::max();
    const Matrix9 stress = {{largest_finite, 0.0, 0.0,
                             0.0, 1.0, 0.0,
                             0.0, 0.0, 1.0}};

    const VirialConversion out = make_ipi_virial(stress, 4.0, 1.0e-10, 1.0e-8);

    EXPECT_FALSE(out.ok);
    EXPECT_NE(std::string::npos, out.message.find("representable"));
}

TEST(SocketFrameTest, NonfiniteStressIsRejected)
{
    Matrix9 stress = {{1.0, 2.0, 3.0, 2.0, 5.0, 6.0, 3.0, 6.0, 9.0}};
    stress[2] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(make_ipi_virial(stress, 4.0, 1.0e-10, 1.0e-8).ok);
    stress[2] = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(make_ipi_virial(stress, 4.0, 1.0e-10, 1.0e-8).ok);
}
