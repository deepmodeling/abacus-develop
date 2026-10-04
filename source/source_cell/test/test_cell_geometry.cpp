#include "source_cell/cell_geometry.h"
#include "source_base/matrix3.h"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

namespace
{
class CellGeometryTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0,
                                         0.0, 0.0, 1.0);
        ASSERT_TRUE(unitcell::make_orthogonal_cell(lattice, 10.0, 1.0e-10, cell, error));
    }

    unitcell::OrthogonalCell cell;
    std::string error;
};

void expect_vector(const ModuleBase::Vector3<double>& actual,
                   const ModuleBase::Vector3<double>& expected)
{
    EXPECT_NEAR(actual.x, expected.x, 1.0e-13);
    EXPECT_NEAR(actual.y, expected.y, 1.0e-13);
    EXPECT_NEAR(actual.z, expected.z, 1.0e-13);
}
} // namespace

TEST_F(CellGeometryTest, BuildsRotatedRectangularCell)
{
    const ModuleBase::Matrix3 lattice(0.0, 1.0, 0.0,
                                     -2.0, 0.0, 0.0,
                                     0.0, 0.0, 3.0);
    ASSERT_TRUE(unitcell::make_orthogonal_cell(lattice, 10.0, 1.0e-10, cell, error));
    EXPECT_DOUBLE_EQ(cell.lengths[0], 10.0);
    EXPECT_DOUBLE_EQ(cell.lengths[1], 20.0);
    EXPECT_DOUBLE_EQ(cell.lengths[2], 30.0);
    const ModuleBase::Vector3<double> expected_origin(-10.0, 5.0, 15.0);
    expect_vector(cell.origin, expected_origin);

    const ModuleBase::Vector3<double> displacement(11.0, 6.0, 16.0);
    const ModuleBase::Vector3<double> position = cell.origin + displacement;
    const ModuleBase::Vector3<double> actual = unitcell::relative_position(position, cell);
    const ModuleBase::Vector3<double> expected(-9.0, -4.0, -14.0);
    expect_vector(actual, expected);
}

TEST_F(CellGeometryTest, WrapsMultipleImagesAndUsesHalfOpenInterval)
{
    cell.origin = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    const ModuleBase::Vector3<double> position(25.0, -25.0, 31.0);
    const ModuleBase::Vector3<double> actual = unitcell::relative_position(position, cell);
    const ModuleBase::Vector3<double> expected(-5.0, -5.0, 1.0);
    expect_vector(actual, expected);
}

TEST_F(CellGeometryTest, UnwrapsWeightedCenterAcrossCellBoundary)
{
    const std::vector<ModuleBase::Vector3<double>> positions = {
        ModuleBase::Vector3<double>(9.8, 9.7, 9.6),
        ModuleBase::Vector3<double>(0.2, 0.3, 0.4)};
    const std::vector<double> weights = {1.0, 3.0};
    ModuleBase::Vector3<double> center;
    ASSERT_TRUE(unitcell::weighted_center(positions, weights, cell, center, error));
    const ModuleBase::Vector3<double> expected(0.1, 0.15, 0.2);
    expect_vector(center, expected);
}

TEST_F(CellGeometryTest, PreservesRelativePositionsUnderWrappedTranslation)
{
    const std::vector<ModuleBase::Vector3<double>> original = {
        ModuleBase::Vector3<double>(9.8, 9.6, 0.1),
        ModuleBase::Vector3<double>(0.3, 0.4, 9.7)};
    const std::vector<double> weights = {2.0, 1.0};
    ModuleBase::Vector3<double> center;
    ASSERT_TRUE(unitcell::weighted_center(original, weights, cell, center, error));
    unitcell::OrthogonalCell first = cell;
    first.origin = center;

    std::vector<ModuleBase::Vector3<double>> translated = original;
    for (auto& position : translated)
    {
        position.x = std::fmod(position.x + 1.4, 10.0);
        position.y = std::fmod(position.y + 1.4, 10.0);
        position.z = std::fmod(position.z + 1.4, 10.0);
    }
    ASSERT_TRUE(unitcell::weighted_center(translated, weights, cell, center, error));
    unitcell::OrthogonalCell second = cell;
    second.origin = center;
    for (std::size_t index = 0; index < original.size(); ++index)
    {
        const ModuleBase::Vector3<double> expected = unitcell::relative_position(original[index], first);
        const ModuleBase::Vector3<double> actual = unitcell::relative_position(translated[index], second);
        expect_vector(actual, expected);
    }
}

TEST_F(CellGeometryTest, RejectsSkewAndDegenerateCellsWithoutChangingResult)
{
    const ModuleBase::Vector3<double> original_origin = cell.origin;
    const ModuleBase::Matrix3 skew(1.0, 0.0, 0.0,
                                  0.1, 1.0, 0.0,
                                  0.0, 0.0, 1.0);
    EXPECT_FALSE(unitcell::make_orthogonal_cell(skew, 10.0, 1.0e-10, cell, error));
    EXPECT_FALSE(error.empty());
    expect_vector(cell.origin, original_origin);
    const ModuleBase::Matrix3 degenerate(1.0, 0.0, 0.0,
                                        0.0, 0.0, 0.0,
                                        0.0, 0.0, 1.0);
    EXPECT_FALSE(unitcell::make_orthogonal_cell(degenerate, 10.0, 1.0e-10, cell, error));
    expect_vector(cell.origin, original_origin);
}

TEST_F(CellGeometryTest, RejectsInvalidScaleToleranceAndNonfiniteLattice)
{
    const ModuleBase::Matrix3 identity(1.0, 0.0, 0.0,
                                      0.0, 1.0, 0.0,
                                      0.0, 0.0, 1.0);
    const double infinity = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(unitcell::make_orthogonal_cell(identity, 0.0, 1.0e-10, cell, error));
    EXPECT_FALSE(unitcell::make_orthogonal_cell(identity, infinity, 1.0e-10, cell, error));
    EXPECT_FALSE(unitcell::make_orthogonal_cell(identity, 10.0, 0.0, cell, error));
    EXPECT_FALSE(unitcell::make_orthogonal_cell(identity, 10.0, infinity, cell, error));
    const ModuleBase::Matrix3 nonfinite(infinity, 0.0, 0.0,
                                       0.0, 1.0, 0.0,
                                       0.0, 0.0, 1.0);
    EXPECT_FALSE(unitcell::make_orthogonal_cell(nonfinite, 10.0, 1.0e-10, cell, error));
    ASSERT_TRUE(unitcell::make_orthogonal_cell(identity, 10.0, 1.0e-10, cell, error));
    EXPECT_TRUE(error.empty());
}

TEST_F(CellGeometryTest, RejectsInvalidCenterInputsWithoutChangingResult)
{
    std::vector<ModuleBase::Vector3<double>> positions = {
        ModuleBase::Vector3<double>(1.0, 2.0, 3.0)};
    const std::vector<double> empty_weights;
    const std::vector<double> zero_weights = {0.0};
    const std::vector<double> weights = {1.0};
    const ModuleBase::Vector3<double> original(4.0, 5.0, 6.0);
    ModuleBase::Vector3<double> center = original;
    EXPECT_FALSE(unitcell::weighted_center(positions, empty_weights, cell, center, error));
    EXPECT_FALSE(unitcell::weighted_center(positions, zero_weights, cell, center, error));
    positions[0].x = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(unitcell::weighted_center(positions, weights, cell, center, error));
    expect_vector(center, original);
    EXPECT_FALSE(error.empty());
}

TEST_F(CellGeometryTest, RejectsOverflowingCenterAccumulation)
{
    const std::vector<ModuleBase::Vector3<double>> positions = {
        ModuleBase::Vector3<double>(1.0, 1.0, 1.0),
        ModuleBase::Vector3<double>(2.0, 2.0, 2.0)};
    const double maximum = std::numeric_limits<double>::max();
    const std::vector<double> weights = {maximum, maximum};
    ModuleBase::Vector3<double> center;
    EXPECT_FALSE(unitcell::weighted_center(positions, weights, cell, center, error));
    EXPECT_FALSE(error.empty());
}

TEST(SlabGeometry, SupportsSkewPeriodicPlaneAndRotatedNormal)
{
    const ModuleBase::Matrix3 lattice(0.0, 2.0, 0.0,
                                     0.0, 1.0, 3.0,
                                     8.0, 0.0, 0.0);
    unitcell::SlabCell cell;
    std::string error;
    ASSERT_TRUE(unitcell::make_slab_cell(lattice, 1.0, 2, 1.0e-6, cell, error));
    EXPECT_DOUBLE_EQ(cell.length, 8.0);
    EXPECT_DOUBLE_EQ(cell.area, 6.0);
    EXPECT_DOUBLE_EQ(cell.normal.x, 1.0);
    const ModuleBase::Vector3<double> point(8.1, 100.0, -200.0);
    EXPECT_NEAR(unitcell::relative_coordinate(point, cell), -3.9, 1.0e-14);
    const std::vector<ModuleBase::Vector3<double>> positions = {
        ModuleBase::Vector3<double>(7.8, 0.0, 0.0),
        ModuleBase::Vector3<double>(0.2, 40.0, 50.0)};
    const std::vector<double> weights = {1.0, 3.0};
    double center = 5.0;
    ASSERT_TRUE(unitcell::weighted_center(positions, weights, cell, center, error));
    EXPECT_NEAR(center, 0.1, 1.0e-14);
}

TEST(SlabGeometry, SelectsAllAxesAndRejectsTiltWithoutChangingResult)
{
    const ModuleBase::Matrix3 lattice(2.0, 0.0, 0.0,
                                     0.0, 3.0, 0.0,
                                     0.0, 0.0, 8.0);
    const double lengths[3] = {2.0, 3.0, 8.0};
    unitcell::SlabCell cell;
    std::string error;
    for (int axis = 0; axis < 3; ++axis)
    {
        ASSERT_TRUE(unitcell::make_slab_cell(lattice, 1.0, axis, 1.0e-6, cell, error));
        EXPECT_DOUBLE_EQ(cell.length, lengths[axis]);
        const double expected_area = 48.0 / lengths[axis];
        EXPECT_DOUBLE_EQ(cell.area, expected_area);
    }
    const ModuleBase::Matrix3 tilted(2.0, 0.0, 0.0,
                                    0.0, 3.0, 0.0,
                                    0.1, 0.0, 8.0);
    EXPECT_FALSE(unitcell::make_slab_cell(tilted, 1.0, 2, 1.0e-6, cell, error));
    EXPECT_DOUBLE_EQ(cell.length, 8.0);
    EXPECT_FALSE(unitcell::make_slab_cell(lattice, 1.0, 3, 1.0e-6, cell, error));
    EXPECT_FALSE(unitcell::make_slab_cell(lattice, 0.0, 2, 1.0e-6, cell, error));
    EXPECT_FALSE(error.empty());
}
