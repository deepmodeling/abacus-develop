#include "source_cell/cell_geometry.h"
#include "source_base/matrix3.h"

#include <gtest/gtest.h>

#include <cmath>

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
        ASSERT_TRUE(unitcell::make_orthogonal_cell(lattice, 10.0, 1.0e-10, cell));
    }

    unitcell::OrthogonalCell cell;
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
    ASSERT_TRUE(unitcell::make_orthogonal_cell(lattice, 10.0, 1.0e-10, cell));
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
    const ModuleBase::Vector3<double> center = unitcell::weighted_center(positions, weights, cell);
    const ModuleBase::Vector3<double> expected(0.1, 0.15, 0.2);
    expect_vector(center, expected);
}

TEST_F(CellGeometryTest, PreservesRelativePositionsUnderWrappedTranslation)
{
    const std::vector<ModuleBase::Vector3<double>> original = {
        ModuleBase::Vector3<double>(9.8, 9.6, 0.1),
        ModuleBase::Vector3<double>(0.3, 0.4, 9.7)};
    const std::vector<double> weights = {2.0, 1.0};
    unitcell::OrthogonalCell first = cell;
    first.origin = unitcell::weighted_center(original, weights, cell);

    std::vector<ModuleBase::Vector3<double>> translated = original;
    for (auto& position : translated)
    {
        const double translated_x = position.x + 1.4;
        const double translated_y = position.y + 1.4;
        const double translated_z = position.z + 1.4;
        position.x = std::fmod(translated_x, 10.0);
        position.y = std::fmod(translated_y, 10.0);
        position.z = std::fmod(translated_z, 10.0);
    }
    unitcell::OrthogonalCell second = cell;
    second.origin = unitcell::weighted_center(translated, weights, cell);
    for (std::size_t index = 0; index < original.size(); ++index)
    {
        const ModuleBase::Vector3<double> expected = unitcell::relative_position(original[index], first);
        const ModuleBase::Vector3<double> actual = unitcell::relative_position(translated[index], second);
        expect_vector(actual, expected);
    }
}

TEST_F(CellGeometryTest, RejectsSkewCell)
{
    const ModuleBase::Matrix3 skew(1.0, 0.0, 0.0,
                                  0.1, 1.0, 0.0,
                                  0.0, 0.0, 1.0);
    EXPECT_FALSE(unitcell::make_orthogonal_cell(skew, 10.0, 1.0e-10, cell));
}

TEST(SlabGeometry, SupportsSkewPeriodicPlaneAndRotatedNormal)
{
    const ModuleBase::Matrix3 lattice(0.0, 2.0, 0.0,
                                     0.0, 1.0, 3.0,
                                     8.0, 0.0, 0.0);
    unitcell::SlabCell cell;
    ASSERT_TRUE(unitcell::make_slab_cell(lattice, 1.0, 2, 1.0e-6, cell));
    EXPECT_DOUBLE_EQ(cell.length, 8.0);
    EXPECT_DOUBLE_EQ(cell.area, 6.0);
    EXPECT_DOUBLE_EQ(cell.normal.x, 1.0);
    const ModuleBase::Vector3<double> point(8.1, 100.0, -200.0);
    EXPECT_NEAR(unitcell::relative_coordinate(point, cell), -3.9, 1.0e-14);
    const std::vector<ModuleBase::Vector3<double>> positions = {
        ModuleBase::Vector3<double>(7.8, 0.0, 0.0),
        ModuleBase::Vector3<double>(0.2, 40.0, 50.0)};
    const std::vector<double> weights = {1.0, 3.0};
    const double center = unitcell::weighted_center(positions, weights, cell);
    EXPECT_NEAR(center, 0.1, 1.0e-14);
}

TEST(SlabGeometry, SelectsAllAxesAndRejectsTilt)
{
    const ModuleBase::Matrix3 lattice(2.0, 0.0, 0.0,
                                     0.0, 3.0, 0.0,
                                     0.0, 0.0, 8.0);
    const double lengths[3] = {2.0, 3.0, 8.0};
    unitcell::SlabCell cell;
    for (int axis = 0; axis < 3; ++axis)
    {
        ASSERT_TRUE(unitcell::make_slab_cell(lattice, 1.0, axis, 1.0e-6, cell));
        EXPECT_DOUBLE_EQ(cell.length, lengths[axis]);
        const double expected_area = 48.0 / lengths[axis];
        EXPECT_DOUBLE_EQ(cell.area, expected_area);
    }
    const ModuleBase::Matrix3 tilted(2.0, 0.0, 0.0,
                                    0.0, 3.0, 0.0,
                                    0.1, 0.0, 8.0);
    EXPECT_FALSE(unitcell::make_slab_cell(tilted, 1.0, 2, 1.0e-6, cell));
}
