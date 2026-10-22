#include "../librpa_stru_symmetry.h"
#include "../librpa_stru_units.h"
#include "../librpa_bz_sampling.h"

#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_base/vector3.h"
#include "source_cell/module_symmetry/symm_rot_spin.h"

#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
TEST(RpaLriSymmetry, WritesUnitaryRowsBeforeAntiunitaryRows)
{
    const std::vector<RpaLriDetail::LibRpaSymmetryOperation> unitary = {
        RpaLriDetail::make_librpa_symmetry_operation(ModuleBase::Matrix3(), ModuleBase::Vector3<double>())};
    const std::vector<RpaLriDetail::LibRpaSymmetryOperation> antiunitary = {
        RpaLriDetail::make_librpa_symmetry_operation(
            ModuleBase::Matrix3(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0),
            ModuleBase::Vector3<double>(0.5, 0.0, 0.0))};

    std::ostringstream output;
    RpaLriDetail::write_librpa_symmetry_rows(output, unitary, antiunitary);
    const std::string text = output.str();
    std::istringstream lines(text);
    std::string line;
    std::vector<std::string> rows;
    while (std::getline(lines, line))
    {
        rows.push_back(line);
    }

    ASSERT_EQ(rows.size(), 3u);
    EXPECT_EQ(rows[0], "2 row");
    std::istringstream unitary_row(rows[1]);
    std::istringstream antiunitary_row(rows[2]);
    int unitary_rotation[9] = {};
    int antiunitary_rotation[9] = {};
    for (int i = 0; i != 9; ++i)
    {
        unitary_row >> unitary_rotation[i];
        antiunitary_row >> antiunitary_rotation[i];
    }
    EXPECT_EQ(antiunitary_rotation[0], -1);
    for (int i = 0; i != 9; ++i)
    {
        const int expected = i == 0 || i == 4 || i == 8 ? 1 : 0;
        EXPECT_EQ(unitary_rotation[i], expected);
    }
    EXPECT_EQ(antiunitary_rotation[4], -1);
    EXPECT_EQ(antiunitary_rotation[8], 1);
    EXPECT_NE(rows[2].find("5.000000000000000e-01"), std::string::npos);
}

TEST(RpaLriSymmetry, WritesSpinRowsWithAntiunitaryFlag)
{
    const std::vector<ModuleBase::Matrix3> unitary = {ModuleBase::Matrix3()};
    const std::vector<ModuleBase::Matrix3> antiunitary = {
        ModuleBase::Matrix3(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0)};
    const auto operations
        = RpaLriDetail::make_librpa_spin_symmetry_operations(ModuleBase::Matrix3(), unitary, antiunitary);

    ASSERT_EQ(operations.size(), 2u);
    EXPECT_EQ(operations[0].antiunitary, 0);
    EXPECT_EQ(operations[1].antiunitary, 1);
    EXPECT_NEAR(operations[0].spin_u[0], 1.0, 1e-12);
    EXPECT_NEAR(operations[0].spin_u[7], 0.0, 1e-12);

    std::ostringstream output;
    RpaLriDetail::write_librpa_spin_symmetry(output, 0, 1, operations);
    const std::string text = output.str();
    EXPECT_NE(text.find("spin_symmetry 0 1\n"), std::string::npos);
    EXPECT_NE(text.find("\n1 "), std::string::npos);
}

TEST(RpaLriSymmetry, UsesPhysicalLatticeConstantForBothUnits)
{
    const auto scales = RpaLriDetail::librpa_stru_unit_scales(5.3);
    EXPECT_DOUBLE_EQ(scales.real_space_bohr, 5.3);
    EXPECT_DOUBLE_EQ(scales.reciprocal_space_bohr_inv, ModuleBase::TWO_PI / 5.3);
    EXPECT_NEAR(scales.real_space_bohr * scales.reciprocal_space_bohr_inv,
                ModuleBase::TWO_PI,
                4.0 * std::numeric_limits<double>::epsilon());
    EXPECT_THROW(RpaLriDetail::librpa_stru_unit_scales(0.0), std::invalid_argument);
}

TEST(RpaLriSymmetry, ConvertsCartesianKPointWithPhysicalScale)
{
    const auto cartesian = RpaLriDetail::librpa_cartesian_kvec(ModuleBase::Vector3<double>(1.0, -2.0, 0.5), 2.0);
    EXPECT_DOUBLE_EQ(cartesian.x, ModuleBase::PI);
    EXPECT_DOUBLE_EQ(cartesian.y, -2.0 * ModuleBase::PI);
    EXPECT_DOUBLE_EQ(cartesian.z, 0.5 * ModuleBase::PI);
    EXPECT_THROW(RpaLriDetail::librpa_cartesian_kvec(ModuleBase::Vector3<double>(), -1.0), std::invalid_argument);
}

TEST(RpaLriSymmetry, ConvertsNonorthogonalLatticeRotationBeforeSpinMapping)
{
    const double half = 0.5;
    const double height = std::sqrt(3.0) * half;
    const ModuleBase::Matrix3 latvec(1.0, 0.0, 0.0, half, height, 0.0, 0.0, 0.0, 1.0);
    // This fractional operation is a 60-degree rotation of the triangular lattice.
    const ModuleBase::Matrix3 fractional_rotation(0.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0, 1.0);
    const auto operations = RpaLriDetail::make_librpa_spin_symmetry_operations(
        latvec, std::vector<ModuleBase::Matrix3>{fractional_rotation}, {});

    ASSERT_EQ(operations.size(), 1u);
    const double half_angle_cos = std::sqrt(3.0) * half;
    EXPECT_NEAR(operations[0].spin_u[0], half_angle_cos, 1e-12);
    EXPECT_NEAR(operations[0].spin_u[1], -half, 1e-12);
    EXPECT_NEAR(operations[0].spin_u[2], 0.0, 1e-12);
    EXPECT_NEAR(operations[0].spin_u[3], 0.0, 1e-12);
    EXPECT_NEAR(operations[0].spin_u[4], 0.0, 1e-12);
    EXPECT_NEAR(operations[0].spin_u[5], 0.0, 1e-12);
    EXPECT_NEAR(operations[0].spin_u[6], half_angle_cos, 1e-12);
    EXPECT_NEAR(operations[0].spin_u[7], half, 1e-12);

    const ModuleBase::Matrix3 cartesian = latvec.Inverse() * fractional_rotation * latvec;
    const auto expected = ModuleSymmetry::SpinRotation::so3_to_su2(cartesian);
    for (int i = 0; i != 4; ++i)
    {
        EXPECT_NEAR(operations[0].spin_u[static_cast<std::size_t>(2 * i)], expected[i].real(), 1e-12);
        EXPECT_NEAR(operations[0].spin_u[static_cast<std::size_t>(2 * i + 1)], expected[i].imag(), 1e-12);
    }
}
} // namespace
