#include "source_estate/pcc_0d.h"
#include "source_base/constants.h"
#include "source_base/matrix3.h"
#include "source_cell/cell_geometry.h"

#include <gtest/gtest.h>

namespace
{
class Pcc0dTest : public testing::Test
{
  protected:
    void SetUp() override
    {
        const ModuleBase::Matrix3 lattice(1.0, 0.0, 0.0,
                                         0.0, 1.0, 0.0,
                                         0.0, 0.0, 1.0);
        ASSERT_TRUE(unitcell::make_orthogonal_cell(lattice, 10.0, 1.0e-10, cell, error));
        ASSERT_TRUE(elecstate::make_pcc_0d_parameters(cell, 1.0e-10, parameters, error));
        moments.charge = 1.2;
        moments.dipole = ModuleBase::Vector3<double>(0.1, -0.2, 0.3);
        moments.second_moment = 0.8;
    }

    unitcell::OrthogonalCell cell;
    elecstate::Pcc0dParameters parameters;
    elecstate::ChargeMoments moments;
    std::string error;
};
} // namespace

TEST_F(Pcc0dTest, SeparatesCubicRestrictionFromOrthogonalGeometry)
{
    EXPECT_DOUBLE_EQ(parameters.length, 10.0);
    const ModuleBase::Matrix3 rectangular(1.0, 0.0, 0.0,
                                         0.0, 2.0, 0.0,
                                         0.0, 0.0, 1.0);
    ASSERT_TRUE(unitcell::make_orthogonal_cell(rectangular, 10.0, 1.0e-10, cell, error));
    EXPECT_FALSE(elecstate::make_pcc_0d_parameters(cell, 1.0e-10, parameters, error));
    EXPECT_DOUBLE_EQ(parameters.length, 10.0);
    EXPECT_FALSE(error.empty());
    EXPECT_FALSE(elecstate::make_pcc_0d_parameters(cell, 0.0, parameters, error));
}

TEST_F(Pcc0dTest, MatchesAnalyticalMonopoleAndNeutralDipoleEnergy)
{
    elecstate::ChargeMoments monopole;
    monopole.charge = 2.0;
    const double actual_monopole = elecstate::pcc_0d_energy(monopole, parameters);
    const double expected_monopole = 0.5 * parameters.madelung * 4.0 / parameters.length;
    EXPECT_NEAR(actual_monopole, expected_monopole, 1.0e-14);

    elecstate::ChargeMoments neutral;
    neutral.dipole = ModuleBase::Vector3<double>(1.0, 2.0, 3.0);
    const double actual_neutral = elecstate::pcc_0d_energy(neutral, parameters);
    const double expected_neutral = 2.0 * ModuleBase::PI * 14.0 / 3000.0;
    EXPECT_NEAR(actual_neutral, expected_neutral, 1.0e-14);
}

TEST_F(Pcc0dTest, BilinearKernelIsSymmetricAndHasConsistentSelfEnergy)
{
    elecstate::ChargeMoments other;
    other.charge = -0.4;
    other.dipole = ModuleBase::Vector3<double>(0.7, 0.5, -0.1);
    other.second_moment = 1.3;
    const double forward = elecstate::pcc_0d_bilinear_energy(moments, other, parameters);
    const double reverse = elecstate::pcc_0d_bilinear_energy(other, moments, parameters);
    EXPECT_NEAR(forward, reverse, 1.0e-14);
    const elecstate::ChargeMoments sum = elecstate::add_charge_moments(moments, other);
    const double sum_energy = elecstate::pcc_0d_energy(sum, parameters);
    const double left_energy = elecstate::pcc_0d_energy(moments, parameters);
    const double right_energy = elecstate::pcc_0d_energy(other, parameters);
    const double expected = left_energy + right_energy + forward;
    EXPECT_NEAR(sum_energy, expected, 1.0e-14);
}

TEST_F(Pcc0dTest, GradientMatchesPotentialFiniteDifferences)
{
    const ModuleBase::Vector3<double> position(0.4, -0.5, 0.6);
    const ModuleBase::Vector3<double> gradient = elecstate::pcc_0d_gradient(moments, position, parameters);
    const double expected[3] = {gradient.x, gradient.y, gradient.z};
    const double step = 1.0e-5;
    for (int axis = 0; axis < 3; ++axis)
    {
        ModuleBase::Vector3<double> plus = position;
        ModuleBase::Vector3<double> minus = position;
        double* plus_component[3] = {&plus.x, &plus.y, &plus.z};
        double* minus_component[3] = {&minus.x, &minus.y, &minus.z};
        *plus_component[axis] += step;
        *minus_component[axis] -= step;
        const double plus_potential = elecstate::pcc_0d_potential(moments, plus, parameters);
        const double minus_potential = elecstate::pcc_0d_potential(moments, minus, parameters);
        const double derivative = (plus_potential - minus_potential) / (2.0 * step);
        EXPECT_NEAR(derivative, expected[axis], 1.0e-10);
    }
}

TEST_F(Pcc0dTest, IonicForceMatchesTotalEnergyFiniteDifferences)
{
    const double charges[2] = {1.3, -0.4};
    ModuleBase::Vector3<double> positions[2] = {
        ModuleBase::Vector3<double>(0.2, -0.3, 0.4),
        ModuleBase::Vector3<double>(-0.5, 0.6, -0.7)};
    ASSERT_TRUE(elecstate::charge_moments(charges, positions, 2, 1.0, moments, error));
    const ModuleBase::Vector3<double> force = elecstate::pcc_0d_force(moments, charges[0], positions[0], parameters);
    const double expected[3] = {force.x, force.y, force.z};
    const double step = 1.0e-5;
    for (int axis = 0; axis < 3; ++axis)
    {
        double* component[3] = {&positions[0].x, &positions[0].y, &positions[0].z};
        const double original = *component[axis];
        *component[axis] = original + step;
        ASSERT_TRUE(elecstate::charge_moments(charges, positions, 2, 1.0, moments, error));
        const double plus_energy = elecstate::pcc_0d_energy(moments, parameters);
        *component[axis] = original - step;
        ASSERT_TRUE(elecstate::charge_moments(charges, positions, 2, 1.0, moments, error));
        const double minus_energy = elecstate::pcc_0d_energy(moments, parameters);
        *component[axis] = original;
        const double numerical_force = -(plus_energy - minus_energy) / (2.0 * step);
        EXPECT_NEAR(numerical_force, expected[axis], 1.0e-10);
    }
}

TEST_F(Pcc0dTest, ElectronPotentialMatchesEnergyDerivative)
{
    const ModuleBase::Vector3<double> position(0.2, -0.3, 0.4);
    const double positive_potential = elecstate::pcc_0d_potential(moments, position, parameters);
    const double electron_potential_ry = -2.0 * positive_potential;
    const double step = 1.0e-5;
    elecstate::ChargeMoments plus = moments;
    plus.charge -= step;
    plus.dipole -= position * step;
    plus.second_moment -= step * position.norm2();
    elecstate::ChargeMoments minus = moments;
    minus.charge += step;
    minus.dipole += position * step;
    minus.second_moment += step * position.norm2();
    const double plus_energy = 2.0 * elecstate::pcc_0d_energy(plus, parameters);
    const double minus_energy = 2.0 * elecstate::pcc_0d_energy(minus, parameters);
    const double derivative = (plus_energy - minus_energy) / (2.0 * step);
    EXPECT_NEAR(derivative, electron_potential_ry, 1.0e-10);
}

TEST_F(Pcc0dTest, EnergyIsIndependentOfMultipoleOrigin)
{
    const double charges[2] = {1.2, -0.3};
    ModuleBase::Vector3<double> positions[2] = {
        ModuleBase::Vector3<double>(0.3, -0.2, 0.1),
        ModuleBase::Vector3<double>(-0.5, 0.6, -0.7)};
    ASSERT_TRUE(elecstate::charge_moments(charges, positions, 2, 1.0, moments, error));
    const double original_energy = elecstate::pcc_0d_energy(moments, parameters);
    const ModuleBase::Vector3<double> shift(0.4, -0.3, 0.2);
    positions[0] -= shift;
    positions[1] -= shift;
    ASSERT_TRUE(elecstate::charge_moments(charges, positions, 2, 1.0, moments, error));
    const double shifted_energy = elecstate::pcc_0d_energy(moments, parameters);
    EXPECT_NEAR(original_energy, shifted_energy, 1.0e-14);
}
