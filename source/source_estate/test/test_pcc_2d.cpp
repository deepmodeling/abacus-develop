#include "source_estate/pcc_2d.h"
#include "source_base/constants.h"
#include "source_cell/cell_geometry.h"

#include <gtest/gtest.h>

TEST(Pcc2d, MatchesChargedSheetGaugeAndNeutralDipoleEnergy)
{
    unitcell::SlabCell cell;
    cell.area = 6.0;
    cell.length = 8.0;
    elecstate::Pcc2dParameters parameters;
    std::string error;
    ASSERT_TRUE(elecstate::make_pcc_2d_parameters(cell, parameters, error));
    elecstate::ChargeMoments moments;
    moments.charge = 2.0;
    const double charged = elecstate::pcc_2d_energy(moments, parameters);
    const double expected_charged = -ModuleBase::PI * 8.0 * 4.0 / 36.0;
    EXPECT_NEAR(charged, expected_charged, 1.0e-14);
    moments.charge = 0.0;
    moments.dipole.x = 2.0;
    const double neutral = elecstate::pcc_2d_energy(moments, parameters);
    const double expected_neutral = 2.0 * ModuleBase::PI * 4.0 / 48.0;
    EXPECT_NEAR(neutral, expected_neutral, 1.0e-14);
    cell.area = 0.0;
    EXPECT_FALSE(elecstate::make_pcc_2d_parameters(cell, parameters, error));
    EXPECT_DOUBLE_EQ(parameters.area, 6.0);
}

TEST(Pcc2d, NormalForceMatchesEnergyDerivativeAndHasNoTangentialComponent)
{
    elecstate::Pcc2dParameters parameters;
    parameters.area = 6.0;
    parameters.length = 8.0;
    const ModuleBase::Vector3<double> normal(0.0, 1.0, 0.0);
    const double charges[2] = {1.3, -0.4};
    ModuleBase::Vector3<double> positions[2] = {
        ModuleBase::Vector3<double>(0.0, 0.2, 0.0),
        ModuleBase::Vector3<double>(0.0, -0.5, 0.0)};
    elecstate::ChargeMoments moments;
    std::string error;
    ASSERT_TRUE(elecstate::charge_moments(charges, positions, 2, 1.0, moments, error));
    const double coordinate = positions[0].y;
    const ModuleBase::Vector3<double> force = elecstate::pcc_2d_force(moments, charges[0], coordinate, normal, parameters);
    EXPECT_DOUBLE_EQ(force.x, 0.0);
    EXPECT_DOUBLE_EQ(force.z, 0.0);
    const double step = 1.0e-5;
    positions[0].y = coordinate + step;
    ASSERT_TRUE(elecstate::charge_moments(charges, positions, 2, 1.0, moments, error));
    const double plus = elecstate::pcc_2d_energy(moments, parameters);
    positions[0].y = coordinate - step;
    ASSERT_TRUE(elecstate::charge_moments(charges, positions, 2, 1.0, moments, error));
    const double minus = elecstate::pcc_2d_energy(moments, parameters);
    const double numerical_force = -(plus - minus) / (2.0 * step);
    EXPECT_NEAR(force.y, numerical_force, 1.0e-10);
}

TEST(Pcc2d, ElectronicPotentialMatchesEnergyDerivative)
{
    elecstate::Pcc2dParameters parameters;
    parameters.area = 6.0;
    parameters.length = 8.0;
    const ModuleBase::Vector3<double> normal(1.0, 0.0, 0.0);
    const double coordinate = 0.3;
    elecstate::ChargeMoments moments;
    moments.charge = 1.2;
    moments.dipole.x = 0.2;
    moments.second_moment = 0.8;
    const double positive_potential = elecstate::pcc_2d_potential(moments, coordinate, normal, parameters);
    const double expected = -2.0 * positive_potential;
    const double step = 1.0e-5;
    elecstate::ChargeMoments plus = moments;
    plus.charge -= step;
    plus.dipole.x -= coordinate * step;
    plus.second_moment -= coordinate * coordinate * step;
    elecstate::ChargeMoments minus = moments;
    minus.charge += step;
    minus.dipole.x += coordinate * step;
    minus.second_moment += coordinate * coordinate * step;
    const double plus_energy = 2.0 * elecstate::pcc_2d_energy(plus, parameters);
    const double minus_energy = 2.0 * elecstate::pcc_2d_energy(minus, parameters);
    const double derivative = (plus_energy - minus_energy) / (2.0 * step);
    EXPECT_NEAR(derivative, expected, 1.0e-10);
}
