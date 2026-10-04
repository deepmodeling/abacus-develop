#include "source_estate/charge_moments.h"

#include <gtest/gtest.h>

#include <limits>

TEST(ChargeMoments, IntegratesSignedDensityWithVolumeElement)
{
    const double charges[2] = {2.0, -1.0};
    const ModuleBase::Vector3<double> positions[2] = {
        ModuleBase::Vector3<double>(1.0, 0.0, 0.0),
        ModuleBase::Vector3<double>(0.0, 2.0, 0.0)};
    elecstate::ChargeMoments moments;
    std::string error;
    ASSERT_TRUE(elecstate::charge_moments(charges, positions, 2, 0.5, moments, error));
    EXPECT_DOUBLE_EQ(moments.charge, 0.5);
    EXPECT_DOUBLE_EQ(moments.dipole.x, 1.0);
    EXPECT_DOUBLE_EQ(moments.dipole.y, -1.0);
    EXPECT_DOUBLE_EQ(moments.dipole.z, 0.0);
    EXPECT_DOUBLE_EQ(moments.second_moment, -1.0);
}

TEST(ChargeMoments, AddsIndependentContributionsAndAcceptsEmptyLocalGrid)
{
    elecstate::ChargeMoments ions;
    ions.charge = 2.0;
    ions.dipole.x = 3.0;
    ions.second_moment = 4.0;
    elecstate::ChargeMoments electrons;
    electrons.charge = -1.0;
    electrons.dipole.y = -2.0;
    electrons.second_moment = -3.0;
    const elecstate::ChargeMoments total = elecstate::add_charge_moments(ions, electrons);
    EXPECT_DOUBLE_EQ(total.charge, 1.0);
    EXPECT_DOUBLE_EQ(total.dipole.x, 3.0);
    EXPECT_DOUBLE_EQ(total.dipole.y, -2.0);
    EXPECT_DOUBLE_EQ(total.second_moment, 1.0);

    std::string error;
    ASSERT_TRUE(elecstate::charge_moments(nullptr, nullptr, 0, 1.0, electrons, error));
    EXPECT_DOUBLE_EQ(electrons.charge, 0.0);
    EXPECT_DOUBLE_EQ(electrons.second_moment, 0.0);
}

TEST(ChargeMoments, RejectsInvalidStorageAndNonfiniteValuesWithoutChangingResult)
{
    double charge = 1.0;
    ModuleBase::Vector3<double> position(1.0, 2.0, 3.0);
    elecstate::ChargeMoments moments;
    moments.charge = 7.0;
    std::string error;
    EXPECT_FALSE(elecstate::charge_moments(nullptr, &position, 1, 1.0, moments, error));
    EXPECT_FALSE(elecstate::charge_moments(&charge, nullptr, 1, 1.0, moments, error));
    EXPECT_FALSE(elecstate::charge_moments(&charge, &position, -1, 1.0, moments, error));
    EXPECT_FALSE(elecstate::charge_moments(&charge, &position, 1, 0.0, moments, error));
    charge = std::numeric_limits<double>::infinity();
    EXPECT_FALSE(elecstate::charge_moments(&charge, &position, 1, 1.0, moments, error));
    charge = 1.0;
    position.x = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(elecstate::charge_moments(&charge, &position, 1, 1.0, moments, error));
    EXPECT_DOUBLE_EQ(moments.charge, 7.0);
    EXPECT_FALSE(error.empty());
}

TEST(ChargeMoments, RejectsOverflowingMoments)
{
    const double charge = 1.0;
    const double maximum = std::numeric_limits<double>::max();
    const ModuleBase::Vector3<double> position(maximum, 0.0, 0.0);
    elecstate::ChargeMoments moments;
    std::string error;
    EXPECT_FALSE(elecstate::charge_moments(&charge, &position, 1, 1.0, moments, error));
    EXPECT_DOUBLE_EQ(moments.charge, 0.0);
}
