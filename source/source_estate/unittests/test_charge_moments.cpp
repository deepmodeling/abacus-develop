#include "source_estate/charge_moments.h"

#include <gtest/gtest.h>

TEST(ChargeMoments, IntegratesAddsAndAcceptsAnEmptyLocalGrid)
{
    const double charges[2] = {2.0, -1.0};
    const ModuleBase::Vector3<double> positions[2] = {
        ModuleBase::Vector3<double>(1.0, 0.0, 0.0),
        ModuleBase::Vector3<double>(0.0, 2.0, 0.0)};
    const elecstate::ChargeMoments moments = elecstate::charge_moments(charges, positions, 2, 0.5);
    EXPECT_DOUBLE_EQ(moments.charge, 0.5);
    EXPECT_DOUBLE_EQ(moments.dipole.x, 1.0);
    EXPECT_DOUBLE_EQ(moments.dipole.y, -1.0);
    EXPECT_DOUBLE_EQ(moments.dipole.z, 0.0);
    EXPECT_DOUBLE_EQ(moments.second_moment, -1.0);

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

    electrons = elecstate::charge_moments(nullptr, nullptr, 0, 1.0);
    EXPECT_DOUBLE_EQ(electrons.charge, 0.0);
    EXPECT_DOUBLE_EQ(electrons.second_moment, 0.0);
}

