#include <gtest/gtest.h>

#include "source_base/mathzone.h"
#include "source_base/matrix.h"

TEST(MathzoneNetForceTest, SumsForceComponents)
{
    ModuleBase::matrix force(3, 3);
    force(0, 0) = 1.0;
    force(0, 1) = -2.0;
    force(0, 2) = 3.0;
    force(1, 0) = -4.0;
    force(1, 1) = 5.0;
    force(1, 2) = -6.0;
    force(2, 0) = 7.0;
    force(2, 1) = -8.0;
    force(2, 2) = 9.0;

    const ModuleBase::Vector3<double> net_force = ModuleBase::net_force_vector(3, force);
    EXPECT_DOUBLE_EQ(net_force.x, 4.0);
    EXPECT_DOUBLE_EQ(net_force.y, -5.0);
    EXPECT_DOUBLE_EQ(net_force.z, 6.0);
}

TEST(MathzoneNetForceTest, RemoveNetForceZeroesTheSum)
{
    ModuleBase::matrix force(3, 3);
    force(0, 0) = 1.0;
    force(0, 1) = -2.0;
    force(0, 2) = 3.0;
    force(1, 0) = -4.0;
    force(1, 1) = 5.0;
    force(1, 2) = -6.0;
    force(2, 0) = 7.0;
    force(2, 1) = -8.0;
    force(2, 2) = 9.0;

    ModuleBase::remove_net_force(3, force);

    const ModuleBase::Vector3<double> net_force = ModuleBase::net_force_vector(3, force);
    EXPECT_NEAR(net_force.x, 0.0, 1.0e-12);
    EXPECT_NEAR(net_force.y, 0.0, 1.0e-12);
    EXPECT_NEAR(net_force.z, 0.0, 1.0e-12);
}
