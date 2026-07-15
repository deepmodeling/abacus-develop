#include "gtest/gtest.h"

#define private public
#include "source_io/module_parameter/parameter.h"
#undef private
#include "setcell.h"
#include "source_md/md_state_view.h"

TEST(MdStateViewTest, UnitCellViewFlattensAndMutatesUnderlyingStorage)
{
    UnitCell ucell;
    Setcell::setupcell(ucell);
    ModuleBase::Vector3<double> force[4];
    for (int i = 0; i < 4; ++i)
    {
        force[i].set(0.1 * (i + 1), 0.2 * (i + 1), 0.3 * (i + 1));
    }

    MdStateView view = MdStateView::from_unitcell(ucell, force);

    ASSERT_EQ(view.size(), 4);
    EXPECT_FALSE(view.uses_distributed_storage());
    EXPECT_EQ(view.global_id(2), 2);
    EXPECT_EQ(view.type(2), 0);
    EXPECT_EQ(view.type_index(2), 2);
    EXPECT_DOUBLE_EQ(view.mass(2), 39.948);
    EXPECT_DOUBLE_EQ(view.cart(1).x, ucell.atoms[0].tau[1].x);
    EXPECT_DOUBLE_EQ(view.frac(1).y, ucell.atoms[0].taud[1].y);
    EXPECT_DOUBLE_EQ(view.vel(3).z, ucell.atoms[0].vel[3].z);
    EXPECT_EQ(view.mbl(0).x, 1);
    EXPECT_DOUBLE_EQ(view.force(3).z, 1.2);

    view.vel(0).x = 9.0;
    view.frac(1).z = 0.25;
    view.mbl(2).y = 0;
    view.force(3).x = 7.0;

    EXPECT_DOUBLE_EQ(ucell.atoms[0].vel[0].x, 9.0);
    EXPECT_DOUBLE_EQ(ucell.atoms[0].taud[1].z, 0.25);
    EXPECT_EQ(ucell.atoms[0].mbl[2].y, 0);
    EXPECT_DOUBLE_EQ(force[3].x, 7.0);
}
