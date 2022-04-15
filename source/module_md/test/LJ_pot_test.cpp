#include "gtest/gtest.h"
#include "setcell.h"
#include "module_md/LJ_potential.h"

class LJ_pot_test : public testing::Test
{
protected:
    ModuleBase::Vector3<double> *force;
    ModuleBase::matrix stress;
    double potential;
    int natom;

    void SetUp()
    {
        UnitCell_pseudo ucell;
        Setcell::setupcell(ucell);

        natom = ucell.nat;
        force = new ModuleBase::Vector3<double> [natom];
        stress.create(3,3);

        Setcell::parameters();
        Grid_Driver grid_neigh(0,0,0);
        Setcell::neighbor(grid_neigh, ucell);

        potential = LJ_potential::Lennard_Jones(
                                ucell,
                                grid_neigh,
                                force,
                                stress);
    }

    void TearDown()
    {
        delete []force;
    }
};

TEST_F(LJ_pot_test, potential)
{
    EXPECT_DOUBLE_EQ(potential, -0.011957818623534381);
}

TEST_F(LJ_pot_test, force)
{
    EXPECT_DOUBLE_EQ(force[0].x, 0.00049817733089377704);
    EXPECT_DOUBLE_EQ(force[0].y, 0.00082237246837022328);
    EXPECT_DOUBLE_EQ(force[0].z, -3.0493186101154812e-20);
    EXPECT_DOUBLE_EQ(force[1].x, -0.00064758615201580339);
    EXPECT_DOUBLE_EQ(force[1].y, -0.00066924999462089304);
    EXPECT_DOUBLE_EQ(force[1].z, -1.8634724839594607e-20);
    EXPECT_DOUBLE_EQ(force[2].x, -0.00035411224839165616);
    EXPECT_DOUBLE_EQ(force[2].y, 0.0008091080910885112);
    EXPECT_DOUBLE_EQ(force[2].z, -1.1858461261560205e-20);
    EXPECT_DOUBLE_EQ(force[3].x, 0.00050352106951368229);
    EXPECT_DOUBLE_EQ(force[3].y, -0.00096223056483784122);
    EXPECT_DOUBLE_EQ(force[3].z, 2.0328790734103208e-20);
}

TEST_F(LJ_pot_test, stress)
{
    EXPECT_DOUBLE_EQ(stress(0,0), 8.0360222227631859e-07);
    EXPECT_DOUBLE_EQ(stress(0,1), 1.7207745586539077e-07);
    EXPECT_DOUBLE_EQ(stress(0,2), 0);
    EXPECT_DOUBLE_EQ(stress(1,0), 1.7207745586539077e-07);
    EXPECT_DOUBLE_EQ(stress(1,1), 1.0630708613186662e-06);
    EXPECT_DOUBLE_EQ(stress(1,2), -1.1858461261560206e-22);
    EXPECT_DOUBLE_EQ(stress(2,0), 0);
    EXPECT_DOUBLE_EQ(stress(2,1), -1.1858461261560206e-22);
    EXPECT_DOUBLE_EQ(stress(2,2), 6.4275429572682057e-07);
}