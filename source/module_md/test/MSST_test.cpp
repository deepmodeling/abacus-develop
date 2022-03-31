#include "gtest/gtest.h"
#include "setcell.h"
#include "module_md/MSST.h"

class MSST_test : public testing::Test
{
protected:
    Verlet *verlet;
    UnitCell_pseudo ucell;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();
        verlet = new MSST(INPUT.mdp, ucell);

        ModuleESolver::ESolver *p_esolver;
        verlet->setup(p_esolver);
    }

    void TearDown()
    {
        delete verlet;
    }
};

TEST_F(MSST_test, setup)
{
    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.0001314186733659715);
    EXPECT_DOUBLE_EQ(verlet->vel[0].y, 7.0985331994796372e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.3947731701005279e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.00015227275651566311);
    EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.00014579875939315496);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 9.5965690649087203e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00013311885204189453);
    EXPECT_DOUBLE_EQ(verlet->vel[2].y, -3.0298400368294885e-06);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -5.3828659173134662e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 0.00011226476889319793);
    EXPECT_DOUBLE_EQ(verlet->vel[3].y, 7.7843267435287586e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -2.8189299775046767e-05);
    
    EXPECT_DOUBLE_EQ(verlet->stress(0,0), 5.9579909955800075e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(0,1), -1.4582038138067117e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(0,2), 1.4889583894898544e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(1,0), -1.4582038138067117e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(1,1), 3.4199108345556597e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(1,2), -1.2389007575245785e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(2,0), 1.4889583894898544e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(2,1), -1.2389007575245785e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(2,2), 1.5964231736437884e-06);
}

TEST_F(MSST_test, first_half)
{
    verlet->first_half();

    EXPECT_DOUBLE_EQ(ucell.lat0, 1.0);
    EXPECT_DOUBLE_EQ(ucell.lat0_angstrom, 0.52917700000000001);
    EXPECT_DOUBLE_EQ(ucell.latvec.e11, 10.0);
    EXPECT_DOUBLE_EQ(ucell.latvec.e12, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e13, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e21, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e22, 10.0);
    EXPECT_DOUBLE_EQ(ucell.latvec.e23, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e31, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e32, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e33, 9.9959581179144905);
    EXPECT_DOUBLE_EQ(ucell.omega, 999.59581179144902);
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.9945728176928519);
    EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.0029442816868202821);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 9.9953814995781549);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2062875654254492);
    EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1939646253791674);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 0.0039673531204680191);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0944925283175158);
    EXPECT_DOUBLE_EQ(verlet->pos[2].y, 9.9998842371692671);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.9957537084439876);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.0046470885642230716);
    EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.3032068557647492);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9968136746863676);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00013127726219846624);
    EXPECT_DOUBLE_EQ(verlet->vel[0].y, 7.121876825065284e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.3947730561390963e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.0001520889345949577);
    EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.00014598873074918282);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 9.596568280810794e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00013321936931429457);
    EXPECT_DOUBLE_EQ(verlet->vel[2].y, -2.8001689685103039e-06);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -5.3828654775006574e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 0.00011240769691879813);
    EXPECT_DOUBLE_EQ(verlet->vel[3].y, 7.7570131467139791e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -2.8189297471809918e-05);
}

TEST_F(MSST_test, second_half)
{
    verlet->first_half();
    verlet->second_half();

    EXPECT_DOUBLE_EQ(ucell.lat0, 1.0);
    EXPECT_DOUBLE_EQ(ucell.lat0_angstrom, 0.52917700000000001);
    EXPECT_DOUBLE_EQ(ucell.latvec.e11, 10.0);
    EXPECT_DOUBLE_EQ(ucell.latvec.e12, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e13, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e21, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e22, 10.0);
    EXPECT_DOUBLE_EQ(ucell.latvec.e23, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e31, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e32, 0.00);
    EXPECT_DOUBLE_EQ(ucell.latvec.e33, 9.9959581179144905);
    EXPECT_DOUBLE_EQ(ucell.omega, 999.59581179144902);
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.9945728176928519);
    EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.0029442816868202821);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 9.9953814995781549);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2062875654254492);
    EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1939646253791674);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 0.0039673531204680191);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0944925283175158);
    EXPECT_DOUBLE_EQ(verlet->pos[2].y, 9.9998842371692671);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.9957537084439876);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.0046470885642230716);
    EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.3032068557647492);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9968136746863676);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00013113585103096098);
    EXPECT_DOUBLE_EQ(verlet->vel[0].y, 7.1452204506509308e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.3953371489538059e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.00015190511267425228);
    EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.00014617870210521068);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 9.600449453585996e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00013331988658669462);
    EXPECT_DOUBLE_EQ(verlet->vel[2].y, -2.5704979001911192e-06);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -5.3850424881082548e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 0.00011255062494439833);
    EXPECT_DOUBLE_EQ(verlet->vel[3].y, 7.7296995498991997e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -2.8200698165338931e-05);
}