#include "gtest/gtest.h"
#include "setcell.h"
#include "module_md/FIRE.h"

class FIRE_test : public testing::Test
{
protected:
    Verlet *verlet;
    UnitCell_pseudo ucell;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();
        verlet = new FIRE(INPUT.mdp, ucell);

        ModuleESolver::ESolver *p_esolver;
        verlet->setup(p_esolver);
    }

    void TearDown()
    {
        delete verlet;
    }
};

TEST_F(FIRE_test, setup)
{
    EXPECT_DOUBLE_EQ(verlet->temperature_, 299.99999999999665);
    EXPECT_DOUBLE_EQ(verlet->stress(0,0), 6.0100555286436806e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(0,1), -1.4746713013791574e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(0,2), 1.5039983732220751e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(1,0), -1.4746713013791574e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(1,1), 3.4437172989317909e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(1,2), -1.251414906590483e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(2,0), 1.5039983732220751e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(2,1), -1.251414906590483e-06);
    EXPECT_DOUBLE_EQ(verlet->stress(2,2), 1.6060561926126463e-06);
}

TEST_F(FIRE_test, first_half)
{
    verlet->first_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.995455294044568);
    EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.003264683323249327);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 9.9994784290476932);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2052136746814073);
    EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.19405131115556);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 0.0035886062145122004);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0947593079696478);
    EXPECT_DOUBLE_EQ(verlet->pos[2].y, 0.00048706739346586134);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.9979870945593055);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.0045717233044145923);
    EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.3021969381277305);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9989458701784848);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00010993118004167345);
    EXPECT_DOUBLE_EQ(verlet->vel[0].y, 7.8968913216100539e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.2616198016939999e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.00012611275970351733);
    EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.00014389190209072655);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 8.6804233262820007e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00012676627812260489);
    EXPECT_DOUBLE_EQ(verlet->vel[2].y, 1.1781596840062159e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -4.8689854212330001e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 0.00011058469846166102);
    EXPECT_DOUBLE_EQ(verlet->vel[3].y, 5.3141392034653857e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -2.5498181033639999e-05);
}

TEST_F(FIRE_test, second_half)
{
    verlet->first_half();
    verlet->second_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.995455294044568);
    EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.003264683323249327);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 9.9994784290476932);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2052136746814073);
    EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.19405131115556);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 0.0035886062145122004);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0947593079696478);
    EXPECT_DOUBLE_EQ(verlet->pos[2].y, 0.00048706739346586134);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.9979870945593055);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.0045717233044145923);
    EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.3021969381277305);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9989458701784848);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00010978976887416819);
    EXPECT_DOUBLE_EQ(verlet->vel[0].y, 7.9202349471957007e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.2616198016939999e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.00012592893778281191);
    EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.00014408187344675441);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 8.6804233262820007e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00012686679539500493);
    EXPECT_DOUBLE_EQ(verlet->vel[2].y, 1.2011267908381344e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -4.8689854212330001e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 0.00011072762648726122);
    EXPECT_DOUBLE_EQ(verlet->vel[3].y, 5.2868256066506055e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -2.5498181033639999e-05);
}