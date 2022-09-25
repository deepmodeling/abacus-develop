#include "gtest/gtest.h"
#include "setcell.h"
#include "module_md/Langevin.h"
#include "module_esolver/esolver_lj.h"

#define doublethreshold 1e-12

class Langevin_test : public testing::Test
{
protected:
    Verlet *verlet;
    UnitCell_pseudo ucell;
    
    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();

        ModuleESolver::ESolver *p_esolver = new ModuleESolver::ESolver_LJ();
        p_esolver->Init(INPUT, ucell);

        verlet = new Langevin(INPUT.mdp, ucell);
        verlet->setup(p_esolver);
    }

    void TearDown()
    {
        delete verlet;
    }
};

TEST_F(Langevin_test, setup)
{
    EXPECT_NEAR(verlet->temperature_, 299.99999999999665, doublethreshold);
    EXPECT_NEAR(verlet->stress(0,0), 6.0100555286436806e-06, doublethreshold);
    EXPECT_NEAR(verlet->stress(0,1), -1.4746713013791574e-06, doublethreshold);
    EXPECT_NEAR(verlet->stress(0,2), 1.5039983732220751e-06, doublethreshold);
    EXPECT_NEAR(verlet->stress(1,0), -1.4746713013791574e-06, doublethreshold);
    EXPECT_NEAR(verlet->stress(1,1), 3.4437172989317909e-06, doublethreshold);
    EXPECT_NEAR(verlet->stress(1,2), -1.251414906590483e-06, doublethreshold);
    EXPECT_NEAR(verlet->stress(2,0), 1.5039983732220751e-06, doublethreshold);
    EXPECT_NEAR(verlet->stress(2,1), -1.251414906590483e-06, doublethreshold);
    EXPECT_NEAR(verlet->stress(2,2), 1.6060561926126463e-06, doublethreshold);
    
    EXPECT_NEAR(verlet->force[0].x, 0.56845590974612858, doublethreshold);
    EXPECT_NEAR(verlet->force[0].y, -0.22894552385324796, doublethreshold);
    EXPECT_NEAR(verlet->force[0].z, 0.30372719169186396, doublethreshold);
    EXPECT_NEAR(verlet->force[1].x, 0.023936136015932585, doublethreshold);
    EXPECT_NEAR(verlet->force[1].y, 0.66317969013068256, doublethreshold);
    EXPECT_NEAR(verlet->force[1].z, -0.46799715716131773, doublethreshold);
    EXPECT_NEAR(verlet->force[2].x, 0.072897549417224428, doublethreshold);
    EXPECT_NEAR(verlet->force[2].y, 0.27055193043529346, doublethreshold);
    EXPECT_NEAR(verlet->force[2].z, -0.12374098810527463, doublethreshold);
    EXPECT_NEAR(verlet->force[3].x, -0.14504625941261556, doublethreshold);
    EXPECT_NEAR(verlet->force[3].y, -0.16104847607824399, doublethreshold);
    EXPECT_NEAR(verlet->force[3].z, 0.17692510112410081, doublethreshold);
}

TEST_F(Langevin_test, first_half)
{
    verlet->first_half();
    
    EXPECT_NEAR(verlet->pos[0].x, 9.9957116654640092, doublethreshold);
    EXPECT_NEAR(verlet->pos[0].y, 0.0016393608896004908, doublethreshold);
    EXPECT_NEAR(verlet->pos[0].z, 0.0049409894499896556, doublethreshold);
    EXPECT_NEAR(verlet->pos[1].x, 5.2079697932877451, doublethreshold);
    EXPECT_NEAR(verlet->pos[1].y, 5.1985329235797453, doublethreshold);
    EXPECT_NEAR(verlet->pos[1].z, 0.0045070523389717319, doublethreshold);
    EXPECT_NEAR(verlet->pos[2].x, 5.0930848914994087, doublethreshold);
    EXPECT_NEAR(verlet->pos[2].y, 0.0011838145470033955, doublethreshold);
    EXPECT_NEAR(verlet->pos[2].z, 4.9932869712840313, doublethreshold);
    EXPECT_NEAR(verlet->pos[3].x, 9.9993642687852322, doublethreshold);
    EXPECT_NEAR(verlet->pos[3].y, 5.2974098983662827, doublethreshold);
    EXPECT_NEAR(verlet->pos[3].z, 5.0029326569457702, doublethreshold);

    EXPECT_NEAR(verlet->vel[0].x, -0.00010372985195918919, doublethreshold);
    EXPECT_NEAR(verlet->vel[0].y, 3.9654243613399205e-05, doublethreshold);
    EXPECT_NEAR(verlet->vel[0].z, 0.00011951681938006538, doublethreshold);
    EXPECT_NEAR(verlet->vel[1].x, 0.00019278008069211768, doublethreshold);
    EXPECT_NEAR(verlet->vel[1].y, -3.5486881587393478e-05, doublethreshold);
    EXPECT_NEAR(verlet->vel[1].z, 0.00010902038261476422, doublethreshold);
    EXPECT_NEAR(verlet->vel[2].x, -0.00016726847568161691, doublethreshold);
    EXPECT_NEAR(verlet->vel[2].y, 2.8635104532351301e-05, doublethreshold);
    EXPECT_NEAR(verlet->vel[2].z, -0.00016238039944434295, doublethreshold);
    EXPECT_NEAR(verlet->vel[3].x, -1.5377602712750055e-05, doublethreshold);
    EXPECT_NEAR(verlet->vel[3].y, -6.2651562458564838e-05, doublethreshold);
    EXPECT_NEAR(verlet->vel[3].z, 7.093757921139429e-05, doublethreshold);
}

TEST_F(Langevin_test, second_half)
{
    verlet->first_half();
    verlet->second_half();
    
    EXPECT_NEAR(verlet->pos[0].x, 9.9933045979909725, doublethreshold);
    EXPECT_NEAR(verlet->pos[0].y, 0.00033862365219131471, doublethreshold);
    EXPECT_NEAR(verlet->pos[0].z, 9.9954281801131337, doublethreshold);
    EXPECT_NEAR(verlet->pos[1].x, 5.1986310958164275, doublethreshold);
    EXPECT_NEAR(verlet->pos[1].y, 5.2027340532086015, doublethreshold);
    EXPECT_NEAR(verlet->pos[1].z, 9.9987348662023798, doublethreshold);
    EXPECT_NEAR(verlet->pos[2].x, 5.0973799076212742, doublethreshold);
    EXPECT_NEAR(verlet->pos[2].y, 0.003868919168865627, doublethreshold);
    EXPECT_NEAR(verlet->pos[2].z, 5.0001845767835942, doublethreshold);
    EXPECT_NEAR(verlet->pos[3].x, 9.9999789728863995, doublethreshold);
    EXPECT_NEAR(verlet->pos[3].y, 5.3031968974372354, doublethreshold);
    EXPECT_NEAR(verlet->pos[3].z, 4.9996952920372832, doublethreshold);
    
    EXPECT_NEAR(verlet->vel[0].x, -0.00011264567421282132, doublethreshold);
    EXPECT_NEAR(verlet->vel[0].y, 7.3217455425250938e-05, doublethreshold);
    EXPECT_NEAR(verlet->vel[0].z, -0.00020991266055077719, doublethreshold);
    EXPECT_NEAR(verlet->vel[1].x, -0.00012678714345051639, doublethreshold);
    EXPECT_NEAR(verlet->vel[1].y, 0.00025321641200267255, doublethreshold);
    EXPECT_NEAR(verlet->vel[1].z, -6.6835251309303343e-05, doublethreshold);
    EXPECT_NEAR(verlet->vel[2].x, 1.0809720799627827e-05, doublethreshold);
    EXPECT_NEAR(verlet->vel[2].y, 0.00025296508641917993, doublethreshold);
    EXPECT_NEAR(verlet->vel[2].z, 1.722136542532692e-07, doublethreshold);
    EXPECT_NEAR(verlet->vel[3].x, -0.00015487279560487432, doublethreshold);
    EXPECT_NEAR(verlet->vel[3].y, 0.00012385428611483992, doublethreshold);
    EXPECT_NEAR(verlet->vel[3].z, 0.00013451165281051973, doublethreshold);
}