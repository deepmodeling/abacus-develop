#include "gtest/gtest.h"
#include "setcell.h"
#include "module_md/NVT_ADS.h"

class NVT_ADS_test : public testing::Test
{
protected:
    Verlet *verlet;
    UnitCell_pseudo ucell;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();
        verlet = new NVT_ADS(INPUT.mdp, ucell);

        ModuleESolver::ESolver *p_esolver;
        verlet->setup(p_esolver);
    }

    void TearDown()
    {
        delete verlet;
    }
};

TEST_F(NVT_ADS_test, setup)
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

TEST_F(NVT_ADS_test, first_half)
{
    verlet->first_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.9945454470992772);
    EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.0029590658162135359);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 9.9994204767196599);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2063192793031234);
    EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1939342598421803);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 0.0039873402383468889);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0944648036273872);
    EXPECT_DOUBLE_EQ(verlet->pos[2].y, 9.9998836061438716);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.997763438399228);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.0046704699702541427);
    EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.303223068197739);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9988287446427613);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00013193932519649473);
    EXPECT_DOUBLE_EQ(verlet->vel[0].y, 7.1576379239356465e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.40179977966e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.00015285605661129458);
    EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.00014672323796402785);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 9.6449148069800003e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00013388999749840003);
    EXPECT_DOUBLE_EQ(verlet->vel[2].y, -2.8154327428808153e-06);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -5.4099838013700003e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 0.0001129732660846002);
    EXPECT_DOUBLE_EQ(verlet->vel[3].y, 7.7962291467652202e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -2.83313122596e-05);
}

TEST_F(NVT_ADS_test, second_half)
{
    verlet->first_half();
    verlet->second_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.9945454470992772);
    EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.0029590658162135359);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 9.9994204767196599);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2063192793031234);
    EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1939342598421803);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 0.0039873402383468889);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0944648036273872);
    EXPECT_DOUBLE_EQ(verlet->pos[2].y, 9.9998836061438716);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.997763438399228);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.0046704699702541427);
    EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.303223068197739);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9988287446427613);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -5.6668160441304611e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].y, 0.00015189506534878691);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -0.00012247221187246948);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, -6.6724169081868237e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.0001963389480223568);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 4.768316146189508e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, 0.0001563891187049906);
    EXPECT_DOUBLE_EQ(verlet->vel[2].y, -0.00016408850308461561);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, 2.1395048951018702e-06);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 7.2200148840519235e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].y, 9.8952839341203403e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, 0.00015843986714134529);
}