#include "gtest/gtest.h"
#include "setcell.h"
#include "module_md/NVT_NHC.h"

class NVT_NHC_test : public testing::Test
{
protected:
    Verlet *verlet;
    UnitCell_pseudo ucell;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();
        verlet = new NVT_NHC(INPUT.mdp, ucell);

        ModuleESolver::ESolver *p_esolver;
        verlet->setup(p_esolver);
    }

    void TearDown()
    {
        delete verlet;
    }
};

TEST_F(NVT_NHC_test, setup)
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

TEST_F(NVT_NHC_test, first_half)
{
    verlet->first_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.9945454470992772);
    EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.0029590658162135398);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 9.9994204767196599);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2063192793031234);
    EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1939342598421803);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 0.0039873402383468932);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0944648036273872);
    EXPECT_DOUBLE_EQ(verlet->pos[2].y, 9.9998836061438716);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.997763438399228);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.0046704699702541487);
    EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.303223068197739);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9988287446427613);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00013193932519649489);
    EXPECT_DOUBLE_EQ(verlet->vel[0].y, 7.157637923935656e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.4017997796600018e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.0001528560566112948);
    EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.00014672323796402804);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 9.6449148069800125e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00013388999749840022);
    EXPECT_DOUBLE_EQ(verlet->vel[2].y, -2.8154327428808196e-06);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -5.4099838013700078e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 0.00011297326608460035);
    EXPECT_DOUBLE_EQ(verlet->vel[3].y, 7.7962291467652311e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -2.8331312259600037e-05);
}

TEST_F(NVT_NHC_test, second_half)
{
    verlet->first_half();
    verlet->second_half();
    
    EXPECT_DOUBLE_EQ(verlet->pos[0].x, 9.9945454470992772);
    EXPECT_DOUBLE_EQ(verlet->pos[0].y, 0.0029590658162135398);
    EXPECT_DOUBLE_EQ(verlet->pos[0].z, 9.9994204767196599);
    EXPECT_DOUBLE_EQ(verlet->pos[1].x, 5.2063192793031234);
    EXPECT_DOUBLE_EQ(verlet->pos[1].y, 5.1939342598421803);
    EXPECT_DOUBLE_EQ(verlet->pos[1].z, 0.0039873402383468932);
    EXPECT_DOUBLE_EQ(verlet->pos[2].x, 5.0944648036273872);
    EXPECT_DOUBLE_EQ(verlet->pos[2].y, 9.9998836061438716);
    EXPECT_DOUBLE_EQ(verlet->pos[2].z, 4.997763438399228);
    EXPECT_DOUBLE_EQ(verlet->pos[3].x, 0.0046704699702541487);
    EXPECT_DOUBLE_EQ(verlet->pos[3].y, 5.303223068197739);
    EXPECT_DOUBLE_EQ(verlet->pos[3].z, 4.9988287446427613);

    EXPECT_DOUBLE_EQ(verlet->vel[0].x, -0.00013179592306188965);
    EXPECT_DOUBLE_EQ(verlet->vel[0].y, 7.1808730720987585e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[0].z, -1.4017786037767368e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[1].x, 0.00015266992839159681);
    EXPECT_DOUBLE_EQ(verlet->vel[1].y, -0.00014691099001811579);
    EXPECT_DOUBLE_EQ(verlet->vel[1].z, 9.6447691088617668e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[2].x, -0.00013398849068181084);
    EXPECT_DOUBLE_EQ(verlet->vel[2].y, -2.585722613500158e-06);
    EXPECT_DOUBLE_EQ(verlet->vel[2].z, -5.4099020770131449e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].x, 0.0001131144853531037);
    EXPECT_DOUBLE_EQ(verlet->vel[3].y, 7.7687981910728375e-05);
    EXPECT_DOUBLE_EQ(verlet->vel[3].z, -2.8330884280818859e-05);
}