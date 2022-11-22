#include "gtest/gtest.h"
#include "setcell.h"
#include "module_md/Nose_Hoover.h"
#include "module_esolver/esolver_lj.h"

#define doublethreshold 1e-12

class NVT_NHC_test : public testing::Test
{
protected:
    MDrun *mdrun;
    UnitCell ucell;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();

        ModuleESolver::ESolver *p_esolver = new ModuleESolver::ESolver_LJ();
        p_esolver->Init(INPUT, ucell);

        mdrun = new Nose_Hoover(INPUT.mdp, ucell);
        mdrun->setup(p_esolver);
    }

    void TearDown()
    {
        delete mdrun;
    }
};

TEST_F(NVT_NHC_test, setup)
{
    EXPECT_NEAR(mdrun->t_current * ModuleBase::Hartree_to_K, 299.99999999999665, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0,0), 6.0100555286436806e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0,1), -1.4746713013791574e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0,2), 1.5039983732220751e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1,0), -1.4746713013791574e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1,1), 3.4437172989317909e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1,2), -1.251414906590483e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2,0), 1.5039983732220751e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2,1), -1.251414906590483e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2,2), 1.6060561926126463e-06, doublethreshold);
}

TEST_F(NVT_NHC_test, first_half)
{
    mdrun->first_half();
    
    EXPECT_NEAR(mdrun->pos[0].x, 9.9945454470992772, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.0029590658162135398, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 9.9994204767196599, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 5.2063192793031234, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 5.1939342598421803, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0039873402383468932, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 5.0944648036273872, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 9.9998836061438716, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 4.997763438399228, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.0046704699702541487, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 5.303223068197739, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 4.9988287446427613, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013193932519649489, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.157637923935656e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.4017997796600018e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.0001528560566112948, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014672323796402804, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6449148069800125e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013388999749840022, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.8154327428808196e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4099838013700078e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011297326608460035, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7962291467652311e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8331312259600037e-05, doublethreshold);
}

TEST_F(NVT_NHC_test, second_half)
{
    mdrun->first_half();
    mdrun->second_half();
    
    EXPECT_NEAR(mdrun->pos[0].x, 9.9945454470992772, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.0029590658162135398, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 9.9994204767196599, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 5.2063192793031234, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 5.1939342598421803, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0039873402383468932, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 5.0944648036273872, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 9.9998836061438716, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 4.997763438399228, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.0046704699702541487, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 5.303223068197739, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 4.9988287446427613, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013179444804958438, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1807927063643868e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.4017629155775739e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015266821976497321, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014690934584344768, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6446611680940509e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013398699113108172, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5856936750278264e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4098415313456568e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011311321941669288, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7687112454931641e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8330567211808204e-05, doublethreshold);
}