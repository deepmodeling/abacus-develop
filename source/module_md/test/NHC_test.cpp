#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "setcell.h"
#include "module_md/Nose_Hoover.h"
#include "module_esolver/esolver_lj.h"

#define doublethreshold 1e-12

/************************************************
 *  unit test of functions in Nose_Hoover.h
 ***********************************************/

/**
 * - Tested Function
 *   - Nose_Hoover::setup
 *     - init before running md, calculate energy, force, and stress of the initial configuration.
 *
 *   - Nose_Hoover::first_half
 *     - the first half of equation of motion, update velocities and positions
 *
 *   - Nose_Hoover::second_half
 *     - the second half of equation of motion, update velocities
 * 
 *   - Nose_Hoover::write_restart
 *     - write the information into files used for MD restarting
 * 
 *   - Nose_Hoover::restart
 *     - restart MD when md_restart is true
 * 
 *   - Nose_Hoover::outputMD
 *     - output MD information such as energy, temperature, and pressure
 */

class NHC_test : public testing::Test
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

        INPUT.mdp.md_pmode = "tri";
        mdrun = new Nose_Hoover(INPUT.mdp, ucell);
        mdrun->setup(p_esolver);
    }

    void TearDown()
    {
        delete mdrun;
    }
};

TEST_F(NHC_test, setup)
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

TEST_F(NHC_test, first_half)
{
    mdrun->first_half();
    
    EXPECT_NEAR(mdrun->pos[0].x, 20.786480524287398, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, -0.17664650841074978, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 12.763557808832209, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 10.417962231096082, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 8.0482704580843034, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0037631015822075175, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 10.326609893752673, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 15.398062039826284, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 6.3799415879158659, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.11327789743365901, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 8.1219079731358974, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 6.380946983836691, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -8.2616217816642025e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 5.3600943884620578e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.1709907223974946e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 8.7469528771814826e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00010625642709930989, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 8.0568608450110721e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -8.5722227110122785e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -4.9154485954659712e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -4.519219457475971e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 8.0868916155623938e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 5.7570931810225902e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.3666506651459607e-05, doublethreshold);
}

TEST_F(NHC_test, second_half)
{
    mdrun->first_half();
    mdrun->second_half();
    
    EXPECT_NEAR(mdrun->pos[0].x, 20.786455241418096, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, -0.1766462896653096, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 12.763542287437854, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 10.417949570368947, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 8.0482606637809138, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0037631021542052384, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, 10.326597330593856, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 15.398043315340637, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 6.3799338269394124, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.11327776681410975, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 8.1218981009532172, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 6.3809392230130584, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -5.6948727377433429e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 4.4909079540334415e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.0924710610576888e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 5.2747513075088019e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -8.5437169183397464e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 7.5166157577425953e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -6.0819449850070865e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -6.6143171452658542e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -4.2161875251693203e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 6.5020664152924765e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 4.7142406788383824e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.2079571715233809e-05, doublethreshold);
}

TEST_F(NHC_test, write_restart)
{
    mdrun->first_half();
    mdrun->second_half();
    mdrun->step_ = 1;
    mdrun->step_rst_ = 2;
    mdrun->write_restart();

    std::ifstream ifs("Restart_md.dat");
    std::string output_str;
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("3"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("4"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("-0.110497   -0.554723   -0.465415   -0.424699   "));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("-0.0116537   -0.0247908   -0.020418   -0.0171934   "));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("0.00391652   0.00331999   0.00198239   -0.000609301   0.000658853   -0.000356508   "));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("4"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("0.537401   2.2182   2.83291   77.7478   "));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr("0.327474   -0.152756   1.07616e-10   3.60329   "));
    ifs.close();
}

TEST_F(NHC_test, restart)
{
    mdrun->restart();
    remove("Restart_md.dat");

    Nose_Hoover* nhc =dynamic_cast<Nose_Hoover*>(mdrun);
    EXPECT_EQ(mdrun->step_rst_, 3);
    EXPECT_EQ(mdrun->mdp.md_tchain, 4);
    EXPECT_EQ(mdrun->mdp.md_pchain, 4);
    EXPECT_EQ(nhc->eta[0], -0.110497);
    EXPECT_EQ(nhc->eta[1], -0.554723);
    EXPECT_EQ(nhc->eta[2], -0.465415);
    EXPECT_EQ(nhc->eta[3], -0.424699);
    EXPECT_EQ(nhc->v_eta[0], -0.0116537);
    EXPECT_EQ(nhc->v_eta[1], -0.0247908);
    EXPECT_EQ(nhc->v_eta[2], -0.020418);
    EXPECT_EQ(nhc->v_eta[3], -0.0171934);
    EXPECT_EQ(nhc->v_omega[0], 0.00391652);
    EXPECT_EQ(nhc->v_omega[1], 0.00331999);
    EXPECT_EQ(nhc->v_omega[2], 0.00198239);
    EXPECT_EQ(nhc->v_omega[3], -0.000609301);
    EXPECT_EQ(nhc->v_omega[4], 0.000658853);
    EXPECT_EQ(nhc->v_omega[5], -0.000356508);
    EXPECT_EQ(nhc->peta[0], 0.537401);
    EXPECT_EQ(nhc->peta[1], 2.2182);
    EXPECT_EQ(nhc->peta[2], 2.83291);
    EXPECT_EQ(nhc->peta[3], 77.7478);
    EXPECT_EQ(nhc->v_peta[0], 0.327474);
    EXPECT_EQ(nhc->v_peta[1], -0.152756);
    EXPECT_EQ(nhc->v_peta[2], 1.07616e-10);
    EXPECT_EQ(nhc->v_peta[3], 3.60329);

}

TEST_F(NHC_test, outputMD)
{
    std::ofstream ofs("running.log");
    mdrun->outputMD(ofs, true);
    ofs.close();

    std::ifstream ifs("running.log");
    std::string output_str;
    getline(ifs,output_str);
    getline(ifs,output_str);
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr(" ------------------------------------------------------------------------------------------------"));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr(" Energy (Ry)         Potential (Ry)      Kinetic (Ry)        Temperature (K)     Pressure (kbar)     "));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr(" -0.015365236        -0.023915637        0.0085504016        300                 1.0846391           "));
    getline(ifs,output_str);
    EXPECT_THAT(output_str,testing::HasSubstr(" ------------------------------------------------------------------------------------------------"));
    ifs.close();
    remove("running.log");
}