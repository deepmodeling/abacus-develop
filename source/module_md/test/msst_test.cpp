#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_esolver/esolver_lj.h"
#include "setcell.h"

#define private public
#define protected public
#include "module_md/msst.h"

#define doublethreshold 1e-12

/************************************************
 *  unit test of functions in msst.h
 ***********************************************/

/**
 * - Tested Function
 *   - MSST::setup
 *     - init before running md, calculate energy, force, and stress of the initial configuration.
 *
 *   - MSST::first_half
 *     - the first half of equation of motion, update velocities and positions
 *
 *   - MSST::second_half
 *     - the second half of equation of motion, update velocities
 *
 *   - MSST::write_restart
 *     - write the information into files used for MD restarting
 *
 *   - MSST::restart
 *     - restart MD when md_restart is true
 *
 *   - MSST::print_md
 *     - output MD information such as energy, temperature, and pressure
 */

class MSST_test : public testing::Test
{
  protected:
    MD_base* mdrun;
    UnitCell ucell;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();

        ModuleESolver::ESolver* p_esolver = new ModuleESolver::ESolver_LJ();
        p_esolver->Init(INPUT, ucell);

        mdrun = new MSST(INPUT.mdp, ucell);
        mdrun->setup(p_esolver, GlobalV::global_readin_dir);
    }

    void TearDown()
    {
        delete mdrun;
    }
};

TEST_F(MSST_test, setup)
{
    EXPECT_NEAR(mdrun->vel[0].x, -0.0001314186733659715, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.0985331994796372e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.3947731701005279e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015227275651566311, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014579875939315496, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.5965690649087203e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013311885204189453, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -3.0298400368294885e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.3828659173134662e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011226476889319793, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7843267435287586e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8189299775046767e-05, doublethreshold);

    EXPECT_NEAR(mdrun->stress(0, 0), 5.9579909955800075e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0, 1), -1.4582038138067117e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0, 2), 1.4889583894898544e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1, 0), -1.4582038138067117e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1, 1), 3.4199108345556597e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1, 2), -1.2389007575245785e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2, 0), 1.4889583894898544e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2, 1), -1.2389007575245785e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2, 2), 1.5964231736437884e-06, doublethreshold);
}

TEST_F(MSST_test, first_half)
{
    mdrun->first_half(GlobalV::ofs_running);

    EXPECT_NEAR(ucell.lat0, 1.0, doublethreshold);
    EXPECT_NEAR(ucell.lat0_angstrom, 0.52917700000000001, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e11, 10.0, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e12, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e13, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e21, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e22, 10.0, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e23, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e31, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e32, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e33, 9.9959581179144905, doublethreshold);
    EXPECT_NEAR(ucell.omega, 999.59581179144902, doublethreshold);

    EXPECT_NEAR(mdrun->pos[0].x, -0.00054271823071484467, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00029442816868202821, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -5.7685149290774873e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00062875654254500096, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00060353746208327032, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0003968957326219519, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00055074716824834991, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, -1.1576283073263842e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00022262503373940808, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00046470885642230719, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00032068557647491725, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00011658554959218052, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013127726219846624, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.121876825065284e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.3947730561390963e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.0001520889345949577, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014598873074918282, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.596568280810794e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013321936931429457, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.8001689685103039e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.3828654775006574e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011240769691879813, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7570131467139791e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8189297471809918e-05, doublethreshold);
}

TEST_F(MSST_test, second_half)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->second_half();
    ;

    EXPECT_NEAR(ucell.lat0, 1.0, doublethreshold);
    EXPECT_NEAR(ucell.lat0_angstrom, 0.52917700000000001, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e11, 10.0, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e12, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e13, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e21, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e22, 10.0, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e23, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e31, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e32, 0.00, doublethreshold);
    EXPECT_NEAR(ucell.latvec.e33, 9.9959581179144905, doublethreshold);
    EXPECT_NEAR(ucell.omega, 999.59581179144902, doublethreshold);

    EXPECT_NEAR(mdrun->pos[0].x, -0.00054271823071484467, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00029442816868202821, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -5.7685149290774873e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00062875654254500096, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00060353746208327032, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0003968957326219519, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00055074716824834991, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, -1.1576283073263842e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00022262503373940808, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00046470885642230719, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00032068557647491725, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00011658554959218052, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013113585103096098, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1452204506509308e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.3953371489538059e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015190511267425228, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014617870210521068, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.600449453585996e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013331988658669462, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5704979001911192e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.3850424881082548e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011255062494439833, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7296995498991997e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8200698165338931e-05, doublethreshold);
}

TEST_F(MSST_test, write_restart)
{
    mdrun->step_ = 1;
    mdrun->step_rst_ = 2;
    mdrun->write_restart(GlobalV::global_out_dir);

    std::ifstream ifs("Restart_md.dat");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("3"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("0.000950045"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("-0.00977662"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("-0.00768262"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("1000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("1.60606e-06"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("0"));
    ifs.close();
}

TEST_F(MSST_test, restart)
{
    mdrun->restart(GlobalV::global_readin_dir);
    remove("Restart_md.dat");

    MSST* msst = dynamic_cast<MSST*>(mdrun);
    EXPECT_EQ(mdrun->step_rst_, 3);
    EXPECT_EQ(msst->omega[mdrun->mdp.msst_direction], -0.00977662);
    EXPECT_EQ(msst->e0, -0.00768262);
    EXPECT_EQ(msst->v0, 1000);
    EXPECT_EQ(msst->p0, 1.60606e-06);
    EXPECT_EQ(msst->lag_pos, 0);
}

TEST_F(MSST_test, print_md)
{
    std::ofstream ofs("running.log");
    mdrun->print_md(ofs, true);
    ofs.close();

    std::ifstream ifs("running.log");
    std::string output_str;
    getline(ifs, output_str);
    getline(ifs, output_str);
    getline(ifs, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr(
            " ------------------------------------------------------------------------------------------------"));
    getline(ifs, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr(
            " Energy (Ry)         Potential (Ry)      Kinetic (Ry)        Temperature (K)     Pressure (kbar)     "));
    getline(ifs, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr(
            " -0.01545074         -0.023915637        0.0084648976        297                 1.0762537           "));
    getline(ifs, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr(
            " ------------------------------------------------------------------------------------------------"));
    ifs.close();
    remove("running.log");
}
