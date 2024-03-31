#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_esolver/esolver_lj.h"
#include "setcell.h"

#define private public
#define protected public
#include "module_md/langevin.h"

#define doublethreshold 1e-12

/************************************************
 *  unit test of functions in langevin.h
 ***********************************************/

/**
 * - Tested Function
 *   - Langevin::setup
 *     - init before running md, calculate energy, force, and stress of the initial configuration.
 *
 *   - Langevin::first_half
 *     - the first half of equation of motion, update velocities and positions
 *
 *   - Langevin::second_half
 *     - the second half of equation of motion, update velocities
 *
 *   - Langevin::write_restart
 *     - write the information into files used for MD restarting
 *
 *   - Langevin::restart
 *     - restart MD when md_restart is true
 *
 *   - Langevin::print_md
 *     - output MD information such as energy, temperature, and pressure
 */

class Langevin_test : public testing::Test
{
  protected:
    MD_base* mdrun;
    UnitCell ucell;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters();

        ModuleESolver::ESolver* p_esolver = new ModuleESolver::ESolver_LJ();
        p_esolver->init(INPUT, ucell);

        mdrun = new Langevin(INPUT.mdp, ucell);
        mdrun->setup(p_esolver, GlobalV::global_readin_dir);
    }

    void TearDown()
    {
        delete mdrun;
    }
};

TEST_F(Langevin_test, setup)
{
    EXPECT_NEAR(mdrun->t_current * ModuleBase::Hartree_to_K, 299.99999999999665, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0, 0), 6.0100555286436806e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0, 1), -1.4746713013791574e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(0, 2), 1.5039983732220751e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1, 0), -1.4746713013791574e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1, 1), 3.4437172989317909e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(1, 2), -1.251414906590483e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2, 0), 1.5039983732220751e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2, 1), -1.251414906590483e-06, doublethreshold);
    EXPECT_NEAR(mdrun->stress(2, 2), 1.6060561926126463e-06, doublethreshold);
}

TEST_F(Langevin_test, first_half)
{
    mdrun->first_half(GlobalV::ofs_running);

    EXPECT_NEAR(mdrun->pos[0].x, -0.00042883345359910814, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00016393608896004904, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, 0.00049409894499896569, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00079697932877452634, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00014670764202547791, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.0004507052338971732, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00069151085005912606, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 0.00011838145470033956, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00067130287159685429, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, -6.3573121476911994e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, -0.00025901016337184232, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, 0.00029326569457701463, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00010372985195918919, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 3.9654243613399205e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, 0.00011951681938006538, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00019278008069211768, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -3.5486881587393478e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 0.00010902038261476422, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00016726847568161691, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, 2.8635104532351301e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -0.00016238039944434295, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, -1.5377602712750055e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, -6.2651562458564838e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, 7.093757921139429e-05, doublethreshold);
}

TEST_F(Langevin_test, second_half)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->second_half();
    ;

    EXPECT_NEAR(mdrun->pos[0].x, -0.00066954020090275205, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 3.3862365219131354e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -0.00045718198868662484, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, -0.0001368904183573199, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, 0.00027340532086011393, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, -0.00012651337976204397, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00026200923787255071, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 0.00038689191688656276, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, 1.8457678359430833e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, -2.1027113600492346e-06, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.000319689743723507, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -3.0470796271690045e-05, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -8.2630969616448438e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 0.0001366029202159129, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -0.00011334362366793093, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 5.9181121902101574e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, 4.0359589497484719e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 6.0216019900454962e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -5.9703272809828887e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, 0.00015656497429546092, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.8392323248176516e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, -4.1390907965075468e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 0.00012448732653877297, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, 0.00011355087370269158, doublethreshold);
}

TEST_F(Langevin_test, write_restart)
{
    mdrun->step_ = 1;
    mdrun->step_rst_ = 2;
    mdrun->write_restart(GlobalV::global_out_dir);

    std::ifstream ifs("Restart_md.dat");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("3"));
    ifs.close();
}

TEST_F(Langevin_test, restart)
{
    mdrun->restart(GlobalV::global_readin_dir);
    remove("Restart_md.dat");

    EXPECT_EQ(mdrun->step_rst_, 3);
}

TEST_F(Langevin_test, print_md)
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
            " -0.015365236        -0.023915637        0.0085504016        300                 1.0846391           "));
    getline(ifs, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr(
            " ------------------------------------------------------------------------------------------------"));
    ifs.close();
    remove("running.log");
}
