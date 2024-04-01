#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_esolver/esolver_lj.h"
#include "setcell.h"

#define private public
#define protected public
#include "module_md/fire.h"

#define doublethreshold 1e-12

/************************************************
 *  unit test of functions in fire.h
 ***********************************************/

/**
 * - Tested Function
 *   - FIRE::setup
 *     - init before running md, calculate energy, force, and stress of the initial configuration.
 *
 *   - FIRE::first_half
 *     - the first half of equation of motion, update velocities and positions
 *
 *   - FIRE::second_half
 *     - the second half of equation of motion, update velocities
 *
 *   - FIRE::write_restart
 *     - write the information into files used for MD restarting
 *
 *   - FIRE::restart
 *     - restart MD when md_restart is true
 *
 *   - FIRE::print_md
 *     - output MD information such as energy, temperature, and pressure
 */

class FIREtest : public testing::Test
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

        mdrun = new FIRE(INPUT.mdp, ucell);
        mdrun->setup(p_esolver, GlobalV::global_readin_dir);
    }

    void TearDown()
    {
        delete mdrun;
    }
};

TEST_F(FIREtest, Setup)
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

TEST_F(FIREtest, FirstHalf)
{
    mdrun->first_half(GlobalV::ofs_running);

    EXPECT_NEAR(mdrun->pos[0].x, -0.00045447059554315662, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00032646833232493271, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -5.215709523063016e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.0005213674681407162, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00059486888444406608, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.00035886062145122004, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00052406920303529794, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 4.8706739346586155e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00020129054406946794, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00045717233044145918, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00021969381277291936, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00010541298215149392, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00010993118004167345, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.8968913216100539e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.2616198016939999e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00012611275970351733, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014389190209072655, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 8.6804233262820007e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00012676627812260489, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, 1.1781596840062159e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -4.8689854212330001e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011058469846166102, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 5.3141392034653857e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.5498181033639999e-05, doublethreshold);
}

TEST_F(FIREtest, SecondHalf)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->second_half();

    EXPECT_NEAR(mdrun->pos[0].x, -0.00045447059554315662, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00032646833232493271, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -5.215709523063016e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.0005213674681407162, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00059486888444406608, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.00035886062145122004, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00052406920303529794, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, 4.8706739346586155e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00020129054406946794, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00045717233044145918, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00021969381277291936, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00010541298215149392, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00010978976887416819, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.9202349471957007e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.2616198016939999e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00012592893778281191, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014408187344675441, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 8.6804233262820007e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00012686679539500493, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, 1.2011267908381344e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -4.8689854212330001e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011072762648726122, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 5.2868256066506055e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.5498181033639999e-05, doublethreshold);
}

TEST_F(FIREtest, WriteRestart)
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
    EXPECT_THAT(output_str, testing::HasSubstr("0.1"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("0"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("-1"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("41.3414"));
    ifs.close();
}

TEST_F(FIREtest, Restart)
{
    mdrun->restart(GlobalV::global_readin_dir);
    remove("Restart_md.dat");

    FIRE* fire = dynamic_cast<FIRE*>(mdrun);
    EXPECT_EQ(mdrun->step_rst_, 3);
    EXPECT_EQ(fire->alpha, 0.1);
    EXPECT_EQ(fire->negative_count, 0);
    EXPECT_EQ(fire->dt_max, -1);
    EXPECT_EQ(mdrun->mdp.md_dt, 41.3414);
}

TEST_F(FIREtest, PrintMD)
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
    for (int i = 0; i < 10; ++i)
    {
        getline(ifs, output_str);
    }
    EXPECT_THAT(output_str, testing::HasSubstr(" LARGEST GRAD (eV/A)  : 0.049479926"));
    ifs.close();
    remove("running.log");
}
