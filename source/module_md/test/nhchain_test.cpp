#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "module_parameter/parameter.h"
#undef private
#define private public
#define protected public
#include "module_esolver/esolver_lj.h"
#include "module_md/nhchain.h"
#include "setcell.h"
#define doublethreshold 1e-12
/************************************************
 *  unit test of functions in nhchain.h
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
 *   - Nose_Hoover::print_md
 *     - output MD information such as energy, temperature, and pressure
 */
class NHC_test : public testing::Test
{
  protected:
    MD_base* mdrun;
    UnitCell ucell;
    Parameter param_in;
    ModuleESolver::ESolver* p_esolver;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters(param_in.input);

        p_esolver = new ModuleESolver::ESolver_LJ();
        p_esolver->before_all_runners(ucell, param_in.inp);

        param_in.input.mdp.md_type = "npt";
        param_in.input.mdp.md_pmode = "tri";
        param_in.input.mdp.md_pfirst = 1;
        param_in.input.mdp.md_plast = 1;
        mdrun = new Nose_Hoover(param_in, ucell);
        mdrun->setup(p_esolver, PARAM.sys.global_readin_dir);
    }

    void TearDown()
    {
        delete mdrun;
        delete p_esolver;
    }
};

TEST_F(NHC_test, setup)
{
    EXPECT_NEAR(mdrun->t_current * ModuleBase::Hartree_to_K, 299.99999999999994, doublethreshold);
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

TEST_F(NHC_test, first_half)
{
    mdrun->first_half(GlobalV::ofs_running);

    EXPECT_NEAR(mdrun->pos[0].x, -0.00035596392702161582, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00026566987683715606, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -6.4082739615824722e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00037007414441809518, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00052501803299631633, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.00044091358349508534, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00036876922955593201, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, -2.6151466573228018e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00024731533582713971, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00035465901216238645, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00028549962273273618, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00012951550805257814, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00010335325828338315, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 6.6973537793984337e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.4644123959592966e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00010943331752057692, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00013283409023334643, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 0.00010075713383789103, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00010717693628353973, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -6.2046899135633754e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.6516254714969195e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00010109687704718878, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.2065242353013738e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.9596755163433345e-05, doublethreshold);
}

TEST_F(NHC_test, second_half)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->second_half();
    ;

    EXPECT_NEAR(mdrun->pos[0].x, -0.00035596392702161582, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00026566987683715606, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -6.4082739615824722e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00037007414441809518, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00052501803299631633, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.00044091358349508534, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00036876922955593201, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, -2.6151466573228018e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00024731533582713971, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00035465901216238645, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00028549962273273618, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00012951550805257814, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -8.4972683205367143e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 6.6834262571392232e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.6287026488367857e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 7.8726485842843947e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00012727726730227848, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 0.00011206092711573642, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -9.0636235945876312e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -9.9771188254262979e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -6.285672943672849e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 9.6882433309157637e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.0420123556394411e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -3.2917171190756263e-05, doublethreshold);
}

TEST_F(NHC_test, write_restart)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->second_half();
    ;
    mdrun->step_ = 1;
    mdrun->step_rst_ = 2;
    mdrun->write_restart(PARAM.sys.global_out_dir);

    std::ifstream ifs("Restart_md.dat");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("3"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("0.000950045"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("4"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("-0.0626326   -0.578523   -0.462472   -0.424503   "));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("-0.00658882   -0.0304055   -0.0188618   -0.0175663   "));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("0.583152   -0.106519   -0.895936   -0.634424   0.627532   -0.473422   "));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("4"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("-6.08823   -0.525329   0.121814   4771.79   "));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("255.853   -0.266732   0   226.197   "));
    ifs.close();
}

TEST_F(NHC_test, restart)
{
    mdrun->restart(PARAM.sys.global_readin_dir);
    remove("Restart_md.dat");

    Nose_Hoover* nhc = dynamic_cast<Nose_Hoover*>(mdrun);
    EXPECT_EQ(mdrun->step_rst_, 3);
    EXPECT_EQ(mdrun->mdp.md_tchain, 4);
    EXPECT_EQ(mdrun->mdp.md_pchain, 4);
    EXPECT_EQ(nhc->eta[0], -0.0626326);
    EXPECT_EQ(nhc->eta[1], -0.578523);
    EXPECT_EQ(nhc->eta[2], -0.462472);
    EXPECT_EQ(nhc->eta[3], -0.424503);
    EXPECT_EQ(nhc->v_eta[0], -0.00658882);
    EXPECT_EQ(nhc->v_eta[1], -0.0304055);
    EXPECT_EQ(nhc->v_eta[2], -0.0188618);
    EXPECT_EQ(nhc->v_eta[3], -0.0175663);
    EXPECT_EQ(nhc->v_omega[0], 0.583152);
    EXPECT_EQ(nhc->v_omega[1], -0.106519);
    EXPECT_EQ(nhc->v_omega[2], -0.895936);
    EXPECT_EQ(nhc->v_omega[3], -0.634424);
    EXPECT_EQ(nhc->v_omega[4], 0.627532);
    EXPECT_EQ(nhc->v_omega[5], -0.473422);
    EXPECT_EQ(nhc->peta[0], -6.08823);
    EXPECT_EQ(nhc->peta[1], -0.525329);
    EXPECT_EQ(nhc->peta[2], 0.121814);
    EXPECT_EQ(nhc->peta[3], 4771.79);
    EXPECT_EQ(nhc->v_peta[0], 255.853);
    EXPECT_EQ(nhc->v_peta[1], -0.266732);
    EXPECT_EQ(nhc->v_peta[2], 0);
    EXPECT_EQ(nhc->v_peta[3], 226.197);
}

TEST_F(NHC_test, print_md)
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
