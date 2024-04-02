#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "module_esolver/esolver_lj.h"
#include "setcell.h"

#define private public
#define protected public
#include "module_md/verlet.h"

#define doublethreshold 1e-12

/************************************************
 *  unit test of functions in verlet.h
 ***********************************************/

/**
 * - Tested Function
 *   - verlet::setup
 *     - init before running md, calculate energy, force, and stress of the initial configuration.
 *
 *   - verlet::first_half
 *     - the first half of equation of motion, update velocities and positions
 *
 *   - verlet::second_half
 *     - the second half of equation of motion, update velocities
 *
 *   - verlet::write_restart
 *     - write the information into files used for MD restarting
 *
 *   - verlet::restart
 *     - restart MD when md_restart is true
 *
 *   - verlet::print_md
 *     - output MD information such as energy, temperature, and pressure
 */

class Verlet_test : public testing::Test
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

        mdrun = new Verlet(INPUT.mdp, ucell);
        mdrun->setup(p_esolver, GlobalV::global_readin_dir);
    }

    void TearDown()
    {
        delete mdrun;
    }
};

TEST_F(Verlet_test, setup)
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

TEST_F(Verlet_test, first_half)
{
    mdrun->first_half(GlobalV::ofs_running);

    EXPECT_NEAR(mdrun->pos[0].x, -0.00054545529007222658, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -5.7952328034033513e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00063192793031220879, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00060657401578200095, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.00039873402383468892, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00055351963726126224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, -1.1639385612741475e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00022365616007718661, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00046704699702541431, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00032230681977380224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00011712553572388214, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013193932519649473, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1576379239356465e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.40179977966e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015285605661129458, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014672323796402785, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6449148069800003e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013388999749840003, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.8154327428808153e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4099838013700003e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.0001129732660846002, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7962291467652202e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.83313122596e-05, doublethreshold);
}

TEST_F(Verlet_test, NVE)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->mdp.md_type = "nve";
    mdrun->second_half();
    ;

    EXPECT_NEAR(mdrun->pos[0].x, -0.00054545529007222658, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -5.7952328034033513e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00063192793031220879, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00060657401578200095, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.00039873402383468892, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00055351963726126224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, -1.1639385612741475e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00022365616007718661, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00046704699702541431, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00032230681977380224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00011712553572388214, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013179791402898947, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1809815495212933e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.40179977966e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015267223469058917, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.00014691320932005571, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6449148069800003e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013399051477080008, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5857616745616307e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4099838013700003e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.0001131161941102004, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7689155499504408e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.83313122596e-05, doublethreshold);
}

TEST_F(Verlet_test, Anderson)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->mdp.md_type = "nvt";
    mdrun->mdp.md_thermostat = "anderson";
    mdrun->second_half();
    ;

    EXPECT_NEAR(mdrun->pos[0].x, -0.00054545529007222658, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -5.7952328034033513e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00063192793031220879, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00060657401578200095, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.00039873402383468892, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00055351963726126224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, -1.1639385612741475e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00022365616007718661, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00046704699702541431, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00032230681977380224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00011712553572388214, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013179791402898947, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1809815495212933e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.40179977966e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 6.9452562329904563e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, 7.321611395307015e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, -8.133446733603267e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, 0.00013239881096711222, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, 0.00030862680563211305, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -0.00012925479702246553, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.0001131161941102004, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7689155499504408e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.83313122596e-05, doublethreshold);
}

TEST_F(Verlet_test, Berendsen)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->mdp.md_type = "nvt";
    mdrun->mdp.md_thermostat = "berendsen";
    mdrun->second_half();
    ;

    EXPECT_NEAR(mdrun->pos[0].x, -0.00054545529007222658, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -5.7952328034033513e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00063192793031220879, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00060657401578200095, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.00039873402383468892, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00055351963726126224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, -1.1639385612741475e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00022365616007718661, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00046704699702541431, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00032230681977380224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00011712553572388214, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013179175250738632, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1806458403162173e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.4017342458487154e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015266509729938552, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.0001469063411619389, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6444639094723906e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013398425074562592, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5856407908091386e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4097308858947404e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.00011311090595462667, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7685523549685863e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8329987777389342e-05, doublethreshold);
}

TEST_F(Verlet_test, rescaling)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->mdp.md_type = "nvt";
    mdrun->mdp.md_thermostat = "rescaling";
    mdrun->second_half();
    ;

    EXPECT_NEAR(mdrun->pos[0].x, -0.00054545529007222658, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -5.7952328034033513e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00063192793031220879, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00060657401578200095, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.00039873402383468892, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00055351963726126224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, -1.1639385612741475e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00022365616007718661, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00046704699702541431, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00032230681977380224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00011712553572388214, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013178559069770653, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1803101154153484e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.4016687089734539e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015265795957447931, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.0001468994726827073, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6440129908834563e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013397798642758268, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5855199014048311e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4094779585946356e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.0001131056175518098, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7681891430058639e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8328663233253657e-05, doublethreshold);
}

TEST_F(Verlet_test, rescale_v)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->mdp.md_type = "nvt";
    mdrun->mdp.md_thermostat = "rescale_v";
    mdrun->second_half();
    ;

    EXPECT_NEAR(mdrun->pos[0].x, -0.00054545529007222658, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].y, 0.00029590658162135359, doublethreshold);
    EXPECT_NEAR(mdrun->pos[0].z, -5.7952328034033513e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].x, 0.00063192793031220879, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].y, -0.00060657401578200095, doublethreshold);
    EXPECT_NEAR(mdrun->pos[1].z, 0.00039873402383468892, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].x, -0.00055351963726126224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].y, -1.1639385612741475e-05, doublethreshold);
    EXPECT_NEAR(mdrun->pos[2].z, -0.00022365616007718661, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].x, 0.00046704699702541431, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].y, 0.00032230681977380224, doublethreshold);
    EXPECT_NEAR(mdrun->pos[3].z, -0.00011712553572388214, doublethreshold);

    EXPECT_NEAR(mdrun->vel[0].x, -0.00013178559069770653, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].y, 7.1803101154153484e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[0].z, -1.4016687089734539e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].x, 0.00015265795957447931, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].y, -0.0001468994726827073, doublethreshold);
    EXPECT_NEAR(mdrun->vel[1].z, 9.6440129908834563e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].x, -0.00013397798642758268, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].y, -2.5855199014048311e-06, doublethreshold);
    EXPECT_NEAR(mdrun->vel[2].z, -5.4094779585946356e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].x, 0.0001131056175518098, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].y, 7.7681891430058639e-05, doublethreshold);
    EXPECT_NEAR(mdrun->vel[3].z, -2.8328663233253657e-05, doublethreshold);
}

TEST_F(Verlet_test, write_restart)
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

TEST_F(Verlet_test, restart)
{
    mdrun->restart(GlobalV::global_readin_dir);
    remove("Restart_md.dat");

    EXPECT_EQ(mdrun->step_rst_, 3);
}

TEST_F(Verlet_test, print_md)
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
