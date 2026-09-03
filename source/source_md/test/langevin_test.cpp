#include "gmock/gmock.h"
#include "gtest/gtest.h"
#define private public
#include "source_io/module_parameter/parameter.h"
#undef private
#define private public
#define protected public
#include "setcell.h"
#include "source_esolver/esolver_lj.h"
#include "source_md/langevin.h"
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
    MDCell* mdcell;
    Parameter param_in;
    ModuleESolver::ESolver* p_esolver;

    void SetUp()
    {
        Setcell::setupcell(ucell);
        Setcell::parameters(param_in.input);

        p_esolver = new ModuleESolver::ESolver_LJ();
        mdcell = new MDCell;
        mdcell->initialize_from_unitcell(ucell,
                                         8.5 * ModuleBase::ANGSTROM_AU,
                                         0.0,
                                         ModuleBase::world_communication_domain());
        p_esolver->before_all_runners(*mdcell, param_in.inp);
        mdrun = new Langevin(param_in, *mdcell);
        mdrun->setup(p_esolver, PARAM.sys.global_readin_dir);
    }

    void TearDown()
    {
        delete mdrun;
        delete mdcell;
        delete p_esolver;
    }
};

TEST_F(Langevin_test, setup)
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

TEST_F(Langevin_test, first_half)
{
    mdrun->first_half(GlobalV::ofs_running);

    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(0)]).x, 0.00012104549072633688, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(0)]).y, 2.6272991724490339e-05, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(0)]).z, 0.0002984728051383459, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(1)]).x, 0.00066077703137157329, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(1)]).y, 0.00017245549939737259, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(1)]).z, -0.00015046260270490386, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(2)]).x, -0.00046755850510571406, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(2)]).y, 0.00030490494761200812, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(2)]).z, -0.00036886672854369307, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(3)]).x, 0.00029624371924650601, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(3)]).y, 0.00013444493932002199, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(3)]).z, 9.0496812405138627e-05, doublethreshold);

    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(0)].vel.x, 2.9279504031204254e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(0)].vel.y, 6.3551327892765255e-06, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(0)].vel.z, 7.2197119023712585e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(1)].vel.x, 0.00015983432044991118, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(1)].vel.y, 4.1714990451188365e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(1)].vel.z, -3.639516314080526e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(2)].vel.x, -0.0001130969939724264, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(2)].vel.y, 7.3752979885251439e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(2)].vel.z, -8.9224594824346108e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(3)].vel.x, 7.1657928931092104e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(3)].vel.y, 3.2520675649909766e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(3)].vel.z, 2.189013211254232e-05, doublethreshold);
}

TEST_F(Langevin_test, second_half)
{
    mdrun->first_half(GlobalV::ofs_running);
    mdrun->second_half();
    ;

    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(0)]).x, 0.00012104549072633688, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(0)]).y, 2.6272991724490339e-05, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(0)]).z, 0.0002984728051383459, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(1)]).x, 0.00066077703137157329, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(1)]).y, 0.00017245549939737259, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(1)]).z, -0.00015046260270490386, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(2)]).x, -0.00046755850510571406, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(2)]).y, 0.00030490494761200812, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(2)]).z, -0.00036886672854369307, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(3)]).x, 0.00029624371924650601, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(3)]).y, 0.00013444493932002199, doublethreshold);
    EXPECT_NEAR(Setcell::fractional_displacement(mdcell->owned_atoms()[static_cast<std::size_t>(3)]).z, 9.0496812405138627e-05, doublethreshold);

    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(0)].vel.x, -2.3049731761587064e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(0)].vel.y, 7.1603385162874621e-06, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(0)].vel.z, 0.00016262437779022168, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(1)].vel.x, 0.0001961773016510733, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(1)].vel.y, 5.8637246942200678e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(1)].vel.z, 4.259822700946159e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(2)].vel.x, -0.00015692223255483009, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(2)].vel.y, 6.7034146380577021e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(2)].vel.z, -0.00017994277784966602, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(3)].vel.x, -3.5963807276704002e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(3)].vel.y, -8.5508938351509974e-05, doublethreshold);
    EXPECT_NEAR(mdcell->owned_atoms()[static_cast<std::size_t>(3)].vel.z, 9.6048301397465443e-05, doublethreshold);
}

TEST_F(Langevin_test, write_restart)
{
    mdrun->step_ = 1;
    mdrun->step_rst_ = 2;
    mdrun->write_restart(PARAM.sys.global_out_dir);

    std::ifstream ifs("Restart_md.txt");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("3"));
    ifs.close();
}

TEST_F(Langevin_test, restart)
{
    mdrun->restart(PARAM.sys.global_readin_dir);
    remove("Restart_md.txt");

    EXPECT_EQ(mdrun->step_rst_, 3);
}

TEST_F(Langevin_test, print_md)
{
    std::ofstream ofs("running_langevin.log");
    mdrun->print_md(ofs, true);
    ofs.close();

    std::ifstream ifs("running_langevin.log");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr(" ELECTRONIC      PART OF STRESS: 0.24609992 kbar"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr(" IONIC (KINETIC) PART OF STRESS: 0.83853919 kbar"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr(" MD PRESSURE (ELECTRONS+IONS)  : 1.0846391 kbar"));
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
    remove("running_langevin.log");
}
