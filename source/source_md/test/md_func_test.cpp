#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_io/module_parameter/parameter.h"
#include "source_io/module_parameter/test_parameters.h"
#define private public
#define protected public
#include "md_test_fixture.h"
#include "source_md/md_func.h"
#undef protected
#undef private

#define doublethreshold 1e-12
/************************************************
 *  unit test of functions in md_func.h
 ***********************************************/

/**
 * - Tested Function
 *   - MD_func::gaussrand
 *     - genarate Gaussian random number
 *
 *   - MD_func::rand_vel
 *     - initialize atomic velocity randomly
 *
 *   - MD_func::get_mass_mbl
 *     - initialize atomic mass and degree of freedom
 *
 *   - MD_func::read_vel
 *     - read atomic velocity from STRU
 *
 *   - MD_func::rescale_vel
 *     - rescale the velocity to the target temperature
 *
 *   - MD_func::compute_stress
 *     - calculate the contribution of classical kinetic energy of atoms to stress
 *
 *   - MD_func::dump_info
 *     - output MD dump information
 *
 *   - MD_func::print_stress
 *     - output stress
 *
 *   - MD_func::current_md_info
 *     - test the current_md_info function with the correct file path
 *
 *   - MD_func::current_step_warning
 *     - test the current_md_info function with an incorrect file path
 */

class MD_func_test : public MdFuncTestFixture
{
};

TEST_F(MD_func_test, gaussrand)
{
    EXPECT_DOUBLE_EQ(MD_func::gaussrand(), 1.1122716058967226);
    EXPECT_DOUBLE_EQ(MD_func::gaussrand(), -0.34532367182326629);
    EXPECT_DOUBLE_EQ(MD_func::gaussrand(), 0.60805637857480721);
}

TEST_F(MD_func_test, readvel)
{
    MD_func::read_vel(ucell, vel);

    EXPECT_DOUBLE_EQ(vel[0].x, -0.0001320807363640);
    EXPECT_DOUBLE_EQ(vel[0].y, 7.13429429835e-05);
    EXPECT_DOUBLE_EQ(vel[0].z, -1.40179977966e-05);
    EXPECT_DOUBLE_EQ(vel[1].x, 0.000153039878532);
    EXPECT_DOUBLE_EQ(vel[1].y, -0.000146533266608);
    EXPECT_DOUBLE_EQ(vel[1].z, 9.64491480698e-05);
    EXPECT_DOUBLE_EQ(vel[2].x, -0.000133789480226);
    EXPECT_DOUBLE_EQ(vel[2].y, -3.0451038112e-06);
    EXPECT_DOUBLE_EQ(vel[2].z, -5.40998380137e-05);
    EXPECT_DOUBLE_EQ(vel[3].x, 0.000112830338059);
    EXPECT_DOUBLE_EQ(vel[3].y, 7.82354274358e-05);
    EXPECT_DOUBLE_EQ(vel[3].z, -2.83313122596e-05);
}

TEST_F(MD_func_test, RescaleVel)
{
    int frozen_freedom = 3;
    for (int i = 0; i < natom; ++i)
    {
        allmass[i] = 39.948 / ModuleBase::AU_to_MASS;
        vel[i].x = 0.1;
        vel[i].y = 0.2;
        vel[i].z = 0.3;
    }
    temperature = 300.0 / ModuleBase::Hartree_to_K;
    MD_func::rescale_vel(natom, temperature, allmass, frozen_freedom, vel);

    EXPECT_DOUBLE_EQ(vel[0].x, 4.579010775069125e-05);
    EXPECT_DOUBLE_EQ(vel[0].y, 9.15802155013825e-05);
    EXPECT_DOUBLE_EQ(vel[0].z, 0.00013737032325207373);
    EXPECT_DOUBLE_EQ(vel[1].x, 4.579010775069125e-05);
    EXPECT_DOUBLE_EQ(vel[1].y, 9.15802155013825e-05);
    EXPECT_DOUBLE_EQ(vel[1].z, 0.00013737032325207373);
    EXPECT_DOUBLE_EQ(vel[2].x, 4.579010775069125e-05);
    EXPECT_DOUBLE_EQ(vel[2].y, 9.15802155013825e-05);
    EXPECT_DOUBLE_EQ(vel[2].z, 0.00013737032325207373);
    EXPECT_DOUBLE_EQ(vel[3].x, 4.579010775069125e-05);
    EXPECT_DOUBLE_EQ(vel[3].y, 9.15802155013825e-05);
    EXPECT_DOUBLE_EQ(vel[3].z, 0.00013737032325207373);
}

TEST_F(MD_func_test, compute_stress)
{
    const ModuleBase::Vector3<double> test_velocity(0.1, 0.2, 0.3);
    MDCell mdcell;
    mdcell.initialize_from_unitcell(ucell, 0.0, 0.0, ModuleBase::world_communication_domain());
    for (LocalAtom& atom : mdcell.mutable_owned_atoms())
    {
        atom.vel = test_velocity;
    }
    MD_func::compute_stress(mdcell, true, virial, stress);
    EXPECT_DOUBLE_EQ(stress(0, 0), 2.9128300667662788);
    EXPECT_DOUBLE_EQ(stress(0, 1), 5.8256601335325575);
    EXPECT_DOUBLE_EQ(stress(0, 2), 8.7384902002988341);
    EXPECT_DOUBLE_EQ(stress(1, 0), 5.8256601335325575);
    EXPECT_DOUBLE_EQ(stress(1, 1), 11.651320267065115);
    EXPECT_DOUBLE_EQ(stress(1, 2), 17.476980400597668);
    EXPECT_DOUBLE_EQ(stress(2, 0), 8.7384902002988358);
    EXPECT_DOUBLE_EQ(stress(2, 1), 17.476980400597672);
    EXPECT_DOUBLE_EQ(stress(2, 2), 26.215470600896506);
}

TEST_F(MD_func_test, dump_info)
{
    MDCell mdcell;
    mdcell.initialize_from_unitcell(ucell, 0.0, 0.0, ModuleBase::world_communication_domain());
    for (LocalAtom& atom : mdcell.mutable_owned_atoms())
    {
        atom.vel = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    }
    MD_func::dump_info(0, TestParameters::sys().global_out_dir, mdcell, param_in, virial);
    std::ifstream ifs("MD_dump");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("MDSTEP:  0"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_CONSTANT: 0.529177000000 Angstrom"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_VECTORS"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  10.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  10.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  10.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("VIRIAL (kbar)"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr("INDEX    LABEL    POSITION (Angstrom)    FORCE (eV/Angstrom)    VELOCITY (Angstrom/fs)"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  0  Ar  0.000000000000  0.000000000000  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  1  Ar  2.751720222021  2.751720222021  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  2  Ar  2.698802525444  0.000000000000  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  3  Ar  0.000000000000  2.804637918599  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    ifs.close();

    // append
    MD_func::dump_info(1, TestParameters::sys().global_out_dir, mdcell, param_in, virial);
    std::ifstream ifs2("MD_dump");
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("MDSTEP:  0"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_CONSTANT: 0.529177000000 Angstrom"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_VECTORS"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  10.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  10.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  10.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("VIRIAL (kbar)"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr("INDEX    LABEL    POSITION (Angstrom)    FORCE (eV/Angstrom)    VELOCITY (Angstrom/fs)"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  0  Ar  0.000000000000  0.000000000000  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  1  Ar  2.751720222021  2.751720222021  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  2  Ar  2.698802525444  0.000000000000  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  3  Ar  0.000000000000  2.804637918599  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    getline(ifs2, output_str);
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("MDSTEP:  1"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_CONSTANT: 0.529177000000 Angstrom"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("LATTICE_VECTORS"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  10.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  10.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  10.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("VIRIAL (kbar)"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(
        output_str,
        testing::HasSubstr("INDEX    LABEL    POSITION (Angstrom)    FORCE (eV/Angstrom)    VELOCITY (Angstrom/fs)"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  0  Ar  0.000000000000  0.000000000000  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  1  Ar  2.751720222021  2.751720222021  0.000000000000  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  2  Ar  2.698802525444  0.000000000000  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    getline(ifs2, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr("  3  Ar  0.000000000000  2.804637918599  2.645884828867  0.000000000000  "
                                   "0.000000000000  0.000000000000  0.000000000000  0.000000000000  0.000000000000"));
    ifs2.close();

    remove("MD_dump");
}

TEST_F(MD_func_test, print_stress)
{
    GlobalV::ofs_running.open("running.log");
    MD_func::print_stress(GlobalV::ofs_running, virial, stress);

    std::ifstream ifs("running.log");
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("ELECTRONIC      PART OF STRESS: 0 kbar"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("IONIC (KINETIC) PART OF STRESS: 0 kbar"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("MD PRESSURE (ELECTRONS+IONS)  : 0 kbar"));
/*
    getline(ifs, output_str);
    getline(ifs, output_str);
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr(" ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><"));
    getline(ifs, output_str);
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr(" MD STRESS (kbar)"));
    getline(ifs, output_str);
    getline(ifs, output_str);
    EXPECT_THAT(output_str,
                testing::HasSubstr(" ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><"));
    getline(ifs, output_str);
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("              0              0              0"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("              0              0              0"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("              0              0              0"));
*/

    ifs.close();
    remove("running.log");
}

TEST_F(MD_func_test, current_md_info_mdcell_accepts_step_only_restart)
{
    std::ofstream file("Restart_md.txt");
    file << 123;
    file.close();

    MDCell mdcell;
    mdcell.initialize_from_unitcell(ucell, 0.0, 0.0, ModuleBase::world_communication_domain());
    int istep = -1;
    double temperature = 0.0;
    MD_func::current_md_info(mdcell, "./", istep, temperature);

    EXPECT_EQ(istep, 123);
    EXPECT_DOUBLE_EQ(temperature, 0.0);
    remove("Restart_md.txt");
}

TEST_F(MD_func_test, global_dof_mdcell)
{
    MDCell mdcell;
    mdcell.initialize_from_unitcell(ucell, 0.0, 0.0, ModuleBase::world_communication_domain());
    EXPECT_EQ(MD_func::global_dof(mdcell), 9);

    for (LocalAtom& atom : mdcell.mutable_owned_atoms())
    {
        atom.mbl.x = 0;
    }
    EXPECT_EQ(MD_func::global_dof(mdcell), 6);
}

TEST_F(MD_func_test, current_step_warning)
{
    MDCell mdcell;
    mdcell.initialize_from_unitcell(ucell, 0.0, 0.0, ModuleBase::world_communication_domain());
    int istep = 0;
    double temperature = 0.0;
    EXPECT_EXIT(MD_func::current_md_info(mdcell, "./", istep, temperature), ::testing::ExitedWithCode(1), "");
}
