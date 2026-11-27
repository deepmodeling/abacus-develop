
#ifndef __MPI
#include "source_cell/dom_decomp.h"
#include "source_md/md_func.h"
#include "source_md/unittests/md_test_fixture.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <unistd.h>

constexpr double doublethreshold = 1e-12;

namespace
{
void initialize_mdcell_from_ucell(MDCell &mdcell, UnitCell &ucell, DomainDecomposition &decomp)
{
    const ModuleBase::CommunicationDomain comm_domain = ModuleBase::world_comm_domain();
    decomp.init(comm_domain, ucell.latvec, ucell.lat0, 0.0, 0.0);
    const std::vector<LocalAtom> owned_atoms = decomp.split_owned_atoms_from_ucell(ucell);
    std::vector<std::string> type_labels;
    std::vector<double> type_masses;
    std::vector<std::int64_t> type_atom_counts;
    for (int it = 0; it < ucell.ntype; ++it)
    {
        type_labels.push_back(ucell.atoms[it].label);
        type_masses.push_back(ucell.atoms[it].mass);
        type_atom_counts.push_back(ucell.atoms[it].na);
    }
    mdcell.initialize_from_owned_atoms(
        ucell.latvec, ucell.GT, ucell.lat0, ucell.omega, ucell.nat, owned_atoms, type_labels, type_masses,
        type_atom_counts, 0.0, comm_domain);
    mdcell.set_backing_unitcell(ucell);
}
} // namespace

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
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    for (LocalAtom &atom : mdcell.owned_atoms())
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

TEST_F(MD_func_test, print_stress)
{
    const std::string log_path = testing::TempDir() + "abacus-md-print-stress.log";
    std::ofstream running_log(log_path.c_str());
    virial(0, 1) = 1.0;
    stress(0, 1) = 2.0;
    MD_func::print_stress(running_log, virial, stress);
    running_log.close();

    std::ifstream ifs(log_path.c_str());
    std::string output_str;
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("ELECTRONIC      PART OF STRESS: 0 kbar"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("IONIC (KINETIC) PART OF STRESS: 0 kbar"));
    getline(ifs, output_str);
    EXPECT_THAT(output_str, testing::HasSubstr("MD PRESSURE (ELECTRONS+IONS)  : 0 kbar"));
    getline(ifs, output_str);
    EXPECT_EQ(output_str, " VIRIAL (kbar):");
    getline(ifs, output_str);
    std::istringstream first_virial_row(output_str);
    double vxx = 0.0;
    double vxy = 0.0;
    double vxz = 0.0;
    first_virial_row >> vxx >> vxy >> vxz;
    ASSERT_FALSE(first_virial_row.fail());
    const double kbar_per_atomic_unit = ModuleBase::HARTREE_SI
                                        / std::pow(ModuleBase::BOHR_RADIUS_SI, 3) * 1.0e-8;
    EXPECT_DOUBLE_EQ(vxx, 0.0);
    EXPECT_NEAR(vxy, kbar_per_atomic_unit, 1.0e-6);
    EXPECT_DOUBLE_EQ(vxz, 0.0);
    getline(ifs, output_str);
    getline(ifs, output_str);
    getline(ifs, output_str);
    EXPECT_EQ(output_str, " MD TOTAL STRESS (kbar):");
    getline(ifs, output_str);
    std::istringstream first_stress_row(output_str);
    double sxx = 0.0;
    double sxy = 0.0;
    double sxz = 0.0;
    first_stress_row >> sxx >> sxy >> sxz;
    ASSERT_FALSE(first_stress_row.fail());
    EXPECT_DOUBLE_EQ(sxx, 0.0);
    const double expected_sxy = 2.0 * kbar_per_atomic_unit;
    EXPECT_NEAR(sxy, expected_sxy, 1.0e-6);
    EXPECT_DOUBLE_EQ(sxz, 0.0);
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
    std::remove(log_path.c_str());
}

TEST_F(MD_func_test, init_vel_rejects_zero_stru_velocities_with_backing_unitcell)
{
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            ucell.atoms[it].vel[ia].set(0.0, 0.0, 0.0);
        }
    }
    MDCell mdcell;
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    double temperature = 300.0 / ModuleBase::Hartree_to_K;
    std::int64_t frozen_freedom = 0;
    EXPECT_EXIT(
        MD_func::init_vel(mdcell, true, false, temperature, frozen_freedom),
        ::testing::ExitedWithCode(1),
        "");
}

TEST_F(MD_func_test, init_vel_rejects_zero_stru_velocities_without_backing_unitcell)
{
    MDCell mdcell = Setcell::setup_mdcell(ucell);
    for (LocalAtom &atom : mdcell.owned_atoms())
    {
        atom.vel.set(0.0, 0.0, 0.0);
    }
    double temperature = 300.0 / ModuleBase::Hartree_to_K;
    std::int64_t frozen_freedom = 0;
    EXPECT_EXIT(
        MD_func::init_vel(mdcell, true, false, temperature, frozen_freedom),
        ::testing::ExitedWithCode(1),
        "");
}

TEST_F(MD_func_test, init_vel_allows_zero_stru_velocities_at_zero_temperature_with_backing_unitcell)
{
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            ucell.atoms[it].vel[ia].set(0.0, 0.0, 0.0);
        }
    }
    MDCell mdcell;
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    double temperature = 0.0;
    std::int64_t frozen_freedom = 0;
    MD_func::init_vel(mdcell, true, false, temperature, frozen_freedom);
    for (const LocalAtom &atom : mdcell.owned_atoms())
    {
        for (int k = 0; k < 3; ++k)
        {
            EXPECT_DOUBLE_EQ(atom.vel[k], 0.0);
        }
    }
}

TEST_F(MD_func_test, init_vel_clears_nonzero_stru_velocities_at_zero_temperature_with_backing_unitcell)
{
    const ModuleBase::Vector3<double> initial_velocity(0.1, -0.2, 0.3);
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            ucell.atoms[it].vel[ia] = initial_velocity;
        }
    }
    MDCell mdcell;
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    double temperature = 0.0;
    std::int64_t frozen_freedom = 0;
    MD_func::init_vel(mdcell, true, false, temperature, frozen_freedom);
    for (const LocalAtom &atom : mdcell.owned_atoms())
    {
        for (int k = 0; k < 3; ++k)
        {
            EXPECT_DOUBLE_EQ(atom.vel[k], 0.0);
        }
    }
}

TEST_F(MD_func_test, init_vel_restart_allows_zero_stru_velocities)
{
    for (int it = 0; it < ucell.ntype; ++it)
    {
        for (int ia = 0; ia < ucell.atoms[it].na; ++ia)
        {
            ucell.atoms[it].vel[ia].set(0.0, 0.0, 0.0);
        }
    }
    MDCell mdcell;
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    double temperature = 300.0 / ModuleBase::Hartree_to_K;
    std::int64_t frozen_freedom = 0;
    // restart = true must not abort even when the STRU velocities are all zero,
    // because the velocities will be restored from the restart file later.
    MD_func::init_vel(mdcell, true, true, temperature, frozen_freedom);
    SUCCEED();
}

TEST_F(MD_func_test, init_vel_scales_nonzero_stru_velocities_to_target_temperature)
{
    MDCell mdcell;
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    double temperature = 300.0 / ModuleBase::Hartree_to_K;
    std::int64_t frozen_freedom = 0;
    MD_func::init_vel(mdcell, true, false, temperature, frozen_freedom);
    double kinetic = 0.0;
    const double current = MD_func::current_temp(kinetic, mdcell, frozen_freedom);
    EXPECT_NEAR(current, temperature, 1e-12);
}

TEST_F(MD_func_test, current_md_info_mdcell_accepts_step_only_restart)
{
    const std::string restart_file = testing::TempDir() + "Restart_md.txt";
    std::ofstream file(restart_file.c_str());
    file << 123;
    file.close();

    MDCell mdcell;
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    int istep = -1;
    double temperature = 0.0;
    const std::string restart_directory = testing::TempDir();
    MD_func::current_md_info(mdcell, restart_directory, istep, temperature);

    EXPECT_EQ(istep, 123);
    EXPECT_DOUBLE_EQ(temperature, 0.0);
    std::remove(restart_file.c_str());
}

TEST_F(MD_func_test, global_dof_mdcell)
{
    MDCell mdcell;
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    EXPECT_EQ(MD_func::global_dof(mdcell), 9);

    for (LocalAtom &atom : mdcell.owned_atoms())
    {
        atom.mbl.x = 0;
    }
    EXPECT_EQ(MD_func::global_dof(mdcell), 6);
}

TEST_F(MD_func_test, current_step_warning)
{
    MDCell mdcell;
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    int istep = 0;
    double temperature = 0.0;
    std::string pattern = testing::TempDir() + "abacus-md-missing-restart-XXXXXX";
    char *created = mkdtemp(&pattern[0]);
    ASSERT_NE(created, nullptr);
    const std::string empty_directory = pattern + "/";
    EXPECT_EXIT(
        MD_func::current_md_info(mdcell, empty_directory, istep, temperature), ::testing::ExitedWithCode(1), "");
    EXPECT_EQ(rmdir(pattern.c_str()), 0);
}

#include "source_base/parallel_cell.h"
#include "source_cell/dom_decomp.h"
#include "source_cell/mdcell.h"
#include "source_esolver/esolver.h"
#include "source_md/md_func.h"

#include <gtest/gtest.h>

namespace
{
class GhostForceSolver : public ModuleESolver::ESolver
{
  public:
    void before_all_runners(BaseCell &, const Input_para &) override {}
    void after_all_runners(BaseCell &) override {}
    void runner(BaseCell &base, const int) override
    {
        MDCell &cell = static_cast<MDCell &>(base);
        EXPECT_TRUE(cell.has_neighbor_search());
        EXPECT_GT(cell.ghost_atoms().size(), 0);
        for (LocalAtom &atom : cell.owned_atoms())
        {
            atom.force.set(4.0, 0.0, 0.0);
        }
        for (LocalAtom &atom : cell.ghost_atoms())
        {
            EXPECT_DOUBLE_EQ(atom.force.norm2(), 0.0);
            atom.force.set(2.0, 0.0, 0.0);
        }
    }
    double cal_energy() override { return 6.0; }
    void cal_force(BaseCell &, ModuleBase::matrix &) override
    {
        ADD_FAILURE() << "Direct MDCell forces must not use cal_force.";
    }
    void cal_stress(BaseCell &, ModuleBase::matrix &stress) override
    {
        stress.zero_out();
        stress(0, 0) = 8.0;
    }
};

TEST(MdForceVirialTest, ReturnsGhostForcesBeforeConvertingUnitsExactlyOnce)
{
    MDCell cell;
    DomainDecomposition decomp;
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    ModuleBase::Matrix3 lattice;
    lattice.e11 = 4.0;
    lattice.e22 = 4.0;
    lattice.e33 = 4.0;
    LocalAtom atom;
    atom.frac.set(0.05, 0.5, 0.5);
    atom.cart = atom.frac * lattice;
    cell.initialize_from_owned_atoms(
        lattice, lattice.Inverse(), 1.0, 64.0, 1, std::vector<LocalAtom>(1, atom), std::vector<std::string>(1, "X"),
        std::vector<double>(1, 1.0), std::vector<std::int64_t>(1, 1), 0.4, domain);
    decomp.init(domain, lattice, 1.0, 0.0, 0.4);
    cell.set_neighbor_cutoff(0.5);
    GhostForceSolver solver;
    double potential = 0.0;
    ModuleBase::matrix virial(3, 3);
    for (int step = 0; step < 2; ++step)
    {
        const bool positions_updated = step > 0;
        MD_func::force_virial(&solver, step, positions_updated, cell, decomp, potential, true, virial);
        EXPECT_DOUBLE_EQ(potential, 3.0);
        EXPECT_DOUBLE_EQ(virial(0, 0), 4.0);
        ASSERT_EQ(cell.owned_atoms().size(), 1);
        EXPECT_DOUBLE_EQ(cell.owned_atoms()[0].force.x, 2.0 + cell.ghost_atoms().size());
        for (const LocalAtom &ghost : cell.ghost_atoms())
        {
            EXPECT_DOUBLE_EQ(ghost.force.norm2(), 0.0);
        }
    }
}
} // namespace

#endif
