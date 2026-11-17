#ifdef __MPI
#include "source_md/md_func.h"
#include "source_base/parallel_cell.h"
#include "source_base/parallel_common.h"
#include "source_base/parallel_global.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace
{
std::string make_dump_directory(const ModuleBase::CommunicationDomain& domain)
{
    std::string directory;
    if (domain.rank() == 0)
    {
        std::string pattern = testing::TempDir() + "abacus-md-dump-XXXXXX";
        char* created = mkdtemp(&pattern[0]);
        if (created != nullptr)
        {
            directory = created;
            directory += "/";
        }
    }
    Parallel_Common::bcast_string(directory);
    return directory;
}

class MdDumpBatchTest : public testing::TestWithParam<int>
{
};

TEST_P(MdDumpBatchTest, OrdersAtomsAcrossRanksAndAppendsCompleteFrames)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    ModuleBase::Matrix3 lattice;
    lattice.e11 = 1.0;
    lattice.e22 = 1.0;
    lattice.e33 = 1.0;
    const ModuleBase::Matrix3 inverse = lattice.Inverse();
    const int owner_count = std::min(domain.size(), 3);
    std::vector<LocalAtom> atoms;
    const int total_atoms = 13;
    for (int index = total_atoms - 1; index >= 0; --index)
    {
        const int owner = (2 * index + 1) % owner_count;
        if (owner != domain.rank())
        {
            continue;
        }
        LocalAtom atom;
        atom.type = index < 6 ? 0 : 1;
        atom.type_index = index - 6 * atom.type;
        const double x = 0.01 * (index + 1);
        atom.frac.set(x, 0.5, 0.5);
        atom.cart = atom.frac;
        atoms.push_back(atom);
    }
    const std::vector<std::string> labels = {"A", "B"};
    const std::vector<double> masses(2, 1.0);
    const std::vector<std::int64_t> counts = {6, 7};
    MDCell cell;
    cell.initialize_from_owned_atoms(lattice, inverse, 1.0, 1.0, total_atoms, atoms, labels, masses, counts, 0.0, domain);
    const std::string directory = make_dump_directory(domain);
    ASSERT_FALSE(directory.empty());
    MD_para parameters;
    ModuleBase::matrix virial(3, 3);
    MD_func::dump_info(0, directory, cell, parameters, false, virial);
    const int batch_bytes = GetParam();
    MD_func::dump_info(1, directory, cell, parameters, false, virial, batch_bytes);
    if (domain.rank() == 0)
    {
        const std::string filename = directory + "MD_dump";
        std::ifstream input(filename.c_str());
        ASSERT_TRUE(input);
        std::string line;
        int frame_count = 0;
        int atom_count = 0;
        std::size_t frame_bytes = 0;
        std::vector<std::size_t> rank_bytes(domain.size(), 0);
        std::vector<std::string> first_frame;
        bool reading_atoms = false;
        while (std::getline(input, line))
        {
            if (line.find("MDSTEP:") == 0)
            {
                if (frame_count > 0)
                {
                    EXPECT_EQ(atom_count, total_atoms);
                }
                ++frame_count;
                atom_count = 0;
                frame_bytes = 0;
                std::fill(rank_bytes.begin(), rank_bytes.end(), 0);
                reading_atoms = false;
            }
            else if (line.find("INDEX") == 0)
            {
                reading_atoms = true;
            }
            else if (reading_atoms && line.empty())
            {
                EXPECT_EQ(atom_count, total_atoms);
                if (batch_bytes == 512)
                {
                    // The whole frame exceeds the injected collective count limit.
                    EXPECT_GT(frame_bytes, 512u);
                }
                if (batch_bytes == 1024)
                {
                    EXPECT_GT(frame_bytes, 1024u);
                    if (domain.size() > 2)
                    {
                        // Reproduce the aggregate-only overflow: each rank fits on its own.
                        for (const std::size_t bytes : rank_bytes)
                        {
                            EXPECT_LE(bytes, 1024u);
                        }
                    }
                }
                reading_atoms = false;
            }
            else if (reading_atoms)
            {
                std::istringstream row(line);
                int index = -1;
                std::string label;
                double x = 0.0;
                row >> index >> label >> x;
                ASSERT_TRUE(row);
                EXPECT_EQ(index, atom_count);
                const std::string expected_label = atom_count < 6 ? "A" : "B";
                EXPECT_EQ(label, expected_label);
                const double expected_x = 0.01 * (atom_count + 1) / ModuleBase::ANGSTROM_AU;
                EXPECT_NEAR(x, expected_x, 1.0e-12);
                if (frame_count == 1)
                {
                    first_frame.push_back(line);
                }
                else
                {
                    ASSERT_LT(atom_count, first_frame.size());
                    EXPECT_EQ(line, first_frame[atom_count]);
                }
                frame_bytes += line.size() + 1;
                const int owner = (2 * index + 1) % owner_count;
                rank_bytes[owner] += line.size() + 1;
                ++atom_count;
            }
        }
        EXPECT_EQ(frame_count, 2);
        EXPECT_EQ(atom_count, total_atoms);
        input.close();
        std::remove(filename.c_str());
        rmdir(directory.c_str());
    }
}

// One-row, several-row, and normal budgets must produce identical complete frames.
INSTANTIATE_TEST_SUITE_P(BufferBudgets, MdDumpBatchTest, testing::Values(1, 512, 1024, 1048576));

TEST(MdDumpMpiTest, BoundedGatherRejectsOversizedBatchesCollectively)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    const int local_bytes = domain.rank() == 0 ? 17 : 0;
    const std::string oversized(local_bytes, 'x');
    EXPECT_THROW(domain.gather_text(oversized, 16), std::length_error);
    if (domain.size() > 2)
    {
        // Every rank fits individually; their aggregate must also be checked on all ranks.
        const std::string local(8, 'x');
        EXPECT_THROW(domain.gather_text(local, 16), std::length_error);
    }
    const std::string empty;
    const std::string gathered = domain.gather_text(empty, 1);
    EXPECT_TRUE(gathered.empty());
}

TEST(MdDumpMpiTest, EmptyRanksCanWriteAnEmptyFrame)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    ModuleBase::Matrix3 lattice;
    lattice.e11 = 1.0;
    lattice.e22 = 1.0;
    lattice.e33 = 1.0;
    const ModuleBase::Matrix3 inverse = lattice.Inverse();
    const std::vector<LocalAtom> atoms;
    const std::vector<std::string> labels(1, "X");
    const std::vector<double> masses(1, 1.0);
    const std::vector<std::int64_t> counts(1, 0);
    MDCell cell;
    cell.initialize_from_owned_atoms(lattice, inverse, 1.0, 1.0, 0, atoms, labels, masses, counts, 0.0, domain);
    const std::string directory = make_dump_directory(domain);
    ASSERT_FALSE(directory.empty());
    MD_para parameters;
    ModuleBase::matrix virial(3, 3);
    MD_func::dump_info(0, directory, cell, parameters, false, virial);
    if (domain.rank() == 0)
    {
        const std::string filename = directory + "MD_dump";
        std::ifstream input(filename.c_str());
        ASSERT_TRUE(input);
        const std::string contents((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
        EXPECT_NE(contents.find("MDSTEP:  0"), std::string::npos);
        EXPECT_NE(contents.find("INDEX"), std::string::npos);
        input.close();
        std::remove(filename.c_str());
        rmdir(directory.c_str());
    }
}
} // namespace

int main(int argc, char** argv)
{
    int process_count = 1;
    int thread_count = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, rank);
    // These application communicators are not created by this focused test.
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    ::testing::InitGoogleTest(&argc, argv);
    const int local_result = RUN_ALL_TESTS();
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    const int result = domain.max(local_result);
    Parallel_Global::finalize_mpi();
    return result;
}

#endif

#ifndef __MPI
#include "source_cell/dom_decomp.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "source_md/unittests/md_test_fixture.h"
#include "source_md/md_func.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <unistd.h>

constexpr double doublethreshold = 1e-12;

namespace
{
void initialize_mdcell_from_ucell(MDCell& mdcell, UnitCell& ucell, DomainDecomposition& decomp)
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
    mdcell.initialize_from_owned_atoms(ucell.latvec,
                                       ucell.GT,
                                       ucell.lat0,
                                       ucell.omega,
                                       ucell.nat,
                                       owned_atoms,
                                       type_labels,
                                       type_masses,
                                       type_atom_counts,
                                       0.0,
                                       comm_domain);
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
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    for (LocalAtom& atom : mdcell.owned_atoms())
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
    std::string pattern = testing::TempDir() + "abacus-md-func-XXXXXX";
    char* directory = mkdtemp(&pattern[0]);
    ASSERT_NE(directory, nullptr);
    const std::string out_dir = pattern + "/";
    const std::string dump_path = out_dir + "MD_dump";
    MDCell mdcell;
    DomainDecomposition decomp;
    initialize_mdcell_from_ucell(mdcell, ucell, decomp);
    for (LocalAtom& atom : mdcell.owned_atoms())
    {
        atom.vel = ModuleBase::Vector3<double>(0.0, 0.0, 0.0);
    }
    MD_func::dump_info(0, out_dir, mdcell, inp.mdp, inp.cal_stress, virial);
    std::ifstream ifs(dump_path.c_str());
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
    MD_func::dump_info(1, out_dir, mdcell, inp.mdp, inp.cal_stress, virial);
    std::ifstream ifs2(dump_path.c_str());
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

    std::remove(dump_path.c_str());
    rmdir(pattern.c_str());
}

TEST_F(MD_func_test, print_stress)
{
    const std::string log_path = testing::TempDir() + "abacus-md-print-stress.log";
    std::ofstream running_log(log_path.c_str());
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

    for (LocalAtom& atom : mdcell.owned_atoms())
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
    char* created = mkdtemp(&pattern[0]);
    ASSERT_NE(created, nullptr);
    const std::string empty_directory = pattern + "/";
    EXPECT_EXIT(MD_func::current_md_info(mdcell, empty_directory, istep, temperature),
                ::testing::ExitedWithCode(1),
                "");
    EXPECT_EQ(rmdir(pattern.c_str()), 0);
}

#include "source_md/md_func.h"
#include "source_cell/mdcell.h"
#include "source_cell/dom_decomp.h"
#include "source_base/parallel_cell.h"
#include "source_esolver/esolver.h"

#include <gtest/gtest.h>

namespace
{
class GhostForceSolver : public ModuleESolver::ESolver
{
public:
    void before_all_runners(BaseCell&, const Input_para&) override {}
    void after_all_runners(BaseCell&) override {}
    void runner(BaseCell& base, const int) override
    {
        MDCell& cell = static_cast<MDCell&>(base);
        EXPECT_TRUE(cell.has_neighbor_search());
        EXPECT_GT(cell.ghost_atoms().size(), 0);
        for (LocalAtom& atom : cell.owned_atoms())
        {
            atom.force.set(4.0, 0.0, 0.0);
        }
        for (LocalAtom& atom : cell.ghost_atoms())
        {
            EXPECT_DOUBLE_EQ(atom.force.norm2(), 0.0);
            atom.force.set(2.0, 0.0, 0.0);
        }
    }
    double cal_energy() override { return 6.0; }
    void cal_force(BaseCell&, ModuleBase::matrix&) override
    {
        ADD_FAILURE() << "Direct MDCell forces must not use cal_force.";
    }
    void cal_stress(BaseCell&, ModuleBase::matrix& stress) override
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
    cell.initialize_from_owned_atoms(lattice, lattice.Inverse(), 1.0, 64.0, 1,
                                     std::vector<LocalAtom>(1, atom),
                                     std::vector<std::string>(1, "X"),
                                     std::vector<double>(1, 1.0),
                                     std::vector<std::int64_t>(1, 1), 0.4, domain);
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
        for (const LocalAtom& ghost : cell.ghost_atoms())
        {
            EXPECT_DOUBLE_EQ(ghost.force.norm2(), 0.0);
        }
    }
}
} // namespace

#endif
