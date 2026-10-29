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
