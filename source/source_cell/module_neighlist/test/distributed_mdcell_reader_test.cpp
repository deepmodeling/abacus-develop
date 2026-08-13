#include <gtest/gtest.h>

#include "source_cell/distributed_mdcell_reader.h"
#include "source_cell/md_cell.h"
#include "source_cell/print_cell.h"
#include "source_base/constants.h"
#include "source_base/parallel_cell.h"
#include "source_base/global_variable.h"
#include "source_cell/module_neighlist/domain_decomposition.h"

#include <cstdio>
#include <cstdint>
#include <fstream>
#include <mpi.h>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>

static_assert(!std::is_copy_constructible<MDCell>::value, "MDCell must not be copy constructible.");
static_assert(!std::is_copy_assignable<MDCell>::value, "MDCell must not be copy assignable.");
static_assert(std::is_move_constructible<MDCell>::value, "MDCell must be move constructible.");

namespace
{
void write_cartesian_stru_case(const std::string& stru_file)
{
    std::ofstream ofs(stru_file.c_str());
    ofs << "ATOMIC_SPECIES\n";
    ofs << "He 4.0026 auto auto\n\n";
    ofs << "LATTICE_CONSTANT\n";
    ofs << "1.0\n\n";
    ofs << "LATTICE_VECTORS\n";
    ofs << "4.0 0.0 0.0\n";
    ofs << "0.0 4.0 0.0\n";
    ofs << "0.0 0.0 4.0\n\n";
    ofs << "ATOMIC_POSITIONS\n";
    ofs << "Cartesian\n\n";
    ofs << "He\n";
    ofs << "0.0\n";
    ofs << "4\n";
    ofs << "0.40 0.40 0.40 m 1 1 1 v 0.01 0.00 0.00\n";
    ofs << "2.40 0.40 0.40 m 1 0 1 v 0.02 0.00 0.00\n";
    ofs << "0.40 2.40 0.40 m 0 1 1 v 0.03 0.00 0.00\n";
    ofs << "2.40 2.40 0.40 m 1 1 0 v 0.04 0.00 0.00\n";
}

ModuleBase::Matrix3 make_lattice()
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 4.0;
    latvec.e12 = 0.0;
    latvec.e13 = 0.0;
    latvec.e21 = 0.0;
    latvec.e22 = 4.0;
    latvec.e23 = 0.0;
    latvec.e31 = 0.0;
    latvec.e32 = 0.0;
    latvec.e33 = 4.0;
    return latvec;
}
} // namespace

TEST(DistributedMDCellReaderTest, ReadOwnedAtomsFromSTRUWithoutUnitCell)
{
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const std::string stru_file = "distributed_mdcell_reader_cartesian.STRU";
    if (world_rank == 0)
    {
        write_cartesian_stru_case(stru_file);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm md_comm = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, world_rank % 2, world_rank, &md_comm);
    const ModuleBase::CommunicationDomain communication_domain(md_comm);

    MdStruFileMetadata stru_metadata;
    MDCell mdcell = DistributedMDCellReader::read_stru(stru_file,
                                                        std::vector<int>{1, 1, 1},
                                                        1.0 * ModuleBase::ANGSTROM_AU,
                                                        0.0,
                                                        stru_metadata,
                                                        communication_domain);

    EXPECT_EQ(mdcell.type_labels().size(), 1U);
    EXPECT_EQ(mdcell.type_labels()[0], "He");
    ASSERT_EQ(mdcell.type_masses().size(), 1U);
    EXPECT_DOUBLE_EQ(mdcell.type_masses()[0], 4.0026);
    ASSERT_EQ(mdcell.type_atom_counts().size(), 1U);
    EXPECT_EQ(mdcell.type_atom_counts()[0], 4);
    ASSERT_EQ(stru_metadata.species.size(), 1U);
    EXPECT_EQ(stru_metadata.species[0].pseudo_file, "auto");
    EXPECT_EQ(stru_metadata.species[0].pseudo_type, "auto");
    EXPECT_EQ(mdcell.nat(), 4);

    DomainDecomposition decomp;
    decomp.init(md_comm, make_lattice(), 1.0, 1.0 * ModuleBase::ANGSTROM_AU, 0.0);

    long long local_count = static_cast<long long>(mdcell.owned_atoms().size());
    long long global_count = 0;
    MPI_Allreduce(&local_count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, md_comm);
    EXPECT_EQ(global_count, 4);

    std::set<std::pair<int, int> > local_ids;
    for (std::size_t iat = 0; iat < mdcell.owned_atoms().size(); ++iat)
    {
        const LocalAtom& atom = mdcell.owned_atoms()[iat];
        EXPECT_EQ(decomp.owner_rank_from_frac(atom.frac), communication_domain.rank());
        local_ids.insert(std::make_pair(atom.type, atom.type_index));
        EXPECT_GE(atom.type, 0);
        EXPECT_DOUBLE_EQ(atom.force.x, 0.0);
        EXPECT_DOUBLE_EQ(atom.force.y, 0.0);
        EXPECT_DOUBLE_EQ(atom.force.z, 0.0);
    }
    EXPECT_EQ(local_ids.size(), mdcell.owned_atoms().size());

    bool saw_v01 = false;
    bool saw_v04 = false;
    for (std::size_t iat = 0; iat < mdcell.owned_atoms().size(); ++iat)
    {
        const LocalAtom& atom = mdcell.owned_atoms()[iat];
        if (atom.type == 0 && atom.type_index == 0)
        {
            saw_v01 = true;
            EXPECT_DOUBLE_EQ(atom.vel.x, 0.01);
            EXPECT_EQ(atom.mbl.x, 1);
            EXPECT_EQ(atom.mbl.y, 1);
            EXPECT_EQ(atom.mbl.z, 1);
            EXPECT_DOUBLE_EQ(atom.mass, 4.0026 / ModuleBase::AU_to_MASS);
        }
        if (atom.type == 0 && atom.type_index == 3)
        {
            saw_v04 = true;
            EXPECT_DOUBLE_EQ(atom.vel.x, 0.04);
            EXPECT_EQ(atom.mbl.x, 1);
            EXPECT_EQ(atom.mbl.y, 1);
            EXPECT_EQ(atom.mbl.z, 0);
        }
    }

    const int saw_flags[2] = {saw_v01 ? 1 : 0, saw_v04 ? 1 : 0};
    int reduced_flags[2] = {0, 0};
    MPI_Allreduce(saw_flags, reduced_flags, 2, MPI_INT, MPI_MAX, md_comm);
    EXPECT_EQ(reduced_flags[0], 1);
    EXPECT_EQ(reduced_flags[1], 1);

    MPI_Comm_free(&md_comm);
}

TEST(DistributedMDCellReaderTest, RestartStruPreservesAtomRecordsAcrossRanks)
{
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ASSERT_GE(size, 4);

    std::vector<LocalAtom> owned_atoms;
    if (rank == 0)
    {
        owned_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(11.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.11, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<int>(1, 0, 1),
                                        1.0,
                                        1,
                                        1,
                                        rank,
                                        false));
    }
    if (rank == 1)
    {
        owned_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(1.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.01, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<int>(0, 1, 1),
                                        1.0,
                                        0,
                                        1,
                                        rank,
                                        false));
    }
    if (rank == 2)
    {
        owned_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(10.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.10, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<int>(1, 1, 0),
                                        1.0,
                                        1,
                                        0,
                                        rank,
                                        false));
    }
    if (rank == 3)
    {
        owned_atoms.push_back(LocalAtom(ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<int>(0, 0, 1),
                                        1.0,
                                        0,
                                        0,
                                        rank,
                                        false));
    }

    ModuleBase::Matrix3 lattice = make_lattice();
    lattice.e11 = 20.0;
    lattice.e22 = 20.0;
    lattice.e33 = 20.0;
    MDCell mdcell(lattice,
                  lattice.Inverse(),
                  1.0,
                  1.0,
                  4,
                  owned_atoms,
                  std::vector<std::string>{"A", "B"},
                  std::vector<double>{1.0, 1.0},
                  std::vector<std::int64_t>{2, 2},
                  0.0,
                  0.0,
                  ModuleBase::world_communication_domain());
    MdStruFileMetadata metadata;
    metadata.species.resize(2);
    const std::string output_file = "distributed_mdcell_restart.STRU";
    mdcell::print_stru_file(mdcell, metadata, output_file);

    MdStruFileMetadata round_trip_metadata;
    MDCell round_trip = DistributedMDCellReader::read_stru(output_file,
                                                            std::vector<int>{1, 1, 1},
                                                            0.1,
                                                            0.0,
                                                            round_trip_metadata,
                                                            ModuleBase::world_communication_domain());
    double local_positions[4] = {0.0, 0.0, 0.0, 0.0};
    double local_velocities[4] = {0.0, 0.0, 0.0, 0.0};
    int local_mbl_x[4] = {0, 0, 0, 0};
    int local_owners[4] = {0, 0, 0, 0};
    for (std::size_t iat = 0; iat < round_trip.owned_atoms().size(); ++iat)
    {
        const LocalAtom& atom = round_trip.owned_atoms()[iat];
        const int index = atom.type == 0 ? static_cast<int>(atom.cart.x)
                                         : 2 + static_cast<int>(atom.cart.x) - 10;
        local_positions[index] = atom.cart.x;
        local_velocities[index] = atom.vel.x;
        local_mbl_x[index] = atom.mbl.x;
        local_owners[index] = 1;
    }
    double global_positions[4] = {0.0, 0.0, 0.0, 0.0};
    double global_velocities[4] = {0.0, 0.0, 0.0, 0.0};
    int global_mbl_x[4] = {0, 0, 0, 0};
    int global_owners[4] = {0, 0, 0, 0};
    MPI_Allreduce(local_positions, global_positions, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_velocities, global_velocities, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_mbl_x, global_mbl_x, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(local_owners, global_owners, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    for (int iat = 0; iat < 4; ++iat) EXPECT_EQ(global_owners[iat], 1);
    EXPECT_DOUBLE_EQ(global_positions[0], 0.0);
    EXPECT_DOUBLE_EQ(global_positions[1], 1.0);
    EXPECT_DOUBLE_EQ(global_positions[2], 10.0);
    EXPECT_DOUBLE_EQ(global_positions[3], 11.0);
    EXPECT_DOUBLE_EQ(global_velocities[0], 0.0);
    EXPECT_DOUBLE_EQ(global_velocities[1], 0.01);
    EXPECT_DOUBLE_EQ(global_velocities[2], 0.10);
    EXPECT_DOUBLE_EQ(global_velocities[3], 0.11);
    EXPECT_EQ(global_mbl_x[0], 0);
    EXPECT_EQ(global_mbl_x[1], 0);
    EXPECT_EQ(global_mbl_x[2], 1);
    EXPECT_EQ(global_mbl_x[3], 1);

    if (rank == 0)
    {
        std::ifstream input(output_file.c_str());
        ASSERT_TRUE(input.good());
        std::vector<double> coordinates_a;
        std::vector<double> coordinates_b;
        std::string line;
        int current_type = -1;
        int atoms_remaining = 0;
        while (std::getline(input, line))
        {
            if (line == "A #label")
            {
                current_type = 0;
                continue;
            }
            if (line == "B #label")
            {
                current_type = 1;
                continue;
            }
            if (current_type >= 0 && line.find("#number of atoms") != std::string::npos)
            {
                std::istringstream count_stream(line);
                count_stream >> atoms_remaining;
                continue;
            }
            if (current_type < 0 || atoms_remaining == 0) continue;
            std::istringstream values(line);
            double x = 0.0;
            values >> x;
            --atoms_remaining;
            if (values)
            {
                if (current_type == 0) coordinates_a.push_back(x);
                else coordinates_b.push_back(x);
            }
        }
        ASSERT_EQ(coordinates_a.size(), 2U);
        ASSERT_EQ(coordinates_b.size(), 2U);
        std::set<double> expected_a = {0.0, 1.0};
        std::set<double> expected_b = {10.0, 11.0};
        EXPECT_EQ(std::set<double>(coordinates_a.begin(), coordinates_a.end()), expected_a);
        EXPECT_EQ(std::set<double>(coordinates_b.begin(), coordinates_b.end()), expected_b);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) std::remove(output_file.c_str());
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &GlobalV::MY_RANK);
    MPI_Comm_size(MPI_COMM_WORLD, &GlobalV::NPROC);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
