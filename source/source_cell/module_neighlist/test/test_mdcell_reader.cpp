#include "source_cell/domain_decomposition.h"
#include <gtest/gtest.h>

#include "source_cell/mdcell_reader.h"
#include "source_cell/mdcell.h"
#include "source_cell/print_cell.h"
#include "source_base/constants.h"
#include "source_base/parallel_cell.h"
#include "source_base/global_variable.h"

#include <cstdio>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>

static_assert(!std::is_copy_constructible<MDCell>::value, "MDCell must not be copy constructible.");
static_assert(!std::is_copy_assignable<MDCell>::value, "MDCell must not be copy assignable.");
static_assert(std::is_default_constructible<MDCell>::value, "MDCell must be default constructible.");
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

TEST(MDCellReaderTest, ReadOwnedAtomsFromSTRUWithoutUnitCell)
{
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const std::string stru_file = "mdcell_reader_cartesian.STRU";
    if (world_rank == 0)
    {
        write_cartesian_stru_case(stru_file);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Comm md_comm = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, world_rank % 2, world_rank, &md_comm);
    ModuleBase::CommunicationDomain comm_domain;
    comm_domain.initialize(md_comm);

    DomainDecomposition reader_decomp;
    MDCell mdcell = MDCellReader::read_stru(stru_file,
                                             std::vector<int>{1, 1, 1},
                                             0.0,
                                             comm_domain, reader_decomp);

    EXPECT_EQ(mdcell.type_labels().size(), 1U);
    EXPECT_EQ(mdcell.type_labels()[0], "He");
    ASSERT_EQ(mdcell.type_masses().size(), 1U);
    EXPECT_DOUBLE_EQ(mdcell.type_masses()[0], 4.0026);
    ASSERT_EQ(mdcell.type_atom_counts().size(), 1U);
    EXPECT_EQ(mdcell.type_atom_counts()[0], 4);
    ASSERT_EQ(mdcell.stru_meta().species.size(), 1U);
    EXPECT_EQ(mdcell.stru_meta().species[0].pseudo_file, "auto");
    EXPECT_EQ(mdcell.stru_meta().species[0].pseudo_type, "auto");
    EXPECT_EQ(mdcell.nat(), 4);

    DomainDecomposition decomp;
    decomp.init(comm_domain, make_lattice(), 1.0, 1.0 * ModuleBase::ANGSTROM_AU, 0.0);

    long long local_count = static_cast<long long>(mdcell.owned_atoms().size());
    long long global_count = 0;
    MPI_Allreduce(&local_count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, md_comm);
    EXPECT_EQ(global_count, 4);

    std::set<std::pair<int, int> > local_ids;
    for (std::size_t iat = 0; iat < mdcell.owned_atoms().size(); ++iat)
    {
        const LocalAtom& atom = mdcell.owned_atoms()[iat];
        EXPECT_EQ(decomp.owner_rank_from_frac(atom.frac), comm_domain.rank());
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

// All inputs describe the same unwrapped fractional positions in a skew cell.
// Expected positions below are independent of the reader's inverse transform.
class MDCellReaderUnitsTest : public testing::TestWithParam<const char*>
{};

TEST_P(MDCellReaderUnitsTest, ConvertsWrapsAndReplicates)
{
    const ModuleBase::CommunicationDomain comm = ModuleBase::world_comm_domain();
    const double lat0 = 2.5;
    const double tolerance = 1.0e-12;
    const double input_frac[4][3] = {{0.1, 0.2, 0.3}, {-0.2, 1.3, 0.4},
                                     {1.4, -0.3, 2.2}, {-1.3, 2.4, -0.2}};
    const double wrapped[4][3] = {{0.1, 0.2, 0.3}, {0.8, 0.3, 0.4},
                                  {0.4, 0.7, 0.2}, {0.7, 0.4, 0.8}};
    const int moves[4][3] = {{1, 1, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 0}};
    const std::string units = GetParam();
    std::ostringstream filename;
    filename << testing::TempDir() << "mdcell_reader_" << units << "_rank" << comm.rank() << ".STRU";
    const std::string stru_file = filename.str();
    double scale = 1.0;
    if (units == "Cartesian_au") scale = lat0;
    if (units == "Cartesian_angstrom") scale = lat0 * ModuleBase::BOHR_TO_A;
    {
        std::ofstream ofs(stru_file.c_str());
        const int local_open_failure = ofs.good() ? 0 : 1;
        const int open_failures = comm.sum(local_open_failure);
        EXPECT_EQ(open_failures, 0) << "Failed to open a temporary STRU on one or more ranks";
        // All ranks must leave together before entering the reader's collectives.
        if (open_failures != 0)
        {
            ofs.close();
            if (local_open_failure == 0) std::remove(stru_file.c_str());
            return;
        }
        ofs << std::setprecision(17);
        ofs << "ATOMIC_SPECIES\nHe 4.0026 auto auto\n\nLATTICE_CONSTANT\n" << lat0;
        ofs << "\n\nLATTICE_VECTORS\n4 0 0\n1 3 0\n0.5 0.75 2\n\nATOMIC_POSITIONS\n";
        ofs << units << "\n\nHe\n0.0\n4\n";
        for (int ia = 0; ia < 4; ++ia)
        {
            const double x = 4.0 * input_frac[ia][0] + input_frac[ia][1] + 0.5 * input_frac[ia][2];
            const double y = 3.0 * input_frac[ia][1] + 0.75 * input_frac[ia][2];
            const double z = 2.0 * input_frac[ia][2];
            ofs << x * scale << " " << y * scale << " " << z * scale
                << " m " << moves[ia][0] << " " << moves[ia][1] << " " << moves[ia][2]
                << " v " << 0.01 * (ia + 1) << " 0 0\n";
        }
    }

    const std::vector<std::vector<int> > replicas = {{1, 1, 1}, {2, 3, 2}};
    for (const std::vector<int>& replica : replicas)
    {
        SCOPED_TRACE(units);
        SCOPED_TRACE(replica[0]);
        DomainDecomposition decomp;
        MDCell cell = MDCellReader::read_stru(stru_file, replica, 0.0, comm, decomp);
        const int expected_count = 4 * replica[0] * replica[1] * replica[2];
        EXPECT_EQ(cell.nat(), expected_count);
        EXPECT_EQ(cell.type_atom_counts().size(), 1U);
        if (cell.type_atom_counts().size() == 1U)
        {
            EXPECT_EQ(cell.type_atom_counts()[0], expected_count);
        }
        EXPECT_DOUBLE_EQ(cell.lat0(), lat0);
        const double expected_volume = 24.0 * lat0 * lat0 * lat0 * expected_count / 4;
        EXPECT_NEAR(cell.omega(), expected_volume, tolerance);
        const ModuleBase::Matrix3& lattice = cell.latvec();
        EXPECT_DOUBLE_EQ(lattice.e11, 4.0 * replica[0]);
        EXPECT_DOUBLE_EQ(lattice.e12, 0.0);
        EXPECT_DOUBLE_EQ(lattice.e13, 0.0);
        EXPECT_DOUBLE_EQ(lattice.e21, 1.0 * replica[1]);
        EXPECT_DOUBLE_EQ(lattice.e22, 3.0 * replica[1]);
        EXPECT_DOUBLE_EQ(lattice.e23, 0.0);
        EXPECT_DOUBLE_EQ(lattice.e31, 0.5 * replica[2]);
        EXPECT_DOUBLE_EQ(lattice.e32, 0.75 * replica[2]);
        EXPECT_DOUBLE_EQ(lattice.e33, 2.0 * replica[2]);
        std::vector<int> owners(expected_count, 0);
        for (const LocalAtom& atom : cell.owned_atoms())
        {
            EXPECT_GE(atom.type_index, 0);
            EXPECT_LT(atom.type_index, expected_count);
            // Report invalid IDs without indexing owners or skipping the reduction.
            if (atom.type_index < 0 || atom.type_index >= expected_count) continue;
            const int id = static_cast<int>(atom.type_index);
            ++owners[id];
            const int ia = id % 4;
            const int image = id / 4;
            const int iz = image % replica[2];
            const int iy = (image / replica[2]) % replica[1];
            const int ix = image / (replica[1] * replica[2]);
            const double x = wrapped[ia][0] + ix;
            const double y = wrapped[ia][1] + iy;
            const double z = wrapped[ia][2] + iz;
            EXPECT_NEAR(atom.frac.x, x / replica[0], tolerance);
            EXPECT_NEAR(atom.frac.y, y / replica[1], tolerance);
            EXPECT_NEAR(atom.frac.z, z / replica[2], tolerance);
            EXPECT_NEAR(atom.cart.x, 4.0 * x + y + 0.5 * z, tolerance);
            EXPECT_NEAR(atom.cart.y, 3.0 * y + 0.75 * z, tolerance);
            EXPECT_NEAR(atom.cart.z, 2.0 * z, tolerance);
            EXPECT_EQ(decomp.owner_rank_from_frac(atom.frac), comm.rank());
            EXPECT_EQ(atom.type, 0);
            EXPECT_DOUBLE_EQ(atom.vel.x, 0.01 * (ia + 1));
            EXPECT_DOUBLE_EQ(atom.vel.y, 0.0);
            EXPECT_DOUBLE_EQ(atom.vel.z, 0.0);
            EXPECT_EQ(atom.mbl.x, moves[ia][0]);
            EXPECT_EQ(atom.mbl.y, moves[ia][1]);
            EXPECT_EQ(atom.mbl.z, moves[ia][2]);
            EXPECT_DOUBLE_EQ(atom.mass, 4.0026 / ModuleBase::AU_to_MASS);
        }
        comm.sum(owners.data(), expected_count);
        for (int count : owners) EXPECT_EQ(count, 1);
    }
    EXPECT_EQ(std::remove(stru_file.c_str()), 0);
}

INSTANTIATE_TEST_SUITE_P(CartesianCoordinates,
                         MDCellReaderUnitsTest,
                         testing::Values("Cartesian", "Cartesian_angstrom", "Cartesian_au"));

TEST(MDCellReaderTest, RestartStruPreservesAtomRecordsAcrossRanks)
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
                                        rank));
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
                                        rank));
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
                                        rank));
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
                                        rank));
    }

    ModuleBase::Matrix3 lattice = make_lattice();
    lattice.e11 = 20.0;
    lattice.e22 = 20.0;
    lattice.e33 = 20.0;
    MDCell mdcell;
    mdcell.initialize_from_owned_atoms(lattice,
                  lattice.Inverse(),
                  1.0,
                  1.0,
                  4,
                  owned_atoms,
                  std::vector<std::string>{"A", "B"},
                  std::vector<double>{1.0, 1.0},
                  std::vector<std::int64_t>{2, 2},
                  0.0,
                  ModuleBase::world_comm_domain());
    StruMeta metadata;
    metadata.species.resize(2);
    const std::string output_file = "distributed_mdcell_restart.STRU";
    mdcell::print_stru_file(mdcell, metadata, output_file);

    DomainDecomposition reader_decomp;
    MDCell round_trip = MDCellReader::read_stru(output_file,
                                                 std::vector<int>{1, 1, 1},
                                                 0.0,
                                                 ModuleBase::world_comm_domain(), reader_decomp);
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
