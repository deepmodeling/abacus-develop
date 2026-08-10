#include <gtest/gtest.h>

#include "source_cell/distributed_mdcell_reader.h"
#include "source_cell/md_cell.h"
#include "source_base/constants.h"
#include "source_base/communication_domain.h"
#include "source_base/global_variable.h"
#include "source_cell/module_neighlist/domain_decomposition.h"

#include <fstream>
#include <mpi.h>
#include <set>
#include <string>

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
