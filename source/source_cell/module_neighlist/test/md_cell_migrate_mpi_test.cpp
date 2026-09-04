#include <gtest/gtest.h>

#include "source_cell/md_cell.h"
#include "source_base/parallel_cell.h"

#include <mpi.h>

#include <cstdint>
#include <mpi.h>
#include <string>
#include <vector>

namespace
{
void ensure_mpi_initialized()
{
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        int provided = 0;
        MPI_Init_thread(NULL, NULL, MPI_THREAD_SINGLE, &provided);
    }
}

ModuleBase::Matrix3 make_lattice()
{
    ModuleBase::Matrix3 latvec;
    latvec.e11 = 1.0;
    latvec.e22 = 1.0;
    latvec.e33 = 1.0;
    return latvec;
}
}

TEST(MdCellMigrateMpiTest, AtomCrossingDomainMigratesToNewOwner)
{
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ASSERT_GE(size, 2);

    const ModuleBase::Matrix3 latvec = make_lattice();
    std::vector<LocalAtom> owned_atoms;
    if (rank < 2)
    {
        const ModuleBase::Vector3<double> frac(rank == 0 ? 0.2 : 0.7, 0.2, 0.2);
        owned_atoms.push_back(LocalAtom(frac,
                                        frac,
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                        ModuleBase::Vector3<int>(1, 1, 1),
                                        1.0,
                                        0,
                                        rank,
                                        rank));
    }
    MDCell mdcell;
    mdcell.initialize_from_owned_atoms(latvec,
                  latvec.Inverse(),
                  1.0,
                  1.0,
                  2,
                  owned_atoms,
                  std::vector<std::string>(1, "X"),
                  std::vector<double>(1, 1.0),
                  std::vector<std::int64_t>(1, 2),
                  0.1,
                  0.0,
                  ModuleBase::world_communication_domain());

    ASSERT_EQ(mdcell.mpi_size(), size);
    if (size == 2)
    {
        ASSERT_EQ(mdcell.nlocal(), 1);
        mdcell.mutable_owned_atoms()[0].vel.x = static_cast<double>(rank + 1);
        mdcell.mutable_owned_atoms()[0].force.y = static_cast<double>(rank + 3);
        mdcell.migrate_owned_atoms();
        ASSERT_EQ(mdcell.nlocal(), 1);
        EXPECT_EQ(mdcell.owned_atoms()[0].owner_rank, rank);
        EXPECT_EQ(mdcell.owned_atoms()[0].vel.x, static_cast<double>(rank + 1));
        EXPECT_EQ(mdcell.owned_atoms()[0].force.y, static_cast<double>(rank + 3));

        if (rank == 0 && mdcell.nlocal() == 1)
        {
            mdcell.mutable_owned_atoms()[0].cart.x = 0.8;
        }
        if (rank == 1 && mdcell.nlocal() == 1)
        {
            mdcell.mutable_owned_atoms()[0].cart.x = 0.3;
        }
        mdcell.migrate_owned_atoms();

        long long local_count = mdcell.nlocal();
        long long global_count = 0;
        MPI_Allreduce(&local_count, &global_count, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        EXPECT_EQ(global_count, 2);

        for (int i = 0; i < mdcell.nlocal(); ++i)
        {
            EXPECT_EQ(mdcell.owned_atoms()[static_cast<std::size_t>(i)].owner_rank, rank);
        }
    }
}

TEST(MdCellMigrateMpiTest, GhostForcesReturnToOwners)
{
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ASSERT_EQ(size, 2);

    const ModuleBase::Matrix3 latvec = make_lattice();
    std::vector<LocalAtom> owned_atoms;
    const ModuleBase::Vector3<double> frac(rank == 0 ? 0.2 : 0.7, 0.2, 0.2);
    owned_atoms.push_back(LocalAtom(frac,
                                    frac,
                                    ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                    ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                    ModuleBase::Vector3<int>(1, 1, 1),
                                    1.0,
                                    0,
                                    rank,
                                    rank));
    MDCell mdcell;
    mdcell.initialize_from_owned_atoms(latvec,
                  latvec.Inverse(),
                  1.0,
                  1.0,
                  2,
                  owned_atoms,
                  std::vector<std::string>(1, "X"),
                  std::vector<double>(1, 1.0),
                  std::vector<std::int64_t>(1, 2),
                  0.6,
                  0.0,
                  ModuleBase::world_communication_domain());

    long long local_copies[2] = {0, 0};
    for (std::size_t iat = 0; iat < mdcell.ghost_atoms().size(); ++iat)
    {
        ++local_copies[mdcell.ghost_atoms()[iat].owner_rank];
    }
    long long global_copies[2] = {0, 0};
    MPI_Allreduce(local_copies, global_copies, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    for (std::size_t iat = 0; iat < mdcell.mutable_ghost_atoms().size(); ++iat)
    {
        LocalAtom& ghost = mdcell.mutable_ghost_atoms()[iat];
        const double value = static_cast<double>(ghost.owner_rank + 1);
        ghost.force.set(value, 2.0 * value, 3.0 * value);
    }
    mdcell.accumulate_ghost_forces();

    ASSERT_EQ(mdcell.nlocal(), 1);
    const double expected = static_cast<double>(global_copies[rank] * (rank + 1));
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[0].force.x, expected);
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[0].force.y, 2.0 * expected);
    EXPECT_DOUBLE_EQ(mdcell.owned_atoms()[0].force.z, 3.0 * expected);
}

TEST(MdCellMigrateMpiTest, SkinUpdatesFixedGhostLayoutBeforeRebuild)
{
    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ASSERT_EQ(size, 2);

    const ModuleBase::Matrix3 latvec = make_lattice();
    const ModuleBase::Vector3<double> frac(rank == 0 ? 0.20 : 0.70, 0.20, 0.20);
    std::vector<LocalAtom> owned_atoms(1, LocalAtom(frac,
                                                     frac,
                                                     ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                                     ModuleBase::Vector3<double>(0.0, 0.0, 0.0),
                                                     ModuleBase::Vector3<int>(1, 1, 1),
                                                     1.0,
                                                     0,
                                                     rank,
                                                     rank));
    MDCell mdcell;
    mdcell.initialize_from_owned_atoms(latvec,
                  latvec.Inverse(),
                  1.0,
                  1.0,
                  2,
                  owned_atoms,
                  std::vector<std::string>(1, "X"),
                  std::vector<double>(1, 1.0),
                  std::vector<std::int64_t>(1, 2),
                  0.1,
                  0.2,
                  ModuleBase::world_communication_domain());

    mdcell.prepare_neighbors();
    mdcell.mutable_owned_atoms()[0].frac.x += rank == 0 ? 0.05 : -0.05;
    mdcell.mutable_owned_atoms()[0].cart = mdcell.mutable_owned_atoms()[0].frac * latvec;
    mdcell.prepare_neighbors();

    ASSERT_EQ(mdcell.nlocal(), 1);
    for (std::size_t i = 0; i < mdcell.ghost_atoms().size(); ++i)
    {
        const LocalAtom& ghost = mdcell.ghost_atoms()[i];
        const double expected_frac = ghost.owner_rank == 0 ? 0.25 : 0.65;
        EXPECT_DOUBLE_EQ(ghost.frac.x, expected_frac);
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
