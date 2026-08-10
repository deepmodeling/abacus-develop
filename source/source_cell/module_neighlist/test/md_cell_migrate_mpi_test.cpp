#include <gtest/gtest.h>

#include "source_cell/md_cell.h"

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
                                        rank,
                                        false));
    }
    MDCell mdcell(latvec,
                  latvec.Inverse(),
                  1.0,
                  1.0,
                  2,
                  owned_atoms,
                  std::vector<std::string>(1, "X"),
                  std::vector<double>(1, 1.0),
                  MPI_COMM_WORLD,
                  0.1,
                  0.0);

    ASSERT_EQ(mdcell.mpi_size(), size);
    if (size == 2)
    {
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

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
}
