#include "source_base/communication_domain.h"
#include "source_base/global_variable.h"
#include "source_base/parallel_comm.h"
#include "source_base/parallel_grid.h"

#include "gtest/gtest.h"

namespace legacy_global = GlobalV;
#include <vector>

TEST(CommunicationDomainTest, ReportsDefaultAndWorldDomains)
{
    const ModuleBase::CommunicationDomain local_domain;
    EXPECT_EQ(local_domain.rank(), 0);
    EXPECT_EQ(local_domain.size(), 1);

    const ModuleBase::CommunicationDomain world_domain = ModuleBase::world_communication_domain();
    EXPECT_EQ(world_domain.communicator(), MPI_COMM_WORLD);
    EXPECT_GE(world_domain.rank(), 0);
    EXPECT_LT(world_domain.rank(), world_domain.size());

    const ModuleBase::CommunicationDomain null_domain(MPI_COMM_NULL);
    EXPECT_EQ(null_domain.communicator(), MPI_COMM_NULL);
    EXPECT_EQ(null_domain.rank(), 0);
    EXPECT_EQ(null_domain.size(), 1);
}

TEST(MPICommGroupTest, DividesWorldIntoEvenGroups)
{
    const ModuleBase::CommunicationDomain world_domain = ModuleBase::world_communication_domain();
    MPICommGroup group(MPI_COMM_WORLD);
    EXPECT_EQ(group.gsize, world_domain.size());
    EXPECT_EQ(group.grank, world_domain.rank());

    const int group_count = world_domain.size() > 1 ? 2 : 1;
    group.divide_group_comm(group_count);

    EXPECT_TRUE(group.is_even);
    EXPECT_EQ(group.ngroups, group_count);
    EXPECT_EQ(group.nprocs_in_group, world_domain.size() / group_count);
    EXPECT_EQ(group.my_group, world_domain.rank() / group.nprocs_in_group);
    EXPECT_EQ(group.rank_in_group, world_domain.rank() % group.nprocs_in_group);
    EXPECT_NE(group.group_comm, MPI_COMM_NULL);
    EXPECT_NE(group.inter_comm, MPI_COMM_NULL);
}

TEST(ParallelGridTest, BroadcastsAndReducesDistributedGrid)
{
    const ModuleBase::CommunicationDomain world_domain = ModuleBase::world_communication_domain();
    const int nx = 2;
    const int ny = 1;
    const int nz = world_domain.size();
    const int local_nz = 1;
    const int local_size = nx * ny * local_nz;

    legacy_global::KPAR = 1;
    legacy_global::MY_POOL = 0;
    legacy_global::NPROC = world_domain.size();
    legacy_global::MY_RANK = world_domain.rank();
    legacy_global::NPROC_IN_POOL = world_domain.size();
    legacy_global::RANK_IN_POOL = world_domain.rank();
    legacy_global::RANK_IN_BPGROUP = world_domain.rank();
    POOL_WORLD = MPI_COMM_WORLD;
    INT_BGROUP = MPI_COMM_WORLD;
    KP_WORLD = MPI_COMM_NULL;

    Parallel_Grid grid;
    grid.init(nx, ny, nz, local_nz, local_size, nz, 1, world_domain.size());
    EXPECT_EQ(grid.get_nx(), nx);
    EXPECT_EQ(grid.get_ny(), ny);
    EXPECT_EQ(grid.get_nz(), nz);
    EXPECT_EQ(grid.get_nrxx(), local_size);

    std::vector<double> global_data(nx * ny * nz, 0.0);
    if (world_domain.rank() == 0)
    {
        for (std::size_t i = 0; i < global_data.size(); ++i)
        {
            global_data[i] = static_cast<double>(i + 1);
        }
    }

    std::vector<double> local_data(local_size, 0.0);
    grid.bcast(global_data.data(), local_data.data(), world_domain.rank(), false);
    for (int ix = 0; ix < nx; ++ix)
    {
        EXPECT_DOUBLE_EQ(local_data[ix], static_cast<double>(ix * nz + world_domain.rank() + 1));
    }

    std::vector<double> stochastic_data(local_size, 0.0);
    grid.bcast(global_data.data(), stochastic_data.data(), world_domain.rank(), true);
    EXPECT_EQ(stochastic_data, local_data);

    grid.reduce_across_pools(local_data.data());
    for (int ix = 0; ix < nx; ++ix)
    {
        EXPECT_DOUBLE_EQ(local_data[ix], static_cast<double>(ix * nz + world_domain.rank() + 1));
    }

    std::vector<double> reduced_data(nx * ny * nz, 0.0);
    grid.reduce(reduced_data.data(), local_data.data(), true);
    if (world_domain.rank() == 0)
    {
        for (std::size_t i = 0; i < reduced_data.size(); ++i)
        {
            EXPECT_DOUBLE_EQ(reduced_data[i], static_cast<double>(i + 1));
        }
    }
}
