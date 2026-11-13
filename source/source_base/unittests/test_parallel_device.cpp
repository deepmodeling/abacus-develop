// parallel_device.cpp: the collectives used by the rVV10 adapter to move grid planes between distributions.
#include "source_base/parallel_device.h"

#include "source_base/parallel_cell.h"
#include "source_base/parallel_comm.h"

#include "gtest/gtest.h"
#include <vector>

TEST(ParallelDevice, AlltoallvExchangesUnevenBlocks)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    MPI_Comm communicator = domain.communicator();
    const int rank = domain.rank();
    MPICommGroup world_group(communicator);
    const int size = world_group.gsize;
    // rank r sends (r + p) % 2 + 1 copies of 100 r + p to rank p
    std::vector<int> send_counts(size);
    std::vector<int> send_offsets(size);
    std::vector<int> recv_counts(size);
    std::vector<int> recv_offsets(size);
    std::vector<double> sent;
    int received_total = 0;
    for (int peer = 0; peer < size; ++peer)
    {
        const double value = 100.0 * rank + peer;
        send_counts[peer] = (rank + peer) % 2 + 1;
        send_offsets[peer] = static_cast<int>(sent.size());
        sent.insert(sent.end(), send_counts[peer], value);
        recv_counts[peer] = send_counts[peer];
        recv_offsets[peer] = received_total;
        received_total += recv_counts[peer];
    }
    std::vector<double> received(received_total);
    Parallel_Common::alltoallv_data(sent.data(),
                                    send_counts.data(),
                                    send_offsets.data(),
                                    received.data(),
                                    recv_counts.data(),
                                    recv_offsets.data(),
                                    communicator);
    for (int peer = 0; peer < size; ++peer)
    {
        const double expected = 100.0 * peer + rank;
        for (int index = 0; index < recv_counts[peer]; ++index)
        {
            const int position = recv_offsets[peer] + index;
            EXPECT_EQ(received[position], expected);
        }
    }
}

TEST(ParallelDevice, AllgatherKeepsRankOrder)
{
    const ModuleBase::CommunicationDomain domain = ModuleBase::world_comm_domain();
    MPI_Comm communicator = domain.communicator();
    const int rank = domain.rank();
    MPICommGroup world_group(communicator);
    const int size = world_group.gsize;
    const int sent[2] = {10 * rank, 10 * rank + 1};
    const int received_count = 2 * size;
    std::vector<int> received(received_count, -1);
    Parallel_Common::allgather_data(sent, 2, received.data(), communicator);
    for (int peer = 0; peer < size; ++peer)
    {
        EXPECT_EQ(received[2 * peer], 10 * peer);
        EXPECT_EQ(received[2 * peer + 1], 10 * peer + 1);
    }
}
