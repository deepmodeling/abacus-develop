// The sub-communicator helpers are only meaningful with MPI, so this test is
// built and registered for the MPI configuration only.
#include "source_base/parallel_common.h"
#include "source_base/parallel_reduce.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#ifdef __MPI
#include <mpi.h>
#endif

namespace
{
void require(const bool condition, const char* message)
{
    if (!condition)
    {
        std::fprintf(stderr, "%s\n", message);
        std::exit(EXIT_FAILURE);
    }
}
} // namespace

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);

    // Split into two groups so the helpers are exercised on a communicator
    // that is neither MPI_COMM_WORLD nor rooted at world rank 0.
    int world_rank = 0;
    int world_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    require(world_size >= 4, "run this test with at least four ranks");

    const int colour = world_rank < world_size / 2 ? 0 : 1;
    MPI_Comm sub = MPI_COMM_NULL;
    MPI_Comm_split(MPI_COMM_WORLD, colour, world_rank, &sub);

    int sub_rank = 0;
    int sub_size = 0;
    MPI_Comm_rank(sub, &sub_rank);
    MPI_Comm_size(sub, &sub_size);
    require(Parallel_Common::my_rank(sub) == sub_rank, "my_rank disagrees with MPI_Comm_rank");
    require(Parallel_Common::nproc(sub) == sub_size, "nproc disagrees with MPI_Comm_size");

    // A root other than zero, which is the case the RI and RPA stages need.
    const int root = sub_size - 1;
    const int scalar_from_root = 41;
    const int scalar_count = 1;
    int scalar = sub_rank == root ? scalar_from_root : 0;
    Parallel_Common::bcast_int(&scalar, scalar_count, sub, root);
    require(scalar == 41, "bcast_int did not reach every rank from a non-zero root");

    // Rule 14: name the payload and the element count before the call.
    const double root_value = 7.5;
    const double other_value = -1.0;
    const double initial_value = sub_rank == root ? root_value : other_value;
    const int value_count = 3;
    std::vector<double> values(static_cast<std::size_t>(value_count), initial_value);
    Parallel_Common::bcast_double(values.data(), value_count, sub, root);
    for (const double value : values)
    {
        require(value == 7.5, "bcast_double did not reach every rank from a non-zero root");
    }

    // Allgatherv must accept an all-empty receive buffer. Every rank sends
    // nothing, so the destination vector's data() may be null.
    std::vector<int> empty_receive;
    std::vector<int> zero_counts(static_cast<std::size_t>(sub_size), 0);
    std::vector<int> zero_displacements(static_cast<std::size_t>(sub_size), 0);
    const int no_send_count = 0;
    Parallel_Reduce::gatherv_int_all(sub, nullptr, no_send_count, empty_receive.data(),
                                     zero_counts.data(), zero_displacements.data());
    require(empty_receive.empty(), "all-empty gatherv produced data");

    // And it must still gather when there is something to receive.
    const int send_value = 100 + sub_rank;
    std::vector<int> receive(static_cast<std::size_t>(sub_size), 0);
    std::vector<int> counts(static_cast<std::size_t>(sub_size), 1);
    std::vector<int> displacements(static_cast<std::size_t>(sub_size), 0);
    for (int irank = 0; irank < sub_size; ++irank)
    {
        displacements[irank] = irank;
    }
    Parallel_Reduce::gatherv_int_all(sub, &send_value, 1, receive.data(),
                                     counts.data(), displacements.data());
    for (int irank = 0; irank < sub_size; ++irank)
    {
        require(receive[irank] == 100 + irank, "gatherv_int_all lost a rank's value");
    }

    // Uneven gather counts: rank 0 receives nothing, the others receive one.
    const bool sends_one = sub_rank > 0;
    const int uneven_send = sends_one ? 200 + sub_rank : 0;
    const int uneven_send_count = sends_one ? 1 : 0;
    std::vector<int> uneven_counts(static_cast<std::size_t>(sub_size), 0);
    std::vector<int> uneven_displacements(static_cast<std::size_t>(sub_size), 0);
    int running_offset = 0;
    for (int irank = 0; irank < sub_size; ++irank)
    {
        uneven_counts[static_cast<std::size_t>(irank)] = irank > 0 ? 1 : 0;
        uneven_displacements[static_cast<std::size_t>(irank)] = running_offset;
        running_offset += uneven_counts[static_cast<std::size_t>(irank)];
    }
    std::vector<int> uneven_receive(static_cast<std::size_t>(running_offset), -1);
    Parallel_Reduce::gatherv_int_all(sub, sends_one ? &uneven_send : nullptr,
                                     uneven_send_count, uneven_receive.data(),
                                     uneven_counts.data(), uneven_displacements.data());
    for (int irank = 1; irank < sub_size; ++irank)
    {
        require(uneven_receive[static_cast<std::size_t>(irank - 1)] == 200 + irank,
                "uneven gatherv misplaced a rank's value");
    }

    // gather_int_all on the same sub-communicator.
    int one_per_rank = 300 + sub_rank;
    std::vector<int> gathered(static_cast<std::size_t>(sub_size), 0);
    Parallel_Reduce::gather_int_all(sub, one_per_rank, gathered.data());
    for (int irank = 0; irank < sub_size; ++irank)
    {
        require(gathered[static_cast<std::size_t>(irank)] == 300 + irank,
                "gather_int_all lost a rank's value");
    }

    // reduce_max_with_rank: the rank supplying the maximum comes back with it. On a
    // tie MPI may report any contributing rank, so only check the value there.
    const double local_value = (sub_rank == 0 || sub_rank == 1) ? 9.5 : 1.0;
    const int local_rank = sub_rank;
    double best_value = local_value;
    int best_rank = local_rank;
    Parallel_Reduce::reduce_max_with_rank(sub, best_value, best_rank);
    require(best_value == 9.5, "reduce_max_with_rank lost the maximum value");
    const bool tied = sub_size > 1;
    if (tied)
    {
        require(best_rank == 0 || best_rank == 1,
                "reduce_max_with_rank reported a rank that did not hold the maximum");
    }

    // The interface must also exist in a serial build. The header keeps its
    // declarations unconditional, so this call compiles whether or not __MPI
    // is defined; only the implementation behind it changes.
    require(Parallel_Common::nproc(sub) >= 1, "nproc must report at least one rank");

    Parallel_Common::barrier(sub);
    MPI_Comm_free(&sub);
    MPI_Finalize();
    if (world_rank == 0)
    {
        std::printf("sub-communicator helper tests passed\n");
    }
#else
    (void)argc;
    (void)argv;
    // The helpers are declared only under __MPI: without it there is no
    // communicator type to pass, so there is nothing to exercise here.
    std::printf("sub-communicator helper tests skipped (no MPI)\n");
#endif
    return EXIT_SUCCESS;
}
