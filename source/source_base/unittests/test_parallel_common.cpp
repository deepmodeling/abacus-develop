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
    int scalar = sub_rank == root ? 41 : 0;
    Parallel_Common::bcast_int(&scalar, 1, sub, root);
    require(scalar == 41, "bcast_int did not reach every rank from a non-zero root");

    std::vector<double> values(3, sub_rank == root ? 7.5 : -1.0);
    Parallel_Common::bcast_double(values.data(), static_cast<int>(values.size()), sub, root);
    for (const double value : values)
    {
        require(value == 7.5, "bcast_double did not reach every rank from a non-zero root");
    }

    // Allgatherv must accept an all-empty receive buffer. Every rank sends
    // nothing, so the destination vector's data() may be null.
    std::vector<int> empty_receive;
    std::vector<int> zero_counts(static_cast<std::size_t>(sub_size), 0);
    std::vector<int> zero_displacements(static_cast<std::size_t>(sub_size), 0);
    Parallel_Reduce::gatherv_int_all(sub, nullptr, 0, empty_receive.data(),
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
