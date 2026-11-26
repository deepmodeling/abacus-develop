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

    // Payloads are derived from the colour, so the two sub-communicators carry
    // different values. An implementation that broadcast on the wrong
    // communicator (MPI_COMM_WORLD, say) would then deliver the other group's
    // value and fail here instead of passing by coincidence.
    const int group_offset = colour * 1000;

    // A root other than zero, which is the case the RI and RPA stages need.
    const int root = sub_size - 1;
    const int scalar_from_root = group_offset + 41;
    const int scalar_count = 1;
    int scalar = sub_rank == root ? scalar_from_root : 0;
    Parallel_Common::bcast_int(&scalar, scalar_count, sub, root);
    require(scalar == group_offset + 41, "bcast_int did not reach every rank from a non-zero root");

    // Rule 14: name the payload and the element count before the call.
    const double root_value = group_offset + 7.5;
    const double other_value = -1.0;
    const double initial_value = sub_rank == root ? root_value : other_value;
    const int value_count = 3;
    std::vector<double> values(static_cast<std::size_t>(value_count), initial_value);
    Parallel_Common::bcast_double(values.data(), value_count, sub, root);
    for (const double value : values)
    {
        require(value == group_offset + 7.5, "bcast_double did not reach every rank from a non-zero root");
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
    const int send_value = group_offset + 100 + sub_rank;
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
        require(receive[irank] == group_offset + 100 + irank, "gatherv_int_all lost a rank's value");
    }

    // Uneven gather counts: rank 0 receives nothing, the others receive one.
    const bool sends_one = sub_rank > 0;
    const int uneven_send = sends_one ? group_offset + 200 + sub_rank : 0;
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
    // Rule 14: the conditional pointer is named before the call.
    const int* uneven_send_buffer = sends_one ? &uneven_send : nullptr;
    Parallel_Reduce::gatherv_int_all(sub, uneven_send_buffer, uneven_send_count,
                                     uneven_receive.data(), uneven_counts.data(),
                                     uneven_displacements.data());
    for (int irank = 1; irank < sub_size; ++irank)
    {
        require(uneven_receive[static_cast<std::size_t>(irank - 1)] == group_offset + 200 + irank,
                "uneven gatherv misplaced a rank's value");
    }

    // gather_int_all on the same sub-communicator.
    int one_per_rank = group_offset + 300 + sub_rank;
    std::vector<int> gathered(static_cast<std::size_t>(sub_size), 0);
    Parallel_Reduce::gather_int_all(sub, one_per_rank, gathered.data());
    for (int irank = 0; irank < sub_size; ++irank)
    {
        require(gathered[static_cast<std::size_t>(irank)] == group_offset + 300 + irank,
                "gather_int_all lost a rank's value");
    }

    // reduce_max_with_rank: the rank supplying the maximum comes back with it.
    // The maximum sits on rank 1 of every group, and its value carries the
    // colour, so a reduction that leaked onto the world communicator would
    // return the other group's number.
    const double max_value = group_offset + 9.5;
    const double plain_value = 1.0;
    const double local_value = sub_rank == 1 ? max_value : plain_value;
    double best_value = local_value;
    int best_rank = sub_rank;
    Parallel_Reduce::reduce_max_with_rank(sub, best_value, best_rank);
    require(best_value == max_value, "reduce_max_with_rank lost the maximum value");
    require(best_rank == 1, "reduce_max_with_rank reported a rank that did not hold the maximum");

    // The remaining communicator-aware reductions, with per-group inputs so
    // the assertions also show which communicator the result came from.
    int sum_scalar = group_offset + sub_rank;
    Parallel_Reduce::reduce_all(sub, sum_scalar);
    const int expected_sum = sub_size * group_offset + sub_size * (sub_size - 1) / 2;
    require(sum_scalar == expected_sum, "reduce_all(int) summed the wrong communicator");

    const int array_count = 3;
    std::vector<double> array_values(static_cast<std::size_t>(array_count), 0.0);
    for (int i = 0; i < array_count; ++i)
    {
        array_values[static_cast<std::size_t>(i)] = group_offset + sub_rank + i;
    }
    Parallel_Reduce::reduce_all(sub, array_values.data(), array_count);
    for (int i = 0; i < array_count; ++i)
    {
        const double expected = sub_size * group_offset + sub_size * (sub_size - 1) / 2 + sub_size * i;
        require(array_values[static_cast<std::size_t>(i)] == expected,
                "reduce_all(double*, int) summed the wrong communicator");
    }

    int max_of_int = group_offset + 10 * sub_rank;
    Parallel_Reduce::reduce_max(sub, max_of_int);
    require(max_of_int == group_offset + 10 * (sub_size - 1),
            "reduce_max(int) returned the wrong maximum");

    int min_of_int = group_offset + 10 * sub_rank;
    Parallel_Reduce::reduce_min(sub, min_of_int);
    require(min_of_int == group_offset, "reduce_min(int) returned the wrong minimum");

    // unsigned long long has its own instantiations; check them too.
    unsigned long long wide_scalar = static_cast<unsigned long long>(group_offset + sub_rank);
    Parallel_Reduce::reduce_all(sub, wide_scalar);
    const unsigned long long expected_wide
        = static_cast<unsigned long long>(sub_size * group_offset + sub_size * (sub_size - 1) / 2);
    require(wide_scalar == expected_wide, "reduce_all(unsigned long long) summed the wrong communicator");

    std::vector<unsigned long long> wide_array(static_cast<std::size_t>(array_count), 0);
    for (int i = 0; i < array_count; ++i)
    {
        wide_array[static_cast<std::size_t>(i)]
            = static_cast<unsigned long long>(group_offset + sub_rank + i);
    }
    Parallel_Reduce::reduce_all(sub, wide_array.data(), array_count);
    for (int i = 0; i < array_count; ++i)
    {
        const unsigned long long expected
            = static_cast<unsigned long long>(sub_size * group_offset + sub_size * (sub_size - 1) / 2
                                              + sub_size * i);
        require(wide_array[static_cast<std::size_t>(i)] == expected,
                "reduce_all(unsigned long long*, int) summed the wrong communicator");
    }

    unsigned long long wide_max = static_cast<unsigned long long>(group_offset + 10 * sub_rank);
    Parallel_Reduce::reduce_max(sub, wide_max);
    require(wide_max == static_cast<unsigned long long>(group_offset + 10 * (sub_size - 1)),
            "reduce_max(unsigned long long) returned the wrong maximum");

    unsigned long long wide_min = static_cast<unsigned long long>(group_offset + 10 * sub_rank);
    Parallel_Reduce::reduce_min(sub, wide_min);
    require(wide_min == static_cast<unsigned long long>(group_offset),
            "reduce_min(unsigned long long) returned the wrong minimum");

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
