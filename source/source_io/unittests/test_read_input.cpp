#include "source_base/parallel_global.h"
#include "source_io/module_parameter/read_input.h"

#include <cstdio>
#include <cstdlib>
#include <string>

namespace
{
void finish_reader()
{
    // Print before finalization so an early nonzero exit cannot make the
    // MPI launcher terminate peers before their shutdown markers are written.
    std::fputs("INPUT reader process reached shutdown\n", stderr);
    std::fflush(stderr);
    // Finalization makes a rank-0-only exit fail the test's timeout even if
    // the MPI launcher would otherwise terminate the remaining ranks for us.
    Parallel_Global::finalize_mpi();
}
} // namespace

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        return 2;
    }
    int process_count = 1;
    int thread_count = 1;
    int rank = 0;
    Parallel_Global::read_pal_param(argc, argv, process_count, thread_count, rank);
    // Input-only tests do not create the legacy calculation communicators.
    POOL_WORLD = MPI_COMM_NULL;
    KP_WORLD = MPI_COMM_NULL;
    INT_BGROUP = MPI_COMM_NULL;
    BP_WORLD = MPI_COMM_NULL;
    GRID_WORLD = MPI_COMM_NULL;
    DIAG_WORLD = MPI_COMM_NULL;
    std::atexit(finish_reader);

    const std::string filename = argv[1];
    ModuleIO::ReadInput reader(rank);
    reader.check_ntype_flag = false;
    Parameter input;
    reader.read_parameters(input, filename);
    return 0;
}
