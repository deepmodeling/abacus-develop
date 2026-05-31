// Stub definitions for global symbols used by ABACUS linked code.
// The benchmark does not exercise the code paths that use these globals,
// but the linker needs them resolved.
#include <mpi.h>
#include <string>

// MPI communicator globals (defined in parallel_global.cpp, used by parallel_reduce.cpp)
MPI_Comm POOL_WORLD  = MPI_COMM_NULL;
MPI_Comm DIAG_WORLD  = MPI_COMM_NULL;
MPI_Comm GRID_WORLD  = MPI_COMM_NULL;
MPI_Comm KP_WORLD    = MPI_COMM_NULL;
MPI_Comm INT_BGROUP  = MPI_COMM_NULL;
MPI_Comm BP_WORLD    = MPI_COMM_NULL;

// Stub for Global_File::close_all_log (used by tool_quit.cpp)
namespace ModuleBase {
namespace Global_File {
void close_all_log(int, bool, const std::string&) {}
} // namespace Global_File

// Stub for GlobalFunc::MAKE_DIR (used by global_file.cpp)
namespace GlobalFunc {
void MAKE_DIR(const std::string&) {}
} // namespace GlobalFunc
} // namespace ModuleBase
