#include <fstream>
//memory.cpp depends on GlobalV::ofs_running and reduce_all
//GPU depends on reduce_pool
#ifdef __MPI
#include "mpi.h"

extern MPI_Comm POOL_WORLD;
namespace Parallel_Reduce
{
    template<typename T> void reduce_all(T& object);
    template<typename T> void reduce_pool(T& object);
    template<typename T> void reduce_min_pool(const int& nproc_in_pool, T& v);
}
#endif

namespace GlobalV
{ 
    extern std::ofstream ofs_running;
}