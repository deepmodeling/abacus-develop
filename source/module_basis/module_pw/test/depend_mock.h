#include <fstream>
//memory.cpp depends on GlobalV::ofs_running and reduce_double_all
//GPU depends on reduce_double_pool
#ifdef __MPI
#include "mpi.h"

extern MPI_Comm POOL_WORLD;
namespace Parallel_Reduce
{
    void reduce_double_all(double &object);
    void reduce_double_pool(double &object);
    void reduce_double_pool(float &object);
}
#endif

namespace GlobalV
{ 
    extern std::ofstream ofs_running;
}