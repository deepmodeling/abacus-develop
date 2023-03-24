#ifdef __MPI
#include "mpi.h"
#endif
#include "depend_mock.h"

namespace GlobalV
{ 
    std::ofstream ofs_running;
}
#ifdef __MPI
MPI_Comm POOL_WORLD;
namespace Parallel_Reduce
{
    void reduce_double_all(double &object){return;};
    void reduce_double_pool(double &object){return;};
    void reduce_double_pool(float &object){return;};
}
#endif