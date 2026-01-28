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
    template<typename T> void reduce_all(T& object) { return; };
    template<typename T> void reduce_pool(T& object) { return; };

    template<>
    void reduce_all<double>(double& object) { return; };
    template<>
    void reduce_pool<double>(double& object) { return; };
    template<>
    void reduce_pool<float>(float& object) { return; };
}
#endif