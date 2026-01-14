#ifdef __MPI
#include "mpi.h"
#endif
#include "depend_mock.h"

namespace GlobalV
{ 
    std::ofstream ofs_running;
    int NPROC_IN_POOL = 1;
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
    template<> void reduce_pool<float>(float& object) { return; };

    template <typename T> void reduce_min_pool(T& object) { return; }
    template <typename T> void reduce_max_pool(T* object, const int n) { return; }

    template<> void reduce_min_pool<double>(double& val) { return; }
    template<> void reduce_max_pool<int>(int* object, const int n) { return; }
}
#endif