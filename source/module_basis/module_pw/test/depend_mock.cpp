#ifdef __MPI
#include "mpi.h"
#endif
#include "depend_mock.h"

namespace GlobalV
{ 
    std::ofstream ofs_running;
}
#ifdef __MPI
//MPI_Comm POOL_WORLD;
//MPI_Comm INTER_POOL = MPI_COMM_NULL;
//MPI_Comm STO_WORLD;
//MPI_Comm PARAPW_WORLD;
//MPI_Comm GRID_WORLD;
//MPI_Comm DIAG_WORLD;
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