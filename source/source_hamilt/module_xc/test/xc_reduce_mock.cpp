#include "../../../source_base/parallel_reduce.h"

namespace Parallel_Reduce
{
template<>
void reduce_all<double>(double& object)
{
#ifdef __MPI
    double swap = object;
    MPI_Allreduce(&swap, &object, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
}
}
