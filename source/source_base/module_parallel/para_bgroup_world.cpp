#include "para_bgroup_world.h"

namespace Parallel
{

ParaBgroupWorld::ParaBgroupWorld()
    : ParaWorld("bdiff_ksame"), my_bndgroup_(0), nbndgroup_(1)
{
}

#ifdef __MPI
ParaBgroupWorld::ParaBgroupWorld(const MPI_Comm& intra_comm, const MPI_Comm& inter_comm, int nbndgroup)
    : ParaWorld("bdiff_ksame", intra_comm), inter_comm_(inter_comm), nbndgroup_(nbndgroup)
{
    if (inter_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank(inter_comm, &my_bndgroup_);
    }
}
#endif

} // namespace Parallel
