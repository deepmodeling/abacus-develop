#include "para_bdiff_ksame_world.h"

namespace Parallel
{

ParaBdiffKsameWorld::ParaBdiffKsameWorld()
    : ParaWorld("bdiff_ksame"), my_bndgroup_(0), nbndgroup_(1)
{
}

#ifdef __MPI
ParaBdiffKsameWorld::ParaBdiffKsameWorld(const MPI_Comm& intra_comm, const MPI_Comm& inter_comm, int nbndgroup)
    : ParaWorld("bdiff_ksame", intra_comm), inter_comm_(inter_comm), nbndgroup_(nbndgroup)
{
    if (inter_comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank(inter_comm, &my_bndgroup_);
    }
}
#endif

void ParaBdiffKsameWorld::reduce_across_bdiff_ksame(double& value) const
{
#ifdef __MPI
    if (inter_comm_ == MPI_COMM_NULL || nbndgroup_ <= 1)
    {
        return;
    }
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_SUM, inter_comm_);
#endif
}

} // namespace Parallel
