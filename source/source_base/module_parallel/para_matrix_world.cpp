#include "para_matrix_world.h"

namespace Parallel
{

ParaMatrixWorld::ParaMatrixWorld()
    : ParaWorld("matrix")
{
    compute_proc_grid();
}

#ifdef __MPI
ParaMatrixWorld::ParaMatrixWorld(const MPI_Comm& comm)
    : ParaWorld("matrix", comm)
{
    compute_proc_grid();
}
#endif

void ParaMatrixWorld::compute_proc_grid()
{
    const int np = size();
    dim0_ = np;
    while (dim1_ = np / dim0_, dim0_ * dim1_ != np)
    {
        --dim0_;
    }
    coord_row_ = rank() / dim1_;
    coord_col_ = rank() % dim1_;
}

} // namespace Parallel
