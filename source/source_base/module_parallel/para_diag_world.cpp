#include "para_diag_world.h"

namespace Parallel
{

ParaDiagWorld::ParaDiagWorld()
    : ParaWorld("diag"), dcolor_(0)
{
}

#ifdef __MPI
ParaDiagWorld::ParaDiagWorld(const MPI_Comm& comm, int dcolor)
    : ParaWorld("diag", comm), dcolor_(dcolor)
{
}
#endif

} // namespace Parallel
