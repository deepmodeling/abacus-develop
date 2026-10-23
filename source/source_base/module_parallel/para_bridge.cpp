#include "para_bridge.h"
#include "para_tag.h"

#ifdef __MPI
#include "source_base/parallel_comm.h"
#endif

namespace Parallel
{

// Temporary bridge: construct a pw-domain ParaWorld from the old globals.
// Delete this file once ParaCollection is wired into driver initialization.
ParaWorld make_pw_world()
{
#ifdef __MPI
    return ParaWorld::make_mpi(ParaTag::pw, POOL_WORLD);
#else
    return ParaWorld::serial(ParaTag::pw);
#endif
}

ParaWorld make_band_world()
{
#ifdef __MPI
    return ParaWorld::make_mpi(ParaTag::bdiff_ksame, BP_WORLD);
#else
    return ParaWorld::serial(ParaTag::bdiff_ksame);
#endif
}

} // namespace Parallel
