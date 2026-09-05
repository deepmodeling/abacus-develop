#include "para_bridge.h"
#include "para_tag.h"

#ifdef __MPI
#include "source_base/global_variable.h"
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

// Reduce-only overload: no k-point distribution data.
ParaKmeshWorld make_kmesh_world()
{
#ifdef __MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    // Any distributed layout (k pools or band groups) may need the
    // world-wide max/min reductions, so build the MPI domain whenever
    // more than one process is running. The sum reduction no-ops for
    // kpar <= 1 on its own.
    if (mpi_initialized && GlobalV::NPROC > 1)
    {
        // Build from globals but skip distribute_kpoints (nkstot=0).
        return ParaKmeshWorld(MPI_COMM_WORLD, GlobalV::KPAR, GlobalV::MY_POOL, 0, 1);
    }
#endif
    return ParaKmeshWorld();
}

// Temporary bridge: construct a kmesh-domain ParaKmeshWorld from the old
// globals. Delete this file once ParaCollection is wired into driver init.
ParaKmeshWorld make_kmesh_world(int nkstot, int nspin)
{
#ifdef __MPI
    // Fall back to a serial single-pool domain when MPI is not initialized
    // (e.g. unit tests linked against the MPI-compiled base library) or when
    // there is only one k-pool, so that no MPI call is made on an unset
    // communicator.
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized && GlobalV::KPAR > 1)
    {
        return ParaKmeshWorld(MPI_COMM_WORLD, GlobalV::KPAR, GlobalV::MY_POOL,
                              nkstot, nspin);
    }
#endif
    return ParaKmeshWorld(nkstot, nspin);
}

// Temporary bridge: construct a bdiff_ksame-domain ParaBdiffKsameWorld from
// the old globals. Delete this file once ParaCollection is wired into driver init.
ParaBdiffKsameWorld make_bdiff_ksame_world()
{
#ifdef __MPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    // NPROC_IN_BNDGROUP stays 0 until divide_pools has run, which also
    // guards unit tests that link the MPI base library without a layout.
    if (mpi_initialized && INT_BGROUP != MPI_COMM_NULL && BP_WORLD != MPI_COMM_NULL
        && GlobalV::NPROC_IN_BNDGROUP > 1)
    {
        int nbndgroup = 1;
        MPI_Comm_size(BP_WORLD, &nbndgroup);
        return ParaBdiffKsameWorld(INT_BGROUP, BP_WORLD, nbndgroup);
    }
#endif
    return ParaBdiffKsameWorld();
}

} // namespace Parallel
