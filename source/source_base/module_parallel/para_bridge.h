#ifndef PARA_BRIDGE_H
#define PARA_BRIDGE_H

#include "para_bdiff_ksame_world.h"
#include "para_kmesh_world.h"
#include "para_world.h"

namespace Parallel
{

/**
 * @brief Temporary bridge: construct a pw-domain ParaWorld from the old
 * global POOL_WORLD (MPI) or as a serial domain (non-MPI).
 *
 * Hides the #ifdef __MPI from call sites so they stay one-liner. Delete
 * this function (and this file) once ParaCollection is wired into driver
 * initialization and callers receive a ParaWorld& from above.
 */
ParaWorld make_pw_world();

/**
 * @brief Temporary bridge: construct a kmesh-domain ParaKmeshWorld from
 * the old globals KP_WORLD / GlobalV::KPAR (MPI) or as a serial domain
 * (non-MPI).
 *
 * Falls back to a serial single-pool domain when MPI is not initialized
 * (e.g. unit tests linked against the MPI-compiled base library) or when
 * there is only one k-point pool, so that no MPI call is made on an
 * unset communicator.
 *
 * @param[in] nkstot  total number of k-points (without spin)
 * @param[in] nspin   number of spin components
 */
ParaKmeshWorld make_kmesh_world(int nkstot, int nspin);

/**
 * @brief Reduce-only overload: construct a kmesh domain for call sites
 * that only need cross-pool reduction (reduce_across_pools etc.) and
 * have no k-point information to pass.
 *
 * The k-point distribution data (nks_pool_, whichpool_, ...) is left
 * empty; calling pool_collection / gather_kvec on the returned object
 * is invalid. Use the (nkstot, nspin) overload when those are needed.
 */
ParaKmeshWorld make_kmesh_world();

/**
 * @brief Temporary bridge: construct a bdiff_ksame-domain
 * ParaBdiffKsameWorld from the old globals INT_BGROUP / BP_WORLD (MPI) or
 * as a serial domain.
 *
 * Falls back to a serial single-band-group domain when MPI is not
 * initialized or the pool layout has not been set up yet (e.g. unit
 * tests), so that no MPI call is made on an unset communicator.
 */
ParaBdiffKsameWorld make_bdiff_ksame_world();

} // namespace Parallel

#endif // PARA_BRIDGE_H
