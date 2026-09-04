#ifndef PARA_BRIDGE_H
#define PARA_BRIDGE_H

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

} // namespace Parallel

#endif // PARA_BRIDGE_H
