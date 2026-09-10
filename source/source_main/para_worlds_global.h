#ifndef PARA_WORLDS_GLOBAL_H
#define PARA_WORLDS_GLOBAL_H

namespace Parallel
{

// Only references and pointers to ParaCollection appear in this header, so a
// forward declaration keeps the include dependency minimal; consumers that
// need the complete type include source_base/module_parallel/para_collection.h.
class ParaCollection;

/**
 * @brief Initialize the process-wide ParaCollection exactly once.
 *
 * Builds the image-level decomposition: MPI_COMM_WORLD is split by image id
 * into one esolver_world per image, and an images_world inter-communicator
 * connecting ranks with the same rank_in_esolver across images. Only these
 * two domains are registered for now; finer domains (pools/diag/rgrid) are
 * still produced by the legacy Parallel_Global path and will be migrated in
 * later steps.
 *
 * With nimage = 1 the esolver_world is a duplicate of MPI_COMM_WORLD, so the
 * existing downstream decomposition on MPI_COMM_WORLD stays bit-identical.
 *
 * Must be called once during startup (from the driver). Calling it again is a
 * programming error and aborts via WARNING_QUIT.
 *
 * @param[in] nproc    total MPI size of MPI_COMM_WORLD
 * @param[in] my_rank  rank of this process in MPI_COMM_WORLD
 * @param[in] nimage   number of images to split into (must be >= 1)
 * @return reference to the initialized collection
 */
ParaCollection& init_global_para_worlds(int nproc, int my_rank, int nimage);

/**
 * @brief Read-only access to the process-wide ParaCollection.
 *
 * Aborts via WARNING_QUIT if init_global_para_worlds has not been called yet.
 */
const ParaCollection& global_para_worlds();

/**
 * @brief Reset the global collection. Test-only; not for production use.
 */
void reset_global_para_worlds_for_test();

} // namespace Parallel

#endif // PARA_WORLDS_GLOBAL_H
