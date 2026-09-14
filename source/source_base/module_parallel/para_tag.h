#ifndef PARA_TAG_H
#define PARA_TAG_H

#include <string>

namespace Parallel
{

/**
 * @brief Domain tag constants for the parallel communication domains.
 *
 * Pool terminology (see parallel_comm.cpp):
 *   - k-pool: one of the KPAR groups of processes that share one subset of
 *     k-points. This split happens first and is independent of bndpar.
 *   - band-pool: one of the BNDPAR sub-groups of a k-pool, holding one
 *     band window ("band group").
 *
 * The tags map to the legacy global communicators as follows:
 *   - esolver      -> one esolver instance (intra-image communicator)
 *   - images       -> cross-image communicator (same rank_in_esolver)
 *   - pw           -> POOL_WORLD (one band-pool)
 *   - kmesh        -> KP_WORLD (links k-pools; only valid when the k-pool split is even)
 *   - bsame_kdiff  -> INT_BGROUP (same band group across k-pools)
 *   - bdiff_ksame  -> BP_WORLD (different band groups inside one k-pool)
 *   - rgrid        -> GRID_WORLD
 *   - diag         -> DIAG_WORLD
 *   - matrix       -> matrix domain
 *   - atom         -> atom domain
 */
namespace ParaTag
{
const std::string esolver = "esolver";
const std::string images = "images";
const std::string pw = "pw";
const std::string kmesh = "kmesh";
const std::string bsame_kdiff = "bsame_kdiff";
const std::string bdiff_ksame = "bdiff_ksame";
const std::string rgrid = "rgrid";
const std::string diag = "diag";
const std::string matrix = "matrix";
const std::string atom = "atom";
} // namespace ParaTag

} // namespace Parallel

#endif // PARA_TAG_H
