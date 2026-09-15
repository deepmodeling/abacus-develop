#ifndef PARA_TAG_H
#define PARA_TAG_H

#include <string>

namespace Parallel
{

/**
 * @brief Domain tag constants for the parallel communication domains.
 *
 * These tags replace raw string literals to avoid typo-induced runtime
 * failures. They map to the legacy global communicators as follows:
 *   - esolver      -> one esolver instance (intra-image communicator)
 *   - images       -> cross-image communicator (same rank_in_esolver)
 *   - pw           -> POOL_WORLD
 *   - kmesh        -> KP_WORLD
 *   - bsame_kdiff  -> INT_BGROUP
 *   - bdiff_ksame  -> BP_WORLD
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
