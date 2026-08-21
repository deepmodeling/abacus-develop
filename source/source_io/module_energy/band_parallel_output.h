#ifndef BAND_PARALLEL_OUTPUT_H_
#define BAND_PARALLEL_OUTPUT_H_

#include "source_base/matrix.h"

namespace ModuleIO
{

/**
 * @brief Reconstruct global band columns from a contiguous band distribution.
 *
 * @param local_matrix Matrix whose columns are locally owned bands.
 * @param global_nbands Expected global number of bands.
 * @return Matrix with all global band columns on every rank in BP_WORLD.
 */
ModuleBase::matrix gather_band_matrix(const ModuleBase::matrix& local_matrix, const int global_nbands);

} // namespace ModuleIO

#endif // BAND_PARALLEL_OUTPUT_H_
