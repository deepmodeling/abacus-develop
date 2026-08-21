#ifndef BAND_PARALLEL_OUTPUT_H_
#define BAND_PARALLEL_OUTPUT_H_

#include "source_base/matrix.h"

namespace ModuleIO
{

/**
 * @brief Obtain a complete band matrix from replicated or distributed input.
 *
 * Band groups may either hold identical complete column ranges, as in SDFT
 * without BPCG, or complementary contiguous ranges, as in BPCG. An empty
 * matrix is valid when the calculation has no deterministic bands.
 *
 * @param local_matrix Complete or locally owned band columns.
 * @param global_nbands Expected global number of bands.
 * @return Complete matrix on every rank in BP_WORLD.
 */
ModuleBase::matrix gather_band_matrix(const ModuleBase::matrix& local_matrix, const int global_nbands);

} // namespace ModuleIO

#endif // BAND_PARALLEL_OUTPUT_H_
