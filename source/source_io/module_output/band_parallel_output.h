#ifndef BAND_PARALLEL_OUTPUT_H_
#define BAND_PARALLEL_OUTPUT_H_

#include "source_base/matrix.h"

#include <vector>

namespace ModuleIO
{

/**
 * @brief Describe replicated or contiguous band-group storage.
 *
 * The layout is constructed collectively on BP_WORLD. Replicated storage
 * designates band group 0 as the owner of every global band. Distributed
 * storage assigns each contiguous range to the band group that stores it.
 */
class BandParallelLayout
{
  public:
    /** Construct and validate the layout collectively on BP_WORLD. */
    BandParallelLayout(const int local_nbands, const int global_nbands);

    /** Return whether every band group stores all global bands. */
    bool bands_are_replicated() const;
    /** Return the global band count. */
    int global_nbands() const;
    /** Return the current band group's local band count. */
    int local_nbands() const;
    /** Return the current rank within BP_WORLD. */
    int band_group() const;
    /** Return the current band group's global starting offset. */
    int local_offset() const;
    /** Return the BP_WORLD rank designated as a global band's owner. */
    int owner_group(const int global_band) const;
    /** Return a global band's local index on its owner. */
    int local_index(const int global_band) const;

  private:
    bool bands_are_replicated_ = true;
    int global_nbands_ = 0;
    int local_nbands_ = 0;
    int band_group_ = 0;
    std::vector<int> band_counts_;
    std::vector<int> band_offsets_;

    friend ModuleBase::matrix gather_band_matrix(const ModuleBase::matrix& local_matrix, const int global_nbands);
};

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
