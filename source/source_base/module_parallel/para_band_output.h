#ifndef PARA_BAND_OUTPUT_H
#define PARA_BAND_OUTPUT_H

#include "source_base/matrix.h"
#include "source_base/module_parallel/para_world.h"

#include <complex>
#include <vector>

namespace Parallel
{

/**
 * @brief Describe replicated or contiguous band-group storage.
 *
 * The layout is constructed collectively on the supplied band domain.
 * Replicated storage designates group 0 as the representative for broadcasting;
 * every group still stores all bands. Distributed
 * storage assigns each contiguous range to the band group that stores it.
 */
class ParaBandOutput
{
  public:
    /**
     * @brief Construct and validate the layout collectively on a valid band domain.
     * @param local_nbands Number of bands stored by this group.
     * @param global_nbands Global band count, identical across the domain.
     * @param band_world Domain connecting corresponding k-point and PW ranks
     * across band groups. Its communicator must outlive this object.
     */
    ParaBandOutput(int local_nbands, int global_nbands, const ParaWorld& band_world);

    /** Return whether every band group stores all global bands. */
    bool bands_are_replicated() const;
    /** Return the global band count. */
    int global_nbands() const;
    /** Return the current band group's local band count. */
    int local_nbands() const;
    /** Return the current rank within the band domain. */
    int band_group() const;
    /** Return the current band group's global starting offset. */
    int local_offset() const;
    /** Return the band-domain rank designated as a global band's owner. */
    int owner_group(const int global_band) const;
    /** Return a global band's local index on its owner. */
    int local_index(const int global_band) const;

    /**
     * @brief Broadcast host data from the group owning a global band.
     *
     * All groups must call in the same band/component order with matching
     * element counts. A serial build validates the band and leaves data unchanged.
     *
     * @param global_band Zero-based global band index.
     * @param data Host buffer to broadcast in place.
     * @param count Number of elements; zero is allowed.
     */
    void bcast_band(int global_band, double* data, int count) const;

    /** @brief Broadcast complex host data with the same contract as the real overload. */
    void bcast_band(int global_band, std::complex<double>* data, int count) const;

    /**
     * @brief Obtain a complete matrix on every rank from local band columns.
     * Replicated input is returned locally without synchronizing replicas.
     * Distributed columns are gathered in band-group order; row counts must
     * match across groups. Empty band sets and empty local shards are valid.
     * @param local_matrix Matrix with local_nbands() columns.
     * @return Complete matrix with global_nbands() columns on every rank.
     */
    ModuleBase::matrix gather_matrix(const ModuleBase::matrix& local_matrix) const;

  private:
    ParaWorld band_world_;
    bool bands_are_replicated_ = true;
    int global_nbands_ = 0;
    int local_nbands_ = 0;
    int band_group_ = 0;
    std::vector<int> band_counts_;
    std::vector<int> band_offsets_;
};

} // namespace Parallel

#endif // PARA_BAND_OUTPUT_H
