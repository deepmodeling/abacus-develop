#ifndef CHG_USPP_H
#define CHG_USPP_H

// Stateless double-grid split/merge helpers extracted from Charge_Mixing.
//
// "dgrid" = double grid, the dense/smooth grid pair used by ultrasoft (USPP)
// and PAW calculations to capture the high-frequency tail of the
// augmentation charge beyond the smooth (soft) plane-wave grid. The dense
// grid (npw_dense) is the union of the smooth grid (npw_smooth) and the
// high-frequency tail (npw_dense - npw_smooth).
//
// These functions do not read Charge_Mixing members or PARAM/GlobalV; all
// inputs are passed explicitly. Memory is managed by the caller through
// std::vector, so no new/delete pair is needed and no clean-up function
// exists.

#include <complex>
#include <vector>

namespace module_charge
{

/**
 * @brief Split dense reciprocal-space data into smooth and high-frequency
 *        parts on the USPP double grid.
 *
 * For each spin channel, the first npw_smooth entries of data_d are copied
 * into data_s and the remaining (npw_dense - npw_smooth) entries are copied
 * into data_hf. No aliasing is performed: both output vectors own their
 * storage and must be pre-sized by the caller.
 *
 * @param data_d     dense input, shape [nspin * npw_dense], non-null
 * @param data_s     smooth output, pre-sized to nspin * npw_smooth
 * @param data_hf    high-frequency output, pre-sized to
 *                   nspin * (npw_dense - npw_smooth); zero-size is allowed
 *                   when npw_dense == npw_smooth
 * @param nspin      number of spin channels, >= 1
 * @param npw_smooth smooth grid npw, >= 0
 * @param npw_dense  dense grid npw, >= npw_smooth
 */
void split_dgrid(const std::complex<double>* data_d,
                 std::vector<std::complex<double>>& data_s,
                 std::vector<std::complex<double>>& data_hf,
                 int nspin,
                 int npw_smooth,
                 int npw_dense);

/**
 * @brief Merge smooth and high-frequency parts back into dense reciprocal-
 *        space data. Inverse of split_dgrid.
 *
 * Vectors are not cleared; the caller may reuse them or let them go out of
 * scope. The output data_d must be pre-allocated by the caller with size
 * nspin * npw_dense.
 *
 * @param data_d     dense output, shape [nspin * npw_dense], non-null
 * @param data_s     smooth input, sized to nspin * npw_smooth
 * @param data_hf    high-frequency input, sized to
 *                   nspin * (npw_dense - npw_smooth)
 * @param nspin      number of spin channels, >= 1
 * @param npw_smooth smooth grid npw, >= 0
 * @param npw_dense  dense grid npw, >= npw_smooth
 */
void merge_dgrid(std::complex<double>* data_d,
                 const std::vector<std::complex<double>>& data_s,
                 const std::vector<std::complex<double>>& data_hf,
                 int nspin,
                 int npw_smooth,
                 int npw_dense);

} // namespace module_charge

#endif // CHG_USPP_H
