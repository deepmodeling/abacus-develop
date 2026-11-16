#ifndef ABACUS_LR_STATE_IO_H
#define ABACUS_LR_STATE_IO_H

#include <complex>
#include <string>

namespace LR
{
/**
 * @brief Write existing LR eigenstate files with full double precision.
 * @param[in] directory Existing output directory including its trailing separator.
 * @param[in] label Spin label used by existing energy/amplitude filenames.
 * @param[in] suffix Existing optional ionic-step filename suffix.
 * @param[in] rank Zero-based rank; energy output is restricted to rank zero.
 * @param[in] energies Excitation energies in Ry, one per state.
 * @param[in] amplitudes State-major local amplitude array.
 * @param[in] local_size Entries per state on this rank.
 * @param[in] state_count Number of states to write.
 * @details Stream precision remains 17 after the decimal in scientific notation.
 * Output calls execute with NDEBUG as well as in assertion-enabled builds.
 */
template <typename T>
void write_excitation_states(
    const std::string &directory,
    const std::string &label,
    const std::string &suffix,
    int rank,
    const double *energies,
    const T *amplitudes,
    int local_size,
    int state_count);
} // namespace LR

#endif
