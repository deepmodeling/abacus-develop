#ifndef ABACUS_VDW_XCNAME_H
#define ABACUS_VDW_XCNAME_H

#include <string>

namespace vdw
{

/**
 * @brief Normalize an XC-functional spelling at the ABACUS/libXC boundary.
 *
 * This function intentionally performs syntax normalization only. It lowercases
 * ASCII names, removes an optional leading "XC_" from each libXC component,
 * trims surrounding whitespace, and writes two-component libXC functionals
 * with ':' as the separator used by libdftd4.
 *
 * Examples:
 *   "PBE"                              -> "pbe"
 *   "XC_HYB_GGA_XC_PBEH"              -> "hyb_gga_xc_pbeh"
 *   "GGA_X_PBE+GGA_C_PBE"             -> "gga_x_pbe:gga_c_pbe"
 *   "XC_GGA_X_PBE+XC_GGA_C_PBE"       -> "gga_x_pbe:gga_c_pbe"
 *
 * Semantic aliases (for example revPBE or damping-specific wB97X aliases)
 * belong to the dispersion-model parameter layer rather than here.
 */
std::string normalize_xc_name(const std::string& input);

} // namespace vdw

#endif // ABACUS_VDW_XCNAME_H
