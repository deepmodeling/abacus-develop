#ifndef LIBXC_R_SUFFIX_H
#define LIBXC_R_SUFFIX_H

#include <string>
#include <vector>

namespace ModuleIO
{

/// @brief Build the channel suffixes appended to the cube file names of one Libxc quantity.
///
/// The suffix encodes the Libxc channel index tuple, which differs per quantity order:
///   - rho / vrho       : s1, s2 (real spin indices)
///   - exc / amag       : empty suffix (scalar quantity)
///   - sigma / vsigma / v2rho2 : (00), (01), (11) for nspin=2
///   - v2rhosigma       : (0_00), (0_01), (0_11), (1_00), (1_01), (1_11)
///   - v2sigma2         : (00_00), (00_01), (00_11), (01_01), (01_11), (11_11)
///   - v3rho3           : (0_0_0), (0_0_1), (0_1_1), (1_1_1)
///   - v4rho4           : (0_0_0_0), (0_0_0_1), (0_0_1_1), (0_1_1_1), (1_1_1_1)
///   - v3rho2sigma      : 9 channels, v3rhosigma2 : 12, v3sigma3 : 10
///   - v4rho3sigma      : 12 channels, v4rho2sigma2 : 18, v4rhosigma3 : 20, v4sigma4 : 15
///
/// @param data_name Quantity label used in the file name (e.g. "v2rhosigma").
/// @param nspin Number of spin channels of the calculation, 1 or 2.
/// @return One suffix per Libxc channel, in Libxc packing order.
/// @throws std::invalid_argument if nspin is not 1 or 2, or data_name is unknown.
std::vector<std::string> libxc_r_channel_suffixes(const std::string& data_name, int nspin);

} // namespace ModuleIO

#endif // LIBXC_R_SUFFIX_H
