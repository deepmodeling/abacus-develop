#ifndef POS_OP_CSR_H
#define POS_OP_CSR_H

#include <string>

class UnitCell;

namespace ModuleIO
{
namespace detail
{
bool lat_r_nonempty(const int nonzero_num[3]);

/**
 * @brief Assemble one r(R) component file from its temporary block payload.
 * @details Text output uses the same self-describing header and CSR banner as
 *          H(R)/S(R)/dH(R) files, while binary output keeps the three-integer
 *          header followed by the raw block records.
 * @param output_filename Final file name of this Cartesian component.
 * @param payload_filename Temporary file holding the CSR block records.
 * @param step One-based ionic step written into the header.
 * @param nlocal Number of localized basis functions.
 * @param output_R_number Number of Bravais lattice vectors written.
 * @param binary Whether the payload and output are binary.
 * @param append Whether to append to the existing output file.
 * @param axis Cartesian component label ("x", "y" or "z").
 * @param ucell Unit cell written into the text header; must not be null.
 * @param context Caller name used in diagnostics.
 */
void assemble_csr(const std::string& output_filename,
                  const std::string& payload_filename,
                  int step,
                  int nlocal,
                  int output_R_number,
                  bool binary,
                  bool append,
                  const std::string& axis,
                  const UnitCell* ucell,
                  const std::string& context);
} // namespace detail
} // namespace ModuleIO

#endif
