#ifndef DFTU_OUTPUT_H
#define DFTU_OUTPUT_H

#include <iosfwd>
#include <string>

class Plus_U_Base;
class UnitCell;

namespace dftu_io
{

/// Output DFT+U information (Hubbard U/J, local occupation matrices) to the
/// running log and, when out_chg is set, to the dm_onsite.txt file.
///
/// Extracted from Plus_U_Base::output as a free function so that IO logic is
/// decoupled from the Plus_U_Base class. The function only reads the
/// Plus_U_Base state via public accessors; no friend declaration needed.
void output(const Plus_U_Base& dftu,
            const UnitCell& ucell,
            bool out_chg,
            const std::string& global_out_dir,
            int nspin,
            int npol);

/// Write local occupation matrices to the given stream.
///
/// Extracted from Plus_U_Base::write_occup_m. When diag is true, eigenvalues
/// and magnetism are also printed; otherwise only raw matrix elements.
/// Caller is responsible for opening/closing the stream.
void write_occup_m(const Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   std::ofstream& ofs,
                   bool diag,
                   int nspin,
                   int npol);

} // namespace dftu_io

#endif
