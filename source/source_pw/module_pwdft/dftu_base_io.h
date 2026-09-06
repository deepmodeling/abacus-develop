#ifndef DFTU_BASE_IO_H
#define DFTU_BASE_IO_H

#include "source_base/matrix.h"
#include "source_estate/occ_matrix.h"

#include <iosfwd>
#include <string>
#include <vector>

class Plus_U_Base;
class UnitCell;

namespace DFTU_BASE
{

/// nested occupation-matrix type used by DFT+U: occ_mat[iat][l][n][spin](m0, m1)
using OccMatData = std::vector<std::vector<std::vector<std::vector<ModuleBase::matrix>>>>;

/// Read the local occupation number matrix from file (rank 0 only).
///
/// The file format matches the output of write_occup_m(). When the file can
/// not be opened, the run quits with an error message that depends on
/// occ_mat_ctrl and init_chg.
void read_occup_m(const UnitCell& ucell,
                  OccupationMatrix& occ,
                  const std::vector<int>& orbital_corr,
                  const int occ_mat_ctrl,
                  const std::string& fn,
                  const std::string& init_chg,
                  int nspin,
                  int npol);

/// Broadcast the local occupation number matrices from rank 0 to all ranks.
///
/// Implemented in dftu_base_io.cpp (only available in MPI builds).
void local_occup_bcast(const UnitCell& ucell,
                       OccupationMatrix& occ,
                       const std::vector<int>& orbital_corr,
                       int nspin,
                       int npol);

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

} // namespace DFTU_BASE

#endif
