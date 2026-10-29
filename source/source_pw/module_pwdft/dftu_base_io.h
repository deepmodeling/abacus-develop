#ifndef DFTU_BASE_IO_H
#define DFTU_BASE_IO_H

#include "source_base/matrix.h"
#include "source_estate/occ_matrix.h"
#include "source_pw/module_pwdft/dftu_types.h"

#include <iosfwd>
#include <string>
#include <vector>

class Plus_U_Base;
class UnitCell;

namespace DFTU_BASE
{

/// nested occupation-matrix type used by DFT+U: occ_mat[iat][l][spin](m0, m1)
using OccMatData = std::vector<std::vector<std::vector<ModuleBase::matrix>>>;

/// Tell whether the given ionic step produces a numbered file.
///
/// @param istep ionic-step index, starting from 0
/// @param cfg frequency configuration
/// @return true when out_freq_ion > 0 and istep is an output ionic step
bool is_ion_step_output_step(int istep, const OccmatOutputCfg& cfg);

/// Tell whether the current electronic step is recorded inside the numbered file.
///
/// @param iter electronic-iteration index, starting from 1
/// @param conv_esolver whether the electronic SCF is converged at this step
/// @param cfg frequency configuration
/// @return true when iter % out_freq_elec == 0, iter == scf_nmax, or converged
bool is_elec_snapshot_trigger(int iter,
                              bool conv_esolver,
                              const OccmatOutputCfg& cfg);

/// Build the name of the per-ionic-step occupation-matrix file.
///
/// @param out_dir output directory (including the trailing separator)
/// @param istep ionic-step index, starting from 0
/// @return occ_matg{istep+1}.txt
std::string gen_ion_step_occ_mat_filename(const std::string& out_dir, int istep);

/// Return the full path of the first candidate file that exists in @p dir.
/// Returns an empty string if none of the candidates exist.
/// Only the calling process probes the filesystem; callers must arrange
/// MPI broadcast of the result if other ranks need it.
std::string find_first_existing_file(const std::string& dir,
                                     const std::vector<std::string>& candidates);

/// Append one electronic-step section to the per-ionic-step file.
///
/// The section records the electronic-step index, the configured charge-
/// density convergence threshold (scf_thr) and the actual residual (drho)
/// of the current electronic step, followed by the occupation matrices and
/// per-atom magnetism. When occmat_ready is false, an "N/A" placeholder is
/// recorded instead of the matrix body (the PW path has no matrix at
/// istep 0 / iter 1 unless it was loaded from file). The caller invokes this
/// function only after drho of the current electronic step is computed.
///
/// @param dftu DFT+U object holding the occupation matrices
/// @param ucell unit cell
/// @param global_out_dir output directory (including the trailing separator)
/// @param step electronic-step metadata and convergence values
/// @param conv_esolver whether the electronic SCF is converged at this step
/// @param occmat_ready whether the occupation matrix of this step exists;
///        false records an "N/A" placeholder instead of the matrix body
/// @param cfg frequency configuration
/// @param soc_layout storage layout of the nspin == 4 occupation matrix
void append_ion_step_snapshot(const Plus_U_Base& dftu,
                              const UnitCell& ucell,
                              const std::string& global_out_dir,
                              const OccmatStepInfo& step,
                              bool conv_esolver,
                              bool occmat_ready,
                              const OccmatOutputCfg& cfg,
                              OccmatSocLayout soc_layout);

/// Read the local occupation number matrix from file (rank 0 only).
///
/// Only the current compact format written to occ_mat.txt is accepted. Legacy
/// Atom=/L=/spin= layouts and malformed or incomplete matrices terminate the
/// run with a diagnostic containing the filename and input line.
///
/// @param soc_layout storage layout of the nspin == 4 occupation matrix;
///        the PW path stores 4 contiguous Pauli blocks (b0, b1, b2, b3) and
///        must reconstruct them from the (n_uu, Re(n_ud), Im(n_ud), n_dd)
///        blocks written to the file; the LCAO path stores a real 2m x 2m
///        spin-basis matrix and discards Im(n_ud).
void read_occup_m(const UnitCell& ucell,
                  OccupationMatrix& occ,
                  const std::vector<int>& l_channel,
                  const std::string& fn,
                  int nspin,
                  int npol,
                  OccmatSocLayout soc_layout);

/// Broadcast the local occupation number matrices from rank 0 to all ranks.
///
/// Implemented in dftu_base_io.cpp (only available in MPI builds).
void local_occup_bcast(const UnitCell& ucell,
                       OccupationMatrix& occ,
                       const std::vector<int>& l_channel);

/// Create (or truncate) the per-ionic-step file occ_matg{istep+1}.txt and
/// write its header (rank 0 only).
///
/// At iter == 1 of an enabled output ionic step, truncates the file and writes
/// its provenance header. Other calls are no-ops. Truncating here guarantees
/// that a rerun can not append snapshots from a previous calculation.
void prepare_ion_step_output(const std::string& global_out_dir,
                             int istep,
                             int iter,
                             const OccmatOutputCfg& cfg);

/// Overwrite occ_mat.txt with the occupation matrix of the current
/// electronic step (rank 0 only).
///
/// occ_mat.txt is a single-section snapshot file: the same provenance
/// header and compact layout as occ_matg{#}.txt, whose section header
/// additionally carries the configured scf_thr and the actual drho of
/// this step. Must be called at the iter_finish stage, after drho is
/// computed; it is the entry file of init_chg=file and NSCF restarts.
/// No-op when cfg.out_occ_mat is false.
///
/// @param dftu DFT+U object holding the occupation matrices
/// @param ucell unit cell
/// @param global_out_dir output directory (including the trailing separator)
/// @param step electronic-step metadata and convergence values
/// @param cfg output configuration; only cfg.out_occ_mat is consulted here
/// @param soc_layout storage layout of the nspin == 4 occupation matrix
void write_latest_occmat(const Plus_U_Base& dftu,
                         const UnitCell& ucell,
                         const std::string& global_out_dir,
                         const OccmatStepInfo& step,
                         const OccmatOutputCfg& cfg,
                         OccmatSocLayout soc_layout);

/// Write local occupation matrices to the given stream.
///
/// Writes the compact snapshot layout accepted by read_occup_m(). The output
/// includes per-atom magnetism, electron counts and the full matrix blocks.
/// soc_layout tells how the nspin == 4 storage is arranged (PW Pauli blocks
/// or LCAO real spin-basis matrix).
/// Caller is responsible for opening/closing the stream.
void write_occup_m(const Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   std::ostream& ofs,
                   OccmatSocLayout soc_layout);

} // namespace DFTU_BASE

#endif
