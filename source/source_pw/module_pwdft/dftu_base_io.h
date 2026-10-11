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

/// nested occupation-matrix type used by DFT+U: occ_mat[iat][l][spin](m0, m1)
using OccMatData = std::vector<std::vector<std::vector<ModuleBase::matrix>>>;

/// Frequency configuration for numbered occupation-matrix output.
///
/// Each output ionic step (out_freq_ion) owns one file occ_matg{#}.txt.
/// Within that file, one section is appended at every electronic step
/// selected by out_freq_elec, scf_nmax or convergence.
struct OccmatOutputCfg
{
    int out_freq_ion;  ///< ionic-step interval; 0 disables numbered files
    int out_freq_elec; ///< electronic-iteration interval recorded inside a numbered file
    int scf_nmax;      ///< maximum number of electronic iterations
    bool out_occ_mat;  ///< master switch of occupation-matrix output (out_occ_mat INPUT parameter)
    /// DFT+U calculation mode (INPUT parameter dft_plus_u). The occupation-
    /// matrix IO is meaningful only when dft_plus_u > 0; otherwise the
    /// occupation matrix is never computed and the writers must early-out
    /// to avoid touching an empty l_channel vector. This honours the
    /// documented contract that out_occ_mat only takes effect for DFT+U.
    int dft_plus_u;
};

/// Text format used by write_occup_m().
enum OccmatTextFormat
{
    OCMAT_FMT_LEGACY,   ///< whitespace-token layout parsed by read_occup_m()
    OCMAT_FMT_READABLE  ///< compact human-readable layout for snapshot files
};

/// Storage layout of the real (2l+1)*npol matrix when nspin == 4.
enum OccmatSocLayout
{
    /// 4 contiguous Pauli-component blocks [charge, sigma_x, sigma_y, sigma_z],
    /// each (2l+1)x(2l+1); used by the PW path
    SOC_LAYOUT_PAULI,
    /// real symmetric matrix in spin basis (m + spin-polarization index);
    /// imaginary parts are not stored; used by the LCAO path
    SOC_LAYOUT_SPIN_BASIS_REAL
};

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
/// istep 0 / iter 1 unless it was loaded from file). The caller is
/// responsible for truncating the file at the first electronic step of the
/// ionic step (see prepare_ion_step_file()) and for invoking this function
/// only after drho of the current electronic step is computed.
///
/// @param dftu DFT+U object holding the occupation matrices
/// @param ucell unit cell
/// @param global_out_dir output directory (including the trailing separator)
/// @param nspin number of spin components (1, 2 or 4)
/// @param npol number of polarizations
/// @param istep ionic-step index, starting from 0
/// @param iter electronic-iteration index, starting from 1
/// @param conv_esolver whether the electronic SCF is converged at this step
/// @param occmat_ready whether the occupation matrix of this step exists;
///        false records an "N/A" placeholder instead of the matrix body
/// @param scf_thr configured charge-density convergence threshold
/// @param drho actual charge-density residual of the current electronic step
/// @param cfg frequency configuration
/// @param soc_layout storage layout of the nspin == 4 occupation matrix
void append_ion_step_snapshot(const Plus_U_Base& dftu,
                              const UnitCell& ucell,
                              const std::string& global_out_dir,
                              int nspin,
                              int npol,
                              int istep,
                              int iter,
                              bool conv_esolver,
                              bool occmat_ready,
                              double scf_thr,
                              double drho,
                              const OccmatOutputCfg& cfg,
                              OccmatSocLayout soc_layout);

/// Read the local occupation number matrix from file (rank 0 only).
///
/// The file format matches the output of write_occup_m(). When the file can
/// not be opened, the run quits with an error message that depends on
/// init_occ_mat and init_chg.
///
/// @param soc_layout storage layout of the nspin == 4 occupation matrix;
///        the PW path stores 4 contiguous Pauli blocks (b0, b1, b2, b3) and
///        must reconstruct them from the (n_uu, Re(n_ud), Im(n_ud), n_dd)
///        blocks written to the file; the LCAO path stores a real 2m x 2m
///        spin-basis matrix and discards Im(n_ud).
void read_occup_m(const UnitCell& ucell,
                  OccupationMatrix& occ,
                  const std::vector<int>& l_channel,
                  const int init_occ_mat,
                  const std::string& fn,
                  const std::string& init_chg,
                  int nspin,
                  int npol,
                  OccmatSocLayout soc_layout);

/// Broadcast the local occupation number matrices from rank 0 to all ranks.
///
/// Implemented in dftu_base_io.cpp (only available in MPI builds).
void local_occup_bcast(const UnitCell& ucell,
                       OccupationMatrix& occ,
                       const std::vector<int>& l_channel,
                       int nspin,
                       int npol);

/// Create (or truncate) the per-ionic-step file occ_matg{istep+1}.txt and
/// write its header (rank 0 only).
///
/// Must be called once at the first electronic step (iter == 1) of an output
/// ionic step, before any append_ion_step_snapshot() section. Truncating here
/// also guarantees that a rerun in the same output directory can not append
/// snapshots of a previous calculation.
void prepare_ion_step_file(const std::string& global_out_dir,
                           const int istep,
                           const OccmatOutputCfg& cfg);

/// Output DFT+U information (Hubbard U/J, local occupation matrices) to the
/// running log.
///
/// When cfg.out_occ_mat is true, cfg.out_freq_ion is positive and istep is an
/// output ionic step, the per-ionic-step file occ_matg{istep+1}.txt is
/// created: at the first electronic step (iter == 1) it is truncated and
/// initialized with a provenance header, then append_ion_step_snapshot()
/// appends one section per recorded electronic step.
///
/// Note: occ_mat.txt is NOT written here. It records the actual
/// charge-density residual drho, which is only known after the electronic
/// solve, so write_latest_occmat() writes it from the iter_finish stage.
///
/// Extracted from Plus_U_Base::output as a free function so that IO logic is
/// decoupled from the Plus_U_Base class. The function only reads the
/// Plus_U_Base state via public accessors; no friend declaration needed.
void output(const Plus_U_Base& dftu,
            const UnitCell& ucell,
            const std::string& global_out_dir,
            int nspin,
            int npol,
            int istep,
            int iter,
            const OccmatOutputCfg& cfg,
            OccmatSocLayout soc_layout);

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
/// @param nspin number of spin components (1, 2 or 4)
/// @param npol number of polarizations
/// @param istep ionic-step index, starting from 0
/// @param iter electronic-iteration index, starting from 1
/// @param scf_thr configured charge-density convergence threshold
/// @param drho actual charge-density residual of the current electronic step
/// @param cfg output configuration; only cfg.out_occ_mat is consulted here
/// @param soc_layout storage layout of the nspin == 4 occupation matrix
void write_latest_occmat(const Plus_U_Base& dftu,
                         const UnitCell& ucell,
                         const std::string& global_out_dir,
                         int nspin,
                         int npol,
                         int istep,
                         int iter,
                         double scf_thr,
                         double drho,
                         const OccmatOutputCfg& cfg,
                         OccmatSocLayout soc_layout);

/// Write local occupation matrices to the given stream.
///
/// When diag is true, eigenvalues and magnetism are also printed; otherwise
/// only raw matrix elements. fmt selects between the legacy token layout
/// (used for running log) and the readable snapshot layout of
/// occ_mat.txt and occ_matg{#}.txt. soc_layout tells how the nspin == 4
/// storage is arranged (PW Pauli blocks or LCAO real spin-basis matrix).
/// Caller is responsible for opening/closing the stream.
void write_occup_m(const Plus_U_Base& dftu,
                   const UnitCell& ucell,
                   std::ostream& ofs,
                   bool diag,
                   int nspin,
                   int npol,
                   OccmatTextFormat fmt,
                   OccmatSocLayout soc_layout);

} // namespace DFTU_BASE

#endif
