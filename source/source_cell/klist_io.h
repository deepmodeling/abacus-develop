#ifndef KLIST_IO_H
#define KLIST_IO_H

#include "source_base/matrix3.h"
#include "source_base/vector3.h"

#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace ModuleSymmetry
{
class Symmetry; // full definition only needed in klist_io.cpp
}

class UnitCell; // full definition only needed in klist_io.cpp

/// this-free helpers extracted from K_Vectors, kept in a separate TU so they
/// can be unit-tested and reused without dragging in the K_Vectors class.
namespace KListIO
{
/// Render the IBZ reduction table ("IBZ" k-point -> originating k-point).
std::string ibz_kpt_table(int nkstot,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                          const std::vector<int>& ibz_index,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz);

/// Render the IBZ weight table (IBZ k-point, weight, multiplicity, origin index).
std::string ibz_wk_table(int nkstot_ibz,
                         const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz,
                         const std::vector<double>& wk_ibz,
                         const std::vector<int>& ibz2bz);

/// Result of line-mode interpolation between special k-points.
struct LineK
{
    std::vector<ModuleBase::Vector3<double>> kpts; ///< interpolated k points
    std::vector<int> segids;                       ///< segment id per k point (ISSUE#3482)
    int nks_total = 0;                             ///< total interpolated k-point count
};

/// Scan `ifk` for the "K_POINTS"/"KPOINTS"/"K" header keyword, skipping any
/// leading comment lines. Returns true with the stream positioned after the
/// header line; returns false if the keyword is not found before EOF.
bool find_kpoints_header(std::ifstream& ifk);

/// Read `nkstot` explicit k points (three coordinates plus a weight per line)
/// from `ifk` into `kvec` and `wk`. The caller is responsible for sizing the
/// arrays (K_Vectors::renew) before calling.
void read_kpt_list(std::ifstream& ifk,
                   int nkstot,
                   std::vector<ModuleBase::Vector3<double>>& kvec,
                   std::vector<double>& wk);

/// Read the special k points and per-point interpolation counts from `ifk`,
/// then linearly interpolate the line-mode k points. Pure function of the
/// stream and `nks_special`; dies via WARNING_QUIT on malformed input.
LineK interp_line(std::ifstream& ifk, int nks_special);

/// Build the EXX k-stars: for every k point, find the symmetry operation
/// (index into `kgmatrix`) that rotates it onto an irreducible k point, and
/// group k points by that IBZ representative. `equal` compares two doubles
/// with the symmetry precision; `epsilon` is the k-restriction tolerance.
/// this-free so the heavy triple loop is isolated and testable.
void build_kstars(const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                  const std::vector<ModuleBase::Matrix3>& kgmatrix,
                  int nrotkm,
                  const std::vector<ModuleBase::Vector3<double>>& kvec_d_ibz,
                  double epsilon,
                  const std::function<bool(double, double)>& equal,
                  std::vector<std::map<int, ModuleBase::Vector3<double>>>& kstars);

/// Append the time-reversal-related k-point symmetry operations into
/// `kgmatrix` (the slots right after the first `nrotkm` operations must be
/// available). For magnetic nspin=4 systems the antiunitary Theta*g coset
/// is appended from `symm.kgmatrix_anti`; otherwise the inverted -g ops are
/// appended unless inversion is already present. Returns the updated total
/// operation count.
int append_time_reversal_ops(const ModuleSymmetry::Symmetry& symm,
                             std::vector<ModuleBase::Matrix3>& kgmatrix,
                             int nrotkm);

/// Flatten k-point arrays into contiguous MPI buffers (x,y,z interleaved).
/// this-free; used on rank 0 before broadcasting in K_Vectors::mpi_k.
void pack_kpts(const std::vector<int>& isk,
               const std::vector<double>& wk,
               const std::vector<ModuleBase::Vector3<double>>& kvec_c,
               const std::vector<ModuleBase::Vector3<double>>& kvec_d,
               const std::vector<ModuleBase::Vector3<double>>& kvec_c_full,
               int nkstot,
               std::vector<int>& isk_aux,
               std::vector<double>& wk_aux,
               std::vector<double>& kvec_c_aux,
               std::vector<double>& kvec_d_aux,
               std::vector<double>& kvec_c_full_aux);

/// Broadcast the EXX k-stars (one (symmetry-index, k-vector) map per IBZ
/// k-point) from `my_rank == 0` to every process. Rank 0 holds the filled
/// maps; other ranks resize and rebuild them from the broadcast. MPI
/// wrappers are compiled as no-ops without __MPI, so the call is safe in
/// serial builds (it simply leaves the rank-0 maps untouched).
void bcast_kstars(std::vector<std::map<int, ModuleBase::Vector3<double>>>& kstars,
                  int nkstot,
                  int my_rank);

/// Scatter the broadcast buffers into this pool's k-point slice, starting at
/// global index `startk`. this-free; mirrors pack_kpts after the broadcast.
void unpack_kpts(const std::vector<int>& isk_aux,
                 const std::vector<double>& wk_aux,
                 const std::vector<double>& kvec_c_aux,
                 const std::vector<double>& kvec_d_aux,
                 const std::vector<double>& kvec_c_full_aux,
                 int nks,
                 int startk,
                 std::vector<int>& isk,
                 std::vector<double>& wk,
                 std::vector<ModuleBase::Vector3<double>>& kvec_c,
                 std::vector<ModuleBase::Vector3<double>>& kvec_d,
                 std::vector<ModuleBase::Vector3<double>>& kvec_c_full);

/// Fill the full-list Cartesian k vectors when only one coordinate set is
/// available: direct coordinates are converted via `reciprocal_vec` when
/// Cartesian points are missing, otherwise the Cartesian points are copied.
/// No-op when both coordinate sets are done. this-free helper called from
/// K_Vectors::set() before IBZ reduction.
void fill_full_kvec(bool kc_done,
                    bool kd_done,
                    int nkstot_nospin,
                    const ModuleBase::Matrix3& reciprocal_vec,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_c,
                    const std::vector<ModuleBase::Vector3<double>>& kvec_d,
                    std::vector<ModuleBase::Vector3<double>>& kvec_c_full);

/// Build the local-to-global k-point index map `ik2iktot` for this pool.
/// In MPI runs the global index is offset by the pool start (with the
/// spin_mult == 2 second half offset by nkstot/2); in serial runs it is the
/// local index itself. `ik2iktot` is resized to `nks` here. this-free helper
/// called from K_Vectors::set() after the pool distribution.
void build_ik2iktot(int my_pool,
                    const std::vector<int>& startk_pool,
                    int spin_mult,
                    int nks,
                    int nkstot,
                    std::vector<int>& ik2iktot);

/// Expand the k-point list for spin-polarized runs (spin_mult == 2): copy
/// coordinates and weights into the second half, tag isk 0 for the first
/// half and 1 for the second, then double nks/nkstot. For spin_mult == 1
/// only isk is zeroed. The running-log output stays in the K_Vectors
/// wrapper. this-free helper backing K_Vectors::set_kup_and_kdw.
void expand_spin_kpoints(int spin_mult,
                         std::vector<ModuleBase::Vector3<double>>& kvec_c,
                         std::vector<ModuleBase::Vector3<double>>& kvec_d,
                         std::vector<double>& wk,
                         std::vector<int>& isk,
                         int& nks,
                         int& nkstot);

/// Overwrite the KPT file with an auto-generated mesh when requested:
/// a single Gamma point if gamma_only_local, or a KSPACING-derived
/// Gamma/Monkhorst-Pack mesh if kspacing[0] > 0 (quits if kspacing[1] or
/// kspacing[2] is non-positive); does nothing otherwise. this-free helper
/// backing K_Vectors::generate_kfile.
void write_auto_kfile(const UnitCell& ucell,
                      const std::string& fn,
                      bool gamma_only_local,
                      const double kspacing[3],
                      const std::string& kmesh_type,
                      const double koffset[3],
                      std::ofstream& ofs_warning);
} // namespace KListIO

#endif // KLIST_IO_H
