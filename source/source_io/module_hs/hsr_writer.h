#ifndef HSR_WRITER_H
#define HSR_WRITER_H

#include "source_base/parallel_blkcyc.h"

#include <complex>
#include <string>
#include <vector>

class UnitCell;

namespace elecstate
{
struct Efermi;
}

namespace hamilt
{
template <typename T>
class HContainer;
} // namespace hamilt

namespace ModuleIO
{

/// Generate filename for spin-dependent HR output.
/// @param numbered whether the ionic step is embedded as a g{step} suffix (out_app_flag == 2)
std::string hsr_gen_fname(const std::string& prefix,
                          const int ispin,
                          const bool numbered,
                          const int istep);

/// Generate filename for spin-dependent HR output in the selected format.
std::string hsr_gen_fname(const std::string& prefix,
                          const int ispin,
                          const bool numbered,
                          const int istep,
                          const int out_type);

/// Generate filename for spin-independent SR output.
std::string sr_gen_fname(const bool numbered, const int istep);

/// Generate filename for spin-independent SR output in the selected format.
std::string sr_gen_fname(const bool numbered, const int istep, const int out_type);

/// Generate filename for derivative matrices (dH/dR, dS/dR).
std::string dhr_gen_fname(const std::string& prefix,
                          const int ispin,
                          const bool numbered,
                          const int istep);

/// Generate filename for the k-space derivative matrices (dH(k), dT(k), ...).
/// A spin suffix is added only for nspin == 2 so the two spin channels do not
/// overwrite each other; a single spin channel keeps the bare name.
/// @param numbered whether the ionic step is embedded as a g{step} suffix (out_app_flag == 2)
std::string dhk_gen_fname(const std::string& prefix,
                          const int ispin,
                          const int nspin,
                          const bool numbered,
                          const int istep);

/// Write a single HContainer to CSR file with header.
/// @param append whether to append to an existing file; the first written step (istep <= 0) truncates
/// @param efermi_eV the Fermi energy in eV for this spin channel
/// @param has_efermi whether to append the Fermi energy to the spin-index line
/// @return true only if the whole record was written and the stream closed
///         cleanly; callers must not announce the file in the running log
///         when this returns false.
template <typename TR>
bool write_hcontainer_csr(const std::string& fname,
                          const UnitCell* ucell,
                          const int precision,
                          hamilt::HContainer<TR>* mat_serial,
                          const int istep,
                          const bool append,
                          const int ispin,
                          const int nspin,
                          const std::string& label,
                          const std::string& representation_note,
                          const double efermi_eV,
                          const bool has_efermi);

/// Write one HContainer record in the native binary CSR format.
template <typename TR>
void write_hcontainer_csr_binary(const std::string& fname,
                                 hamilt::HContainer<TR>* mat_serial,
                                 const int istep,
                                 const bool append);

/// Write H(R) and S(R) in CSR format, unified with write_dmr interface.
/// @param out_app_flag output mode: 0 = fixed filename, overwritten; 1 = fixed filename, appended;
///                     2 = per-step g{step} filename.
template <typename TR>
void write_hsr(const std::vector<hamilt::HContainer<TR>*>& hr_vec,
               const hamilt::HContainer<TR>* sr,
               const UnitCell* ucell,
               const int out_type,
               const int precision,
               const Parallel_Blkcyclic& paraV,
               const int out_app_flag,
               const bool gamma_only,
               const int* iat2iwt,
               const int nat,
               const int istep,
               const std::string& global_out_dir,
               const elecstate::Efermi& eferm,
               std::ofstream& ofs_running);

} // namespace ModuleIO

#endif // HSR_WRITER_H
