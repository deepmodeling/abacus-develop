#ifndef HSR_WRITER_H
#define HSR_WRITER_H

#include "source_base/parallel_2d.h"

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
std::string hsr_gen_fname(const std::string& prefix,
                          const int ispin,
                          const bool append,
                          const int istep);

/// Generate filename for spin-dependent HR output in the selected format.
std::string hsr_gen_fname(const std::string& prefix,
                          const int ispin,
                          const bool append,
                          const int istep,
                          const int out_type);

/// Generate filename for spin-independent SR output.
std::string sr_gen_fname(const bool append, const int istep);

/// Generate filename for spin-independent SR output in the selected format.
std::string sr_gen_fname(const bool append, const int istep, const int out_type);

/// Generate filename for derivative matrices (dH/dR, dS/dR).
std::string dhr_gen_fname(const std::string& prefix,
                          const int ispin,
                          const bool append,
                          const int istep);

/// Write a single HContainer to CSR file with header.
/// @param efermi_eV the Fermi energy in eV for this spin channel
/// @param has_efermi whether to append the Fermi energy to the spin-index line
template <typename TR>
void write_hcontainer_csr(const std::string& fname,
                          const UnitCell* ucell,
                          const int precision,
                          hamilt::HContainer<TR>* mat_serial,
                          const int istep,
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
template <typename TR>
void write_hsr(const std::vector<hamilt::HContainer<TR>*>& hr_vec,
               const hamilt::HContainer<TR>* sr,
               const UnitCell* ucell,
               const int out_type,
               const int precision,
               const Parallel_2D& paraV,
               const bool append,
               const bool gamma_only,
               const int* iat2iwt,
               const int nat,
               const int istep,
               const std::string& global_out_dir,
               const elecstate::Efermi& eferm,
               std::ofstream& ofs_running);

} // namespace ModuleIO

#endif // HSR_WRITER_H
