#ifndef HSR_WRITER_H
#define HSR_WRITER_H

#include "source_base/parallel_2d.h"

#include <complex>
#include <string>
#include <vector>

class UnitCell;

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
template <typename TR>
void write_hcontainer_csr(const std::string& fname,
                          const UnitCell* ucell,
                          const int precision,
                          hamilt::HContainer<TR>* mat_serial,
                          const int istep,
                          const int ispin,
                          const int nspin,
                          const std::string& label,
                          const std::string& representation_note);

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
               const std::string& global_out_dir);

} // namespace ModuleIO

#endif // HSR_WRITER_H
