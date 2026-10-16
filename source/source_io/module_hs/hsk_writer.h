#ifndef HSK_WRITER_H
#define HSK_WRITER_H

#include <fstream>
#include <string>
#include <vector>

#include "source_base/parallel_2d.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_hamilt/hamilt.h"

// mohan note: this file holds the H(k)/S(k) writer
// (formerly write_hsk in write_hs.h / write_hs.hpp).

namespace ModuleIO
{

template <typename T>
void write_hsk(const std::string& global_out_dir,
               const int nspin,
               const int nks,
               const int nkstot,
               const std::vector<int>& ik2iktot,
               const std::vector<int>& isk,
               hamilt::Hamilt<T>* p_hamilt,
               const Parallel_Orbitals& pv,
               const bool gamma_only,
               const bool out_app_flag,
               const int istep,
               const int out_type,
               const int precision,
               const int nlocal,
               const std::string& ks_solver,
               const int drank,
               std::ofstream& ofs_running);

} // namespace ModuleIO

#endif // HSK_WRITER_H
