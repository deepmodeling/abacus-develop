#ifndef DHS_SPARSE_WRITER_H
#define DHS_SPARSE_WRITER_H

#include "source_lcao/lcao_hs_arrays.h"

#include <string>

class Parallel_Orbitals;

namespace ModuleIO
{
void save_dH_sparse(const int& istep,
                    const Parallel_Orbitals& pv,
                    LCAO_HS_Arrays& HS_Arrays,
                    const double& sparse_thr,
                    const bool& binary,
                    const std::string& fileflag,
                    const int precision,
                    const std::string& global_out_dir,
                    const std::string& global_matrix_dir,
                    const std::string& calculation,
                    const bool out_app_flag,
                    const int nspin,
                    const int nlocal);
} // namespace ModuleIO

#endif
