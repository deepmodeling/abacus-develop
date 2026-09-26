#ifndef WRITE_HS_SPARSE_H
#define WRITE_HS_SPARSE_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/lcao_hs_arrays.h"

#include <cstddef>
#include <map>
#include <set>
#include <string>

namespace ModuleIO
{
using RCoordinate = Abfs::Vector3_Order<int>;

template <typename T>
using SparseRBlock = std::map<size_t, std::map<size_t, T>>;

template <typename T>
using SparseRMatrix = std::map<RCoordinate, SparseRBlock<T>>;

struct SparseWriteOptions
{
    std::string filename;
    std::string label;
    double threshold = 0.0;
    bool binary = false;
    int precision = 16;
    int istep = -1;
    bool reduce = true;
    std::string temp_dir;
    // Runtime flags that decide file-open mode (append on md restart).
    // Must be provided explicitly by the caller instead of reading PARAM.
    std::string calculation;
    bool out_app_flag = false;
};

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

template <typename Tdata>
void save_sparse(const SparseRMatrix<Tdata>& smat,
                 const std::set<RCoordinate>& all_R_coor,
                 const Parallel_Orbitals& pv,
                 const SparseWriteOptions& options);
} // namespace ModuleIO

#endif
