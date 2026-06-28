#ifndef WRITE_HS_SPARSE_H
#define WRITE_HS_SPARSE_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/LCAO_HS_arrays.hpp"

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
};

void save_dH_sparse(const int& istep,
                    const Parallel_Orbitals& pv,
                    LCAO_HS_Arrays& HS_Arrays,
                    const double& sparse_thr,
                    const bool& binary,
                    const std::string& fileflag = "h",
                    const int precision = 16);

template <typename Tdata>
void save_sparse(const SparseRMatrix<Tdata>& smat,
                 const std::set<RCoordinate>& all_R_coor,
                 const Parallel_Orbitals& pv,
                 const SparseWriteOptions& options);
} // namespace ModuleIO

#endif
