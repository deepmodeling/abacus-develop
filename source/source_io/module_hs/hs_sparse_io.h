#ifndef HS_SPARSE_IO_H
#define HS_SPARSE_IO_H

#include "source_basis/module_ao/parallel_orbitals.h"
#include "source_lcao/module_ri/abfs_vector3_order.h"

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
    // Runtime flags that decide file-open mode (append on md restart).
    // Must be provided explicitly by the caller instead of reading PARAM.
    std::string calculation;
    bool out_app_flag = false;
};

template <typename Tdata>
void save_sparse(const SparseRMatrix<Tdata>& smat,
                 const std::set<RCoordinate>& all_R_coor,
                 const Parallel_Orbitals& pv,
                 const SparseWriteOptions& options);
} // namespace ModuleIO

#endif
