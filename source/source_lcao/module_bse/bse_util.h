#pragma once

#include <ATen/core/tensor.h>
#include "source_psi/psi.h"
#include "source_base/tool_title.h"
#include "source_base/global_variable.h"
#include "source_base/module_external/blas_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_2d.h"
#include "source_basis/module_ao/parallel_orbitals.h"
#include <RI/global/Tensor.h>
#include <array>
#include <map>
namespace BSE_Util
{
/// ================ info ==================
inline void print_mem_estimate(const std::string& name,
                               const size_t& local_size,
                               const size_t& type_size)
{
    double mem_MB = static_cast<double>(local_size * type_size) / (1024.0 * 1024.0);
    double mem_GB = mem_MB / 1024.0;

    GlobalV::ofs_running << "Allocating " << name << ", memory size: ";
    if (mem_GB >= 1.0) {
        GlobalV::ofs_running << std::fixed << std::setprecision(5) << mem_GB << " GB";
    } else {
        GlobalV::ofs_running << std::fixed << std::setprecision(5) << mem_MB << " MB";
    }
    GlobalV::ofs_running << std::endl;
}

/// ================ RI ==================
using TA = int;
using TC = std::array<int, 3>;
using TAC = std::pair<int, TC>;
template <typename T>
using TLRI = std::map<TA, std::map<TAC, RI::Tensor<T>>>;
template <typename T>
bool move_R_tensor(TLRI<T>& tensor_map, const TA ia, const TA ja, const TC& R_original, const TC& R_new)
{
    auto it_a = tensor_map.find(ia);
    if (it_a == tensor_map.end()) return false;

    auto& map_ja = it_a->second;

    const TAC old_key{ja, R_original};
    auto it_b = map_ja.find(old_key);
    if (it_b == map_ja.end()) return false;

    const TAC new_key{ja, R_new};
    assert(map_ja.count(new_key) == 0 && "Error: target R tensor already exists in move_R_tensor!");

    map_ja[new_key] = std::move(it_b->second);
    map_ja.erase(it_b);
    return true;
}

/// ================ Container ===============
using DAT = container::DataType;
using DEV = container::DeviceType;

/// @brief Struct to get DataType enum for different tensor data types, not used currently
template <typename T>
using DAT_ENUM = container::DataTypeToEnum<T>;

/// =============== Algorithm ===================

/// ================ DM_onebase ===================
#ifdef __MPI
/// @brief calculate the 2d-block transition density matrix in AO basis
/// \f[ \tilde{\rho}_{\mu\mu}=c_{j,\mu}c^*_{b,\nu} \f]
template<typename T>
container::Tensor cal_dm_trans_onebase_pblas(
    const psi::Psi<T>& c,
    const Parallel_2D& pc,
    const int& ik,
    const int& naos,
    const int& imo1, // imo1 → imo2, for excitation imo1[1,nocc], imo2[nocc+1,nocc+nvirt]
    const int& imo2,
    const Parallel_Orbitals& pmat,
    const T& factor = (T)1.0);

#endif
    
/// @brief calculate the transition density matrix in AO basis
template<typename T>
container::Tensor cal_dm_trans_onebase_blas(
    const psi::Psi<T>& c,
    const int& ik,
    const int& naos,
    const int& imo1, // imo1 → imo2, for excitation imo1[1,nocc], imo2[nocc+1,nocc+nvirt]
    const int& imo2,
    const T& factor = (T)1.0);

}
