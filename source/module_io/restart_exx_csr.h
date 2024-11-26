#pragma once
#include "module_base/abfs-vector3_order.h"
#include "module_cell/unitcell.h"
#include "module_ri/serialization_cereal.h"
#include <RI/global/Tensor.h>
#include <map>

namespace ModuleIO
{
    using TC = std::array<int, 3>;
    using TAC = std::pair<int, TC>;

    /// read Hexxs in CSR format
    template<typename Tdata>
    void read_Hexxs_csr(const std::string& file_name, const UnitCell& ucell,
        const int nspin, const int nbasis,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs);

    /// read Hexxs in cereal format
    template<typename Tdata>
    void read_Hexxs_cereal(const std::string& file_name,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs);

    /// write Hexxs in CSR format
    template<typename Tdata>
    void write_Hexxs_csr(const std::string& file_name, const UnitCell& ucell,
        const std::map<int, std::map<TAC, RI::Tensor<Tdata>>>& Hexxs);

    /// calculate CSR sparse matrix from the global matrix stored with RI::Tensor
    /// the return type is same as SR_sparse,  HR_sparse, etc.
    template<typename Tdata>
    std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, Tdata>>>
        calculate_RI_Tensor_sparse(const double& sparse_threshold,
            const std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs,
            const UnitCell& ucell);
}

#include "module_io/restart_exx_csr.hpp"