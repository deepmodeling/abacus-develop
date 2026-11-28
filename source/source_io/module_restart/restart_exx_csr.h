#ifndef ABACUS_SOURCE_IO_MODULE_RESTART_RESTART_EXX_CSR_H
#define ABACUS_SOURCE_IO_MODULE_RESTART_RESTART_EXX_CSR_H

#include "source_cell/unitcell.h"
#include "source_lcao/module_ri/ri_util.h"
#include "source_lcao/module_ri/abfs_vector3_order.h"

#include <RI/global/Tensor.h>
#include <map>

namespace ModuleIO
{
using TC = std::array<int, 3>;
using TAC = std::pair<int, TC>;

/// Build the EXX restart CSR file name: hrs{ispin+1}[g{istep+1}]_exx.csr
std::string exx_csr_filename(const int ispin, const bool numbered, const int istep);

/// read Hexxs in CSR format; also extracts E_exx from the description line
template <typename Tdata>
void read_Hexxs_csr(const std::string& file_dir,
                    const UnitCell& ucell,
                    const int nspin,
                    const int nbasis,
                    std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs,
                    double& eexx);

/// read Hexxs in cereal format
template <typename Tdata>
void read_Hexxs_cereal(const std::string& file_name,
                       std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs);

/// write Hexxs in CSR format; embeds E_exx in the description line
template <typename Tdata>
void write_Hexxs_csr(const std::string& out_dir,
                     const UnitCell& ucell,
                     const std::vector<std::map<int, std::map<TAC, RI::Tensor<Tdata>>>>& Hexxs,
                     const double eexx,
                     const int precision,
                     const int istep,
                     const int out_app_flag);

/// calculate CSR sparse matrix from the global matrix stored with RI::Tensor
/// the return type is same as SR_sparse,  HR_sparse, etc.
template <typename Tdata>
std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, Tdata>>> calculate_RI_Tensor_sparse(
    const double& sparse_threshold,
    const std::map<int, std::map<TAC, RI::Tensor<Tdata>>>& Hexxs,
    const UnitCell& ucell);
} // namespace ModuleIO

#endif // ABACUS_SOURCE_IO_MODULE_RESTART_RESTART_EXX_CSR_H
