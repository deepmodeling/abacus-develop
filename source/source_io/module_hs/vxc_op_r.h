#ifndef VXC_OP_R_H
#define VXC_OP_R_H

#include "source_lcao/module_operator_lcao/veff_lcao.h"
#include "source_lcao/module_dftu/dftu_nao_op_legacy.h"
#include "source_lcao/spar_hsr.h"
#include "source_lcao/module_ri/abfs_vector3_order.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_hamilt/module_xc/exx_info.h"

#include <complex>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#ifdef __EXX
#include "source_lcao/module_operator_lcao/op_exx_lcao.h"
#include "source_lcao/module_ri/ri_2d_comm.h"
#endif

namespace ModuleIO
{

/// @brief Write the Vxc matrix in real space (R), useful for GW calculation
/// including terms: local/semi-local XC, EXX, DFTU
/// @tparam TK K-point data type (double or std::complex<double>)
/// @tparam TR Real-space data type (double or std::complex<double>)
template <typename TK, typename TR>
void write_Vxc_R(const int nspin,
                 const Parallel_Orbitals* pv,
                 const UnitCell& ucell,
                 Structure_Factor& sf,
                 surchem& solvent,
                 const ModulePW::PW_Basis& rho_basis,
                 const ModulePW::PW_Basis& rhod_basis,
                 const ModuleBase::matrix& vloc,
                 const Charge& chg,
                 const K_Vectors& kv,
                 const std::vector<double>& orb_cutoff,
                 Grid_Driver& gd,
                 const std::string& global_out_dir,
                 bool cal_exx,
                 double hybrid_alpha,
                 bool real_number
#ifdef __EXX
                 ,
                 const std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>* Hexxd,
                 const std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>* Hexxc
#endif
                 ,
                 const double sparse_thr,
                 std::ofstream& ofs_running);

} // namespace ModuleIO

#endif // VXC_OP_R_H
