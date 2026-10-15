#ifndef VXC_OP_MAT_H
#define VXC_OP_MAT_H

#include "source_psi/psi.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_hamilt/module_xc/exx_info.h"
#include "source_lcao/module_operator_lcao/veff_lcao.h"
#include "source_lcao/module_dftu/dftu_nao_op_legacy.h"
#include "source_io/module_hs/vxc_op_tools.h"

#include <complex>
#include <map>
#include <string>
#include <vector>

#ifdef __EXX
#include "source_lcao/module_operator_lcao/op_exx_lcao.h"
#endif

namespace hamilt
{
template <typename T>
class HContainer;
} // namespace hamilt

namespace ModuleIO
{

/// @brief Write the Vxc matrix in KS orbital representation, useful for GW calculation
/// including terms: local/semi-local XC, EXX, DFTU
/// @tparam TK K-point data type (double or std::complex<double>)
/// @tparam TR Real-space data type (double or std::complex<double>)
template <typename TK, typename TR>
void write_Vxc(const int nspin,
               const int nbasis,
               const int drank,
               const Parallel_Orbitals* pv,
               const psi::Psi<TK>& psi,
               const UnitCell& ucell,
               Structure_Factor& sf,
               surchem& solvent,
               const ModulePW::PW_Basis& rho_basis,
               const ModulePW::PW_Basis& rhod_basis,
               const ModuleBase::matrix& vloc,
               const Charge& chg,
               const K_Vectors& kv,
               const std::vector<double>& orb_cutoff,
               const ModuleBase::matrix& wg,
               Grid_Driver& gd,
               const bool dft_plus_u,
               const bool gamma_only,
               const std::string& global_out_dir,
               const int out_ndigits,
               const std::string& ks_solver,
               bool cal_exx,
               const Exx_Info& exx_info
#ifdef __EXX
               ,
               std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<double>>>>* Hexxd,
               std::vector<std::map<int, std::map<hamilt::TAC, RI::Tensor<std::complex<double>>>>>* Hexxc
#endif
);

} // namespace ModuleIO

#endif // VXC_OP_MAT_H
