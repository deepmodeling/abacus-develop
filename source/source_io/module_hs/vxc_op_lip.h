#ifndef VXC_OP_LIP_H
#define VXC_OP_LIP_H

#include "source_psi/psi.h"
#include "source_cell/unitcell.h"
#include "source_cell/klist.h"
#include "source_estate/module_pot/potential_new.h"
#include "source_pw/module_pwdft/op_pw_veff.h"

#include <complex>
#include <string>
#include <vector>

#ifdef __EXX
#include "source_lcao/module_ri/exx_lip.h"
#endif

namespace ModuleIO
{

/// @brief Write the Vxc matrix in KS orbital representation for LIP (LCAO-in-PW), useful for GW calculation
/// including terms: local/semi-local XC and EXX
/// @tparam FPTYPE Floating point type (float or double)
template <typename FPTYPE>
void write_Vxc_LIP(int nspin,
                   int naos,
                   int drank,
                   const psi::Psi<std::complex<FPTYPE>>& psi_pw,
                   const UnitCell& ucell,
                   Structure_Factor& sf,
                   surchem& solvent,
                   const ModulePW::PW_Basis_K& wfc_basis,
                   const ModulePW::PW_Basis& rho_basis,
                   const ModulePW::PW_Basis& rhod_basis,
                   const ModuleBase::matrix& vloc,
                   const Charge& chg,
                   const K_Vectors& kv,
                   const ModuleBase::matrix& wg,
                   const bool gamma_only,
                   const std::string& global_out_dir,
                   const int out_ndigits,
                   const std::string& ks_solver,
                   bool cal_exx,
                   double hybrid_alpha
#ifdef __EXX
                   ,
                   const Exx_Lip<std::complex<FPTYPE>>& exx_lip
#endif
);

} // namespace ModuleIO

#endif // VXC_OP_LIP_H
