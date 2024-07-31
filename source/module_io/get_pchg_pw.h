#ifndef GET_PCHG_PW_H
#define GET_PCHG_PW_H

#include "module_base/module_device/device.h"
#include "module_cell/unitcell.h"
#include "module_elecstate/elecstate.h"
#include "module_psi/psi.h"
#include "module_pw/pw_basis.h"
#include "module_pw/pw_basis_k.h"

#include <string>
#include <vector>

namespace pchg
{
template <typename Device>
void get_pchg_pw(const std::vector<int>& bands_to_print,
                 const int nbands,
                 const int nspin,
                 const int nx,
                 const int ny,
                 const int nz,
                 const int nxyz,
                 const int nks,
                 const std::vector<int>& isk,
                 const std::vector<double>& wk,
                 const int pw_big_bz,
                 const int pw_big_nbz,
                 const int ngmc,
                 const UnitCell* ucell,
                 psi::Psi<std::complex<double>>* psi,
                 ModulePW::PW_Basis* pw_rhod,
                 ModulePW::PW_Basis_K* pw_wfc,
                 const Device* ctx,
                 const std::string& global_out_dir,
                 const bool if_separate_k);
} // namespace pchg

#endif // GET_PCHG_PW_H
