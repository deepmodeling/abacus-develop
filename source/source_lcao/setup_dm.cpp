#include "source_lcao/setup_dm.h"

#include "source_estate/cal_dm.h"
#include "source_base/timer.h"
#include "source_estate/module_dm/cal_dm_psi.h"
#include "source_hamilt/module_xc/xc_functional.h"
#include "source_lcao/module_deltaspin/spin_constrain.h"
#include "source_pw/module_pwdft/global.h"
#include "source_io/module_parameter/parameter.h"

#include "source_lcao/module_gint/gint_interface.h"

#include <vector>

namespace LCAO_domain
{

// change init_dm to allocate_dm, mohan 2025-10-31
template <typename TK>
void Setup_DM<TK>::allocate_dm(const K_Vectors* kv, const Parallel_Orbitals* pv, const int nspin)
{
    const int nspin_dm = nspin == 2 ? 2 : 1;
    this->dm = new elecstate::DensityMatrix<TK, double>(pv, nspin_dm, kv->kvec_d, kv->get_nks() / nspin_dm);
}

template class Setup_DM<double>;               // Gamma_only case
template class Setup_DM<std::complex<double>>; // multi-k case

} // namespace elecstate
