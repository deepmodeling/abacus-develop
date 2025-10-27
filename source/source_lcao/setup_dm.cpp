#include "source_lcao/setup_dm.h"

#include "cal_dm.h"
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

template <typename TK>
void Setup_DM<TK>::init_DM(const K_Vectors* kv, const Parallel_Orbitals* paraV, const int nspin)
{
    const int nspin_dm = nspin == 2 ? 2 : 1;
    this->dm = new DensityMatrix<TK, double>(paraV, nspin_dm, kv->kvec_d, kv->get_nks() / nspin_dm);
}

template class Setup_DM<double>;               // Gamma_only case
template class Setup_DM<std::complex<double>>; // multi-k case

} // namespace elecstate
