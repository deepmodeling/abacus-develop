#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_ELECSTATE_ELECSTATE_LCAO_TDDFT_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_ELECSTATE_ELECSTATE_LCAO_TDDFT_H

#include "elecstate.h"
#include "elecstate_lcao.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"

namespace elecstate
{
class ElecStateLCAO_TDDFT : public ElecStateLCAO<std::complex<double>>
{
  public:
    ElecStateLCAO_TDDFT(Charge* chg_in,
                        const K_Vectors* klist_in,
                        int nks_in,
                        Gint_k* gint_k_in, // mohan add 2024-04-01
                        ModulePW::PW_Basis* rhopw_in,
                        ModulePW::PW_Basis_Big* bigpw_in)
    {
        init_ks(chg_in, klist_in, nks_in, rhopw_in, bigpw_in);
        this->gint_k = gint_k_in;
        this->classname = "ElecStateLCAO_TDDFT";
    }
    void psiToRho_td(const psi::Psi<std::complex<double>>& psi);
    void calculate_weights_td();
};

} // namespace elecstate

#endif
