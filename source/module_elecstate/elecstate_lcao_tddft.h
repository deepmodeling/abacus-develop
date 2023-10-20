#ifndef ELECSTATELCAOTDDFT_H
#define ELECSTATELCAOTDDFT_H

#include "elecstate.h"
#include "elecstate_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"

namespace elecstate
{
class ElecStateLCAO_TDDFT : public ElecStateLCAO<std::complex<double>>
{
  public:
    ElecStateLCAO_TDDFT(Charge* chg_in ,
                        const K_Vectors* klist_in ,
                        int nks_in,
                        Local_Orbital_Charge* loc_in ,
                        LCAO_Hamilt* uhm_in ,
                        Local_Orbital_wfc* lowf_in ,
                        ModulePW::PW_Basis* rhopw_in ,
                        ModulePW::PW_Basis_Big* bigpw_in )
    {
        init_ks(chg_in, klist_in, nks_in, rhopw_in, bigpw_in);        
        this->loc = loc_in;
        this->uhm = uhm_in;
        this->lowf = lowf_in;
        this->classname = "ElecStateLCAO_TDDFT";
    }
    void psiToRho_td(const psi::Psi<std::complex<double>>& psi);
    void calculate_weights_td();
};

} // namespace elecstate

#endif