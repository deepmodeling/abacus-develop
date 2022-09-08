#ifndef ELECSTATELCAOTDDFT_H
#define ELECSTATELCAOTDDFT_H

#include "elecstate.h"
#include "elecstate_lcao.h"
#include "src_lcao/LCAO_hamilt.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/local_orbital_wfc.h"

namespace elecstate
{

class ElecStateLCAO_TDDFT : public ElecStateLCAO
{
  public:
    ElecStateLCAO_TDDFT(Charge* chg_in = nullptr,
                        const K_Vectors* klist_in = nullptr,
                        int nks_in = 1,
                        int nbands_in = 1,
                        Local_Orbital_Charge* loc_in = nullptr,
                        LCAO_Hamilt* uhm_in = nullptr,
                        Local_Orbital_wfc* lowf_in = nullptr)
    {
        init(chg_in, klist_in, nks_in, nbands_in);
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