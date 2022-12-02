#ifndef ELECSTATEPW_SDFT_H
#define ELECSTATEPW_SDFT_H
#include "elecstate_pw.h"
namespace elecstate
{
    class ElecStatePW_SDFT : public ElecStatePW<double>
    {
    public:
        ElecStatePW_SDFT(
            ModulePW::PW_Basis_K *wfc_basis_in, 
            Charge* chg_in, 
            K_Vectors *pkv_in) : ElecStatePW(wfc_basis_in, chg_in, pkv_in)
        {
            this->classname = "ElecStatePW_SDFT";
        }
        virtual void psiToRho(const psi::Psi<std::complex<double>>& psi) override; 
    };
}
#endif