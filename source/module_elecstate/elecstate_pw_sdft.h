#ifndef ELECSTATEPW_SDFT_H
#define ELECSTATEPW_SDFT_H
#include "elecstate_pw.h"
namespace elecstate
{
    class ElecStatePW_SDFT : public ElecStatePW<std::complex<double>>
    {
    public:
      ElecStatePW_SDFT(ModulePW::PW_Basis_K* wfc_basis_in,
                       Charge* chg_in,
                       K_Vectors* pkv_in,
                       UnitCell* ucell_in,
                       pseudopot_cell_vnl* ppcell_in,
                       ModulePW::PW_Basis* rhodpw_in,
                       ModulePW::PW_Basis* rhopw_in,
                       ModulePW::PW_Basis_Big* bigpw_in)
          : ElecStatePW(wfc_basis_in, chg_in, pkv_in, ucell_in, ppcell_in, rhodpw_in, rhopw_in, bigpw_in)
      {
          this->classname = "ElecStatePW_SDFT";
      }
        virtual void psiToRho(const psi::Psi<std::complex<double>>& psi) override; 
    };
}
#endif