#ifndef POTSURCHEM_H
#define POTSURCHEM_H

#include "pot_base.h"
#include "module_surchem/surchem.h"

namespace elecstate
{

class PotSurChem : public PotBase
{
    public:
    //constructor for exchange-correlation potential
    //meta-GGA should input matrix of kinetic potential, it is optional
    PotSurChem(
        const ModulePW::PW_Basis* rho_basis_in,
        surchem* surchem_in
    ):surchem_(surchem_in)
    {
        this->rho_basis_ = rho_basis_in;
        this->dynamic_mode = true;
        this->fixed_mode = false;
    }

    void cal_v_eff(
        const Charge* chg, 
        const UnitCell_pseudo* ucell, 
        ModuleBase::matrix& v_eff)override
    {
            surchem_->v_correction(
                *ucell, 
                const_cast<ModulePW::PW_Basis *>(this->rho_basis_), 
                v_eff.nr, 
                chg->rho);
    }

    private:
    surchem* surchem_;
};

}

#endif