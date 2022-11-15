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
        const double* vlocal_in,
        surchem* surchem_in
    ):vlocal(vlocal_in), surchem_(surchem_in)
    {
        this->rho_basis_ = rho_basis_in;
        this->dynamic_mode = true;
        this->fixed_mode = false;
    }

    void cal_v_eff(
        const Charge* chg, 
        const UnitCell* ucell, 
        ModuleBase::matrix& v_eff)override
    {
            v_eff += surchem_->v_correction(
                *ucell, 
                const_cast<ModulePW::PW_Basis *>(this->rho_basis_), 
                v_eff.nr, 
                chg->rho,
                this->vlocal);
    }

    private:
    surchem* surchem_ = nullptr;
    const double* vlocal = nullptr;
};

}

#endif