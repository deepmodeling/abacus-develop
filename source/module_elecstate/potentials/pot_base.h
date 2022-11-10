#ifndef POTBASE_H
#define POTBASE_H

#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include "module_pw/pw_basis.h"
#include "src_pw/charge.h"
#include "module_cell/unitcell_pseudo.h"


namespace elecstate
{

class PotBase
{
    public:
    PotBase(){};
    virtual ~PotBase(){};

    virtual void cal_v_eff(
        const Charge* chg, 
        const UnitCell_pseudo* ucell, 
        ModuleBase::matrix& v_eff)
    {return;}

    virtual void cal_fixed_v(double *vl_pseudo){return;}

    bool fixed_mode = 0;
    bool dynamic_mode = 0;
    
    protected:

    const ModulePW::PW_Basis *rho_basis_ = nullptr;
};

}

#endif