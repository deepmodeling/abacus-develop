#ifndef POTBASE_H
#define POTBASE_H

#include "module_base/complexmatrix.h"
#include "module_base/matrix.h"
#include "module_cell/unitcell.h"
#include "module_basis/module_pw/pw_basis.h"
#include "module_elecstate/module_charge/charge.h"

namespace elecstate
{
/** This class is the base class of Potential module
 1. Main class Potential is derived from it.
 2. components of potentials on real space grids can derived from it and will be registered into Potential.
    a. cal_fixed_v() is a virtual function, it can be override to contribute potentials which do not change with Charge object.
    b. cal_v_eff() is a virtual function, it can be override to contribute potentials which change with Charge object.
    c. fixed_mode should be set "true" if you want Potential class call cal_fixed_v()
    d. dynamic_mode should be set "true" if you want Potential class call cal_v_eff()
    e. rho_basis_ is needed to provide number of real space grids(nrxx) and number of spin(nspin) and FFT(real<->recip) interface
*/
class PotBase
{
  public:
    PotBase(){};
    virtual ~PotBase(){};

    virtual void cal_v_eff(const Charge* chg, const UnitCell* ucell, ModuleBase::matrix& v_eff)
    {
        return;
    }

    virtual void cal_fixed_v(double* vl_pseudo)
    {
        return;
    }

    bool fixed_mode = 0;
    bool dynamic_mode = 0;

  protected:
    const ModulePW::PW_Basis* rho_basis_ = nullptr;
    const ModulePW::PW_Basis* rho_basis_smooth_ = nullptr;
};

} // namespace elecstate

#endif