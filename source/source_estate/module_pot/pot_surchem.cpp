#include "pot_surchem.h"

namespace elecstate
{

PotSurChem::PotSurChem(const ModulePW::PW_Basis* rho_basis_in,
                       Structure_Factor* structure_factors_in,
                       const double* vlocal_in,
                       surchem* surchem_in)
    : vlocal(vlocal_in), surchem_(surchem_in)
{
    this->rho_basis_ = rho_basis_in;
    this->structure_factors_ = structure_factors_in;
    this->dynamic_mode = true;
    this->fixed_mode = false;
}

PotSurChem::~PotSurChem()
{
    if (this->allocated)
    {
        this->surchem_->clear();
    }
}

void PotSurChem::cal_v_eff(const Charge* const chg, const UnitCell* const ucell, ModuleBase::matrix& v_eff)
{
    if (!this->allocated)
    {
        this->surchem_->allocate(this->rho_basis_->nrxx, v_eff.nr);
        this->allocated = true;
    }
    ModuleBase::matrix v_sol_correction(v_eff.nr, this->rho_basis_->nrxx);
    this->surchem_->v_correction(*ucell,
                                 *chg->pgrid,
                                 const_cast<ModulePW::PW_Basis*>(this->rho_basis_),
                                 v_eff.nr,
                                 chg->rho,
                                 this->vlocal,
                                 this->structure_factors_,
                                 v_sol_correction);
    v_eff += v_sol_correction;
}

} // namespace elecstate
