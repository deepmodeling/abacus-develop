#ifndef POTSURCHEM_H
#define POTSURCHEM_H

#include "source_hamilt/module_surchem/surchem.h"
#include "pot_base.h"

namespace elecstate
{

class PotSurChem : public PotBase
{
  public:
    // constructor for exchange-correlation potential
    // meta-GGA should input matrix of kinetic potential, it is optional
    PotSurChem(const ModulePW::PW_Basis* rho_basis_in,
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
    ~PotSurChem()
    {
        if (this->allocated)
        {
            this->surchem_->clear();
        }
    }

    // Passing an explicit output matrix makes the lifetime and allocation explicit and avoids hidden allocations.
    void cal_v_eff(const Charge*const chg, const UnitCell*const ucell, ModuleBase::matrix& v_eff) override
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

  private:
    surchem* surchem_ = nullptr;
    Structure_Factor* structure_factors_ = nullptr;
    const double* vlocal = nullptr;
    bool allocated = false;
};

} // namespace elecstate

#endif