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
               surchem* surchem_in);
    ~PotSurChem();

    // Passing an explicit output matrix makes the lifetime and allocation explicit and avoids hidden allocations.
    void cal_v_eff(const Charge* const chg, const UnitCell* const ucell, ModuleBase::matrix& v_eff) override;

  private:
    surchem* surchem_ = nullptr;
    Structure_Factor* structure_factors_ = nullptr;
    const double* vlocal = nullptr;
    bool allocated = false;
};

} // namespace elecstate

#endif
