#ifndef VEFFLCAO_H
#define VEFFLCAO_H
#include "module_base/timer.h"
#include "module_elecstate/potentials/potential_new.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __VEFFTEMPLATE
#define __VEFFTEMPLATE

template <class T>
class Veff : public T
{
};

#endif

template <typename T>
class Veff<OperatorLCAO<T>> : public OperatorLCAO<T>
{
  public:
    Veff<OperatorLCAO<T>>(Gint_k* GK_in,
                          Local_Orbital_Charge* loc_in,
                          LCAO_Matrix* LM_in,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                          elecstate::Potential* pot_in,
                          std::vector<double>* HR_pointer_in,
                          std::vector<T>* HK_pointer_in)
        : GK(GK_in),
          loc(loc_in),
          pot(pot_in),
          HR_pointer(HR_pointer_in),
          HK_pointer(HK_pointer_in),
          OperatorLCAO<T>(kvec_d_in)
    {
        this->LM = LM_in;
        this->cal_type = lcao_gint;
    }
    Veff<OperatorLCAO<T>>(Gint_Gamma* GG_in,
                          Local_Orbital_Charge* loc_in,
                          LCAO_Matrix* LM_in,
                          elecstate::Potential* pot_in,
                          std::vector<double>* HR_pointer_in,
                          std::vector<T>* HK_pointer_in,
                          std::vector<ModuleBase::Vector3<double>> kvec_d_in)
        : GG(GG_in), loc(loc_in), pot(pot_in), HR_pointer(HR_pointer_in), HK_pointer(HK_pointer_in), OperatorLCAO<T>(std::vector<ModuleBase::Vector3<double>>{ModuleBase::Vector3<double>(0,0,0)})
    {
        this->LM = LM_in;
        this->cal_type = lcao_gint;
    }

    ~Veff<OperatorLCAO<T>>();

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:
    // used for k-dependent grid integration.
    Gint_k* GK = nullptr;

    // used for gamma only algorithms.
    Gint_Gamma* GG = nullptr;

    // Charge calculating method in LCAO base and contained grid base calculation: DM_R, DM, pvpR_reduced
    Local_Orbital_Charge* loc = nullptr;

    std::vector<double>* HR_pointer = nullptr;

    std::vector<T>* HK_pointer = nullptr;

    elecstate::Potential* pot = nullptr;

    bool allocated_pvpR = false;
};

} // namespace hamilt
#endif