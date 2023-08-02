#ifndef DEEPKSLCAO_H
#define DEEPKSLCAO_H
#include "module_base/timer.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __DEEPKSTEMPLATE
#define __DEEPKSTEMPLATE

template <class T>
class DeePKS : public T
{
};

#endif

template <typename T>
class DeePKS<OperatorLCAO<T>> : public OperatorLCAO<T>
{
  public:
    DeePKS<OperatorLCAO<T>>(Local_Orbital_Charge* loc_in,
                            LCAO_Matrix* LM_in,
                            const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                            std::vector<double>* HR_pointer_in,
                            std::vector<T>* HK_pointer_in,
                            const int& nks_in)
        : loc(loc_in),
          HR_pointer(HR_pointer_in),
          HK_pointer(HK_pointer_in),
          nks(nks_in),
          OperatorLCAO<T>(LM_in, kvec_d_in)
    {
        this->cal_type = lcao_deepks;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:
    Local_Orbital_Charge* loc;

    std::vector<double>* HR_pointer = nullptr;

    std::vector<T>* HK_pointer = nullptr;

    bool HR_fixed_done = false;

    const int& nks;
};

} // namespace hamilt
#endif