#ifndef NONLOCALLCAO_H
#define NONLOCALLCAO_H
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_gen_fixedH.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __NONLOCALTEMPLATE
#define __NONLOCALTEMPLATE

template <class T>
class Nonlocal : public T
{
};

#endif

template <typename T>
class Nonlocal<OperatorLCAO<T>> : public OperatorLCAO<T>
{
  public:
    Nonlocal<OperatorLCAO<T>>(LCAO_gen_fixedH* genH_in,
                              LCAO_Matrix* LM_in,
                              const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                              std::vector<double>* HR_pointer_in,
                              std::vector<T>* HK_pointer_in)
        : genH(genH_in), HR_pointer(HR_pointer_in), HK_pointer(HK_pointer_in), OperatorLCAO<T>(LM_in, kvec_d_in)
    {
        this->cal_type = lcao_fixed;
    }

    virtual void contributeHR() override;

  private:
    // use overlap matrix to generate fixed Hamiltonian
    LCAO_gen_fixedH* genH = nullptr;

    std::vector<double>* HR_pointer = nullptr;

    std::vector<T>* HK_pointer = nullptr;

    bool HR_fixed_done = false;
};

} // namespace hamilt
#endif