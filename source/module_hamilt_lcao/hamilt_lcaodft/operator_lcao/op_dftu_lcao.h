#ifndef OPDFTULCAO_H
#define OPDFTULCAO_H
#include "module_base/timer.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __OPDFTUTEMPLATE
#define __OPDFTUTEMPLATE

template <class T>
class OperatorDFTU : public T
{
};

#endif

template <typename T>
class OperatorDFTU<OperatorLCAO<T>> : public OperatorLCAO<T>
{
  public:
    OperatorDFTU<OperatorLCAO<T>>(LCAO_Matrix* LM_in,
                                  const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                  std::vector<double>* HR_pointer_in,
                                  std::vector<T>* HK_pointer_in,
                                  const std::vector<int>& isk_in)
        : HR_pointer(HR_pointer_in), HK_pointer(HK_pointer_in), isk(isk_in), OperatorLCAO<T>(LM_in, kvec_d_in)
    {
        this->cal_type = lcao_dftu;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:
    std::vector<double>* HR_pointer = nullptr;

    std::vector<T>* HK_pointer = nullptr;

    bool HR_fixed_done = false;

    const std::vector<int>& isk;
};
} // namespace hamilt
#endif