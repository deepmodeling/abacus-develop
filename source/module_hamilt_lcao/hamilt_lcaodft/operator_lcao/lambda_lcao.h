#ifndef LAMBDA_LCAO_H
#define LAMBDA_LCAO_H

#include "module_base/timer.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __OPLAMBDATEMPLATE
#define __OPLAMBDATEMPLATE

template <class T>
class OperatorLambda : public T
{
};

#endif

template <typename T>
class OperatorLambda<OperatorLCAO<T>> : public OperatorLCAO<T>
{
  public:
    OperatorLambda<OperatorLCAO<T>>(LCAO_Matrix* LM_in,
                                    std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                    std::vector<double>* HR_pointer_in,
                                    std::vector<T>* HK_pointer_in)
        : HR_pointer(HR_pointer_in), HK_pointer(HK_pointer_in), OperatorLCAO<T>(LM_in, kvec_d_in)
    {
        this->cal_type = lcao_constrained_m;
    }

    ~OperatorLambda<OperatorLCAO<T>>();

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:
    std::vector<double>* HR_pointer = nullptr;
    std::vector<T>* HK_pointer = nullptr;
};

} // namespace hamilt

#endif // LAMBDA_LCAO_H