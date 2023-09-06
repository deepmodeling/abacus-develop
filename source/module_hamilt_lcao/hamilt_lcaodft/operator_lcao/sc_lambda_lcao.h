#ifndef SC_LAMBDA_LCAO_H
#define SC_LAMBDA_LCAO_H

#include "module_base/timer.h"
#include "module_base/vector3.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __OPLAMBDATEMPLATE
#define __OPLAMBDATEMPLATE

template <class T>
class OperatorScLambda : public T
{
};

#endif

template<typename T>
class OperatorScLambda<OperatorLCAO<T>> : public OperatorLCAO<T>
{
  public:
    OperatorScLambda<OperatorLCAO<T>>(LCAO_Matrix* LM_in,
                                    const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                    std::vector<double>* HR_pointer_in,
                                    std::vector<T>* HK_pointer_in)
        : HR_pointer(HR_pointer_in), HK_pointer(HK_pointer_in), OperatorLCAO<T>(LM_in, kvec_d_in)
    {
        this->cal_type = lcao_sc_lambda;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:
    std::vector<double>* HR_pointer = nullptr;
    std::vector<T>* HK_pointer = nullptr;
};

} // namespace hamilt

#endif // SC_LAMBDA_LCAO_H