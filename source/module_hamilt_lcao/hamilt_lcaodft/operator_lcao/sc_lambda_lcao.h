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

template<typename TK, typename TR>
class OperatorScLambda<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    OperatorScLambda<OperatorLCAO<TK, TR>>(LCAO_Matrix* LM_in,
                                    const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                    hamilt::HContainer<TR>* hR_in,
                                    std::vector<TK>* hK_in,
                                    const std::vector<int>& isk_in)
        : isk(isk_in), OperatorLCAO<TK, TR>(LM_in, kvec_d_in, hR_in, hK_in)
    {
        this->cal_type = lcao_sc_lambda;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;
  private:

    const std::vector<int>& isk;
};

} // namespace hamilt

#endif // SC_LAMBDA_LCAO_H
