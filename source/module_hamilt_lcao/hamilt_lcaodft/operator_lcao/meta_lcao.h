#ifndef METALCAO_H
#define METALCAO_H
#include "module_base/timer.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __METATEMPLATE
#define __METATEMPLATE

template <class T>
class Meta : public T
{
};

#endif

template <typename TK, typename TR>
class Meta<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    Meta<OperatorLCAO<TK, TR>>(
                          HS_Matrix_K<TK>* hsk_in,
                          const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                          HContainer<TR>* hR_in)
        : OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in)
    {
        this->cal_type = calculation_type::lcao_gint;
    }

    ~Meta<OperatorLCAO<TK, TR>>(){};

    virtual void contributeHR() override{}//do nothing now

    virtual void contributeHk(int ik) override{};//do nothing now

  private:
};

} // namespace hamilt
#endif