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

template <typename TK, typename TR>
class OperatorDFTU<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    OperatorDFTU<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                                  const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                  hamilt::HContainer<TR>* hR_in,
                                  const std::vector<int>& isk_in)
        : isk(isk_in), OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in)
    {
        this->cal_type = calculation_type::lcao_dftu;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:

    bool HR_fixed_done = false;

    const std::vector<int>& isk;
};
} // namespace hamilt
#endif