#ifndef OPEXXLCAO_H
#define OPEXXLCAO_H
#include "module_base/timer.h"
#include "module_cell/klist.h"
#include "operator_lcao.h"

namespace hamilt
{

#ifndef __OPEXXTEMPLATE
#define __OPEXXTEMPLATE

template <class T>
class OperatorEXX : public T
{
};

#endif

template <typename TK, typename TR>
class OperatorEXX<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    OperatorEXX<OperatorLCAO<TK, TR>>(LCAO_Matrix* LM_in,
                                 hamilt::HContainer<TR>* hR_in,
                                 std::vector<TK>* hK_in,
                                 const K_Vectors& kv_in)
        : kv(kv_in), OperatorLCAO<TK, TR>(LM_in, kv_in.kvec_d, hR_in, hK_in)
    {
        this->cal_type = lcao_exx;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:

    bool HR_fixed_done = false;

    const K_Vectors& kv;
};

} // namespace hamilt
#endif