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

template <typename T>
class OperatorEXX<OperatorLCAO<T>> : public OperatorLCAO<T>
{
  public:
    OperatorEXX<OperatorLCAO<T>>(LCAO_Matrix* LM_in,
                                 std::vector<double>* HR_pointer_in,
                                 std::vector<T>* HK_pointer_in,
                                 const K_Vectors& kv_in)
        : HR_pointer(HR_pointer_in), HK_pointer(HK_pointer_in), kv(kv_in), OperatorLCAO<T>(LM_in, kv_in.kvec_d)
    {
        this->cal_type = lcao_exx;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:
    std::vector<double>* HR_pointer = nullptr;

    std::vector<T>* HK_pointer = nullptr;

    bool HR_fixed_done = false;

    const K_Vectors& kv;
};

} // namespace hamilt
#endif