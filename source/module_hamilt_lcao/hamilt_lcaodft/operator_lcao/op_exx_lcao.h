#ifndef OPEXXLCAO_H
#define OPEXXLCAO_H
#include "module_base/timer.h"
#include "module_cell/klist.h"
#include "operator_lcao.h"
#ifdef __EXX
#include <RI/global/Tensor.h>
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
    using TAC = std::pair<int, std::array<int, 3>>;
public:
    OperatorEXX<OperatorLCAO<TK, TR>>(LCAO_Matrix* LM_in,
        hamilt::HContainer<TR>* hR_in,
        std::vector<TK>* hK_in,
        const K_Vectors& kv_in,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd_in = nullptr,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc_in = nullptr,
        int symbol_in = 0)
        : kv(kv_in), Hexxd(Hexxd_in), Hexxc(Hexxc_in), symbol(symbol_in) ,OperatorLCAO<TK, TR>(LM_in, kv_in.kvec_d, hR_in, hK_in)
    {
        this->cal_type = lcao_exx;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:

      bool HR_fixed_done = false;

      std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd = nullptr;
      std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc = nullptr;

      // test
      int symbol;

      const K_Vectors& kv;
};

} // namespace hamilt
#endif
#endif