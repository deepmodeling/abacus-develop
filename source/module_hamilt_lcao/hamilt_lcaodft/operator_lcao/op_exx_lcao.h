#ifndef OPEXXLCAO_H
#define OPEXXLCAO_H

#ifdef __EXX

#include <RI/global/Tensor.h>
#include "operator_lcao.h"
#include "module_cell/klist.h"

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
        int* two_level_step_in = nullptr,
        const bool restart_in = false);

    virtual void contributeHk(int ik) override;

  private:

      bool HR_fixed_done = false;

      std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd = nullptr;
      std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc = nullptr;

      /// @brief  the step of the outer loop.
      /// nullptr: no dependence on the number of two_level_step, contributeHk will do enerything normally.
      /// 0: the first outer loop. If restart, contributeHk will directly add Hexx to Hloc. else, do nothing.
      /// >0: not the first outer loop. contributeHk will do enerything normally.
      int* two_level_step = nullptr;
      /// @brief if restart, read and save Hexx, and directly use it during the first outer loop.
      bool restart = false;

      void add_loaded_Hexx(const int ik);
      const K_Vectors& kv;
};

} // namespace hamilt
#endif // __EXX
#include "op_exx_lcao.hpp"
#endif // OPEXXLCAO_H