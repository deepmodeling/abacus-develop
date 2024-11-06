#ifndef OPEXXLCAO_H
#define OPEXXLCAO_H

#ifdef __EXX

#include <RI/global/Tensor.h>
#include <RI/ri/Cell_Nearest.h>
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
enum Add_Hexx_Type { R, k };
template <typename TK, typename TR>
class OperatorEXX<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
    using TAC = std::pair<int, std::array<int, 3>>;
public:
    OperatorEXX<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
        hamilt::HContainer<TR>* hR_in,
        const K_Vectors& kv_in,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<double>>>>* Hexxd_in = nullptr,
        std::vector<std::map<int, std::map<TAC, RI::Tensor<std::complex<double>>>>>* Hexxc_in = nullptr,
        Add_Hexx_Type add_hexx_type_in = Add_Hexx_Type::R,
        const int istep_in = 0,
        int* two_level_step_in = nullptr,
        const bool restart_in = false);

    virtual void contributeHk(int ik) override;
    virtual void contributeHR() override;

  private:
      Add_Hexx_Type add_hexx_type = Add_Hexx_Type::R;
      int current_spin = 0;
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

      const int istep = 0; // the ion step

      void add_loaded_Hexx(const int ik);

      const K_Vectors& kv;

      // if k points has no shift, use cell_nearest to reduce the memory cost
      RI::Cell_Nearest<int, int, 3, double, 3> cell_nearest;
      bool use_cell_nearest = true;

      /// @brief Hexxk for all k-points, only for the 1st scf loop ofrestart load
      std::vector<std::vector<double>> Hexxd_k_load;
      std::vector<std::vector<std::complex<double>>> Hexxc_k_load;
};

} // namespace hamilt
#endif // __EXX
#include "op_exx_lcao.hpp"
#endif // OPEXXLCAO_H