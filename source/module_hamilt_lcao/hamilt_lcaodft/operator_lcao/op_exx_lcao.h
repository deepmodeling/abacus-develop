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
        int* two_level_step = nullptr,
        const bool restart = false)
    : kv(kv_in), Hexxd(Hexxd_in), Hexxc(Hexxc_in), OperatorLCAO<TK, TR>(LM_in, kv_in.kvec_d, hR_in, hK_in),
        two_level_step(two_level_step), restart(restart)
    {
        this->cal_type = lcao_exx;
        if (restart)
        {///  Now only Hexx depends on DM, so we can directly read Hexx to reduce the computational cost.
        /// If other operators depends on DM, we can also read DM and then calculate the operators to save the memory to store operator terms.
            assert(two_level_step != nullptr);
            /// read in Hexx
            if (std::is_same<TK, double>::value)
            {
                this->LM->Hexxd_k_load.resize(this->kv.nks);
                for (int ik = 0; ik < this->kv.nks; ik++)
                {
                    this->LM->Hexxd_k_load[ik].resize(this->LM->ParaV->get_local_size(), 0.0);
                    restart = GlobalC::restart.load_disk("Hexx", ik,
                        this->LM->ParaV->get_local_size(), this->LM->Hexxd_k_load[ik].data(), false);
                    if (!restart) break;
                }
            }
            else
            {
                this->LM->Hexxc_k_load.resize(this->kv.nks);
                for (int ik = 0; ik < this->kv.nks; ik++)
                {
                    this->LM->Hexxc_k_load[ik].resize(this->LM->ParaV->get_local_size(), 0.0);
                    restart = GlobalC::restart.load_disk("Hexx", ik,
                        this->LM->ParaV->get_local_size(), this->LM->Hexxc_k_load[ik].data(), false);
                    if (!restart) break;
                }
            }
            if (!restart) std::cout << "WARNING: Hexx not found, restart from the non-exx loop. \n \
                If the loaded charge density is EXX-solved, this may lead to poor convergence." << std::endl;
            GlobalC::restart.info_load.load_H_finish = restart;
        }
    }

    virtual void contributeHR() override;

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
#endif
#endif