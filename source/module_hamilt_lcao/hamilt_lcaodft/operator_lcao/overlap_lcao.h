#ifndef OVERLAPLCAO_H
#define OVERLAPLCAO_H
#include "operator_lcao.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_gen_fixedH.h"

namespace hamilt
{

#ifndef __OVERLAPTEMPLATE
#define __OVERLAPTEMPLATE

template<class T> class Overlap : public T 
{};

#endif

template<typename T>
class Overlap<OperatorLCAO<T>> : public OperatorLCAO<T> 
{
    public:
      Overlap<OperatorLCAO<T>>(LCAO_gen_fixedH* genH_in,
                               LCAO_Matrix* LM_in,
                               const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                               std::vector<double>* SR_pointer_in,
                               std::vector<T>* SK_pointer_in)
          : genH(genH_in), SR_pointer(SR_pointer_in), SK_pointer(SK_pointer_in), OperatorLCAO<T>(kvec_d_in)
      {
          this->LM = LM_in;
          this->cal_type = lcao_fixed;
      }

    virtual void contributeHR() override;

    private:
    // use overlap matrix to generate fixed Hamiltonian
    LCAO_gen_fixedH* genH = nullptr;

    std::vector<double>* SR_pointer = nullptr;

    std::vector<T>* SK_pointer = nullptr;

    bool SR_fixed_done = false;

};

}
#endif