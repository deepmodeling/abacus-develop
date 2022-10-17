#ifndef DEEPKSLCAO_H
#define DEEPKSLCAO_H
#include "operator_lcao.h"
#include "module_base/timer.h"
#include "src_lcao/local_orbital_charge.h"

namespace hamilt
{

#ifndef __DEEPKSTEMPLATE
#define __DEEPKSTEMPLATE

template<class T> class DeePKS : public T 
{};

#endif

template<typename T>
class DeePKS<OperatorLCAO<T>> : public OperatorLCAO<T> 
{
    public:

    DeePKS<OperatorLCAO<T>>(
        Local_Orbital_Charge* loc_in,
        LCAO_Matrix* LM_in,
        std::vector<double>* HR_pointer_in,
        std::vector<T>* HK_pointer_in):loc(loc_in), HR_pointer(HR_pointer_in), HK_pointer(HK_pointer_in)
    {
        this->LM = LM_in;
        this->cal_type = lcao_deepks;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

    private:

    Local_Orbital_Charge* loc;

    std::vector<double>* HR_pointer = nullptr;

    std::vector<T>* HK_pointer = nullptr;

    bool HR_fixed_done = false;

};

}
#endif