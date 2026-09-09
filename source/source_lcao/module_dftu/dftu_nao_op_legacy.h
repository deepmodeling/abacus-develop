#ifndef OPDFTULCAO_H
#define OPDFTULCAO_H

#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_lcao/module_dftu/dftu_nao.h" // mohan add 20251107

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
                                  const UnitCell& ucell_in,
                                  Plus_U* dftu_in,
                                  const std::vector<int>& isk_in,
                                  const int npol_in)
        : isk(isk_in), npol(npol_in), OperatorLCAO<TK, TR>(hsk_in, kvec_d_in, hR_in)
    {
        this->cal_type = calculation_type::lcao_dftu;
        this->dftu = dftu_in;
        this->ucell = &ucell_in;
    }

    virtual void contributeHR() override;

    virtual void contributeHk(int ik) override;

  private:

    Plus_U *dftu;

    const UnitCell* ucell = nullptr;

    bool HR_fixed_done = false;

    const std::vector<int>& isk;

    const int npol;
};
} // namespace hamilt
#endif
