#ifndef OPDFTULCAO_H
#define OPDFTULCAO_H

#include "source_lcao/module_operator_lcao/operator_lcao.h"
#include "source_pw/module_pwdft/dftu_base.h" // mohan add 20251107

namespace hamilt
{

template <class T>
class OperatorDFTU : public T
{
};

template <typename TK, typename TR>
class OperatorDFTU<OperatorLCAO<TK, TR>> : public OperatorLCAO<TK, TR>
{
  public:
    OperatorDFTU<OperatorLCAO<TK, TR>>(HS_Matrix_K<TK>* hsk_in,
                                  const std::vector<ModuleBase::Vector3<double>>& kvec_d_in,
                                  hamilt::HContainer<TR>* hR_in,
                                  const UnitCell& ucell_in,
                                  Plus_U_Base* dftu_in,
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

    Plus_U_Base *dftu;

    const UnitCell* ucell = nullptr;

    const std::vector<int>& isk;

    const int npol;
};
} // namespace hamilt
#endif
