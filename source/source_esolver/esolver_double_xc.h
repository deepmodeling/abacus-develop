#ifndef ESOLVER_DOUBLE_XC_H
#define ESOLVER_DOUBLE_XC_H

#include "source_esolver/esolver_ks_lcao.h"

namespace ModuleESolver
{
// used in deepks, run target and base xc functional simultaneously
template <typename TK, typename TR>
class ESolver_DoubleXC : public ESolver_KS_LCAO<TK, TR>
{
  public:
    ESolver_DoubleXC();
    ~ESolver_DoubleXC();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

    void after_all_runners(UnitCell& ucell) override;

    void runner(UnitCell& ucell, const int istep) override;
};
} // namespace ModuleESolver
#endif
