#ifndef ESOLVER_DM2RHO_H
#define ESOLVER_DM2RHO_H

#include "module_esolver/esolver_ks_lcao.h"

#include <memory>

namespace ModuleESolver
{

template <typename TK, typename TR>
class ESolver_DM2rho : public ESolver_KS_LCAO<TK, TR>
{
  public:
    ESolver_DM2rho();
    ~ESolver_DM2rho();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

    void after_all_runners(UnitCell& ucell) override;

    void runner(UnitCell& ucell, const int istep) override;
};
} // namespace ModuleESolver
#endif
