#ifndef ESOLVER_GETS_H
#define ESOLVER_GETS_H

#include "module_basis/module_nao/two_center_bundle.h"
#include "module_cell/unitcell.h"
#include "module_esolver/esolver_ks.h"

#include <memory>

namespace ModuleESolver
{
template <typename TK, typename TR>
class ESolver_GetS : public ESolver_KS<TK>
{
  public:
    ESolver_GetS();
    ~ESolver_GetS();

    void before_all_runners(const Input_para& inp, UnitCell& ucell) override;

    void after_all_runners() override;

    void runner(const int istep, UnitCell& ucell) override;

  protected:
    // we will get rid of this class soon, don't use it, mohan 2024-03-28
    Record_adj RA;

    // 2d block - cyclic distribution info
    Parallel_Orbitals pv;

    TwoCenterBundle two_center_bundle_;

    //   // temporary introduced during removing GlobalC::ORB
    LCAO_Orbitals orb_;
};
} // namespace ModuleESolver
#endif
