#ifndef ESOLVER_GETS_H
#define ESOLVER_GETS_H

#include "source_basis/module_nao/two_center_bundle.h"
#include "source_cell/unitcell.h"
#include "source_esolver/esolver_ks.h"

#include <memory>

namespace ModuleESolver
{

class ESolver_GetS : public ESolver_KS
{
  public:
    ESolver_GetS();
    ~ESolver_GetS();

    void before_all_runners(BaseCell& basecell, const Input_para& inp) override;

    void after_all_runners(BaseCell& basecell) override;

    void runner(BaseCell& basecell, const int istep) override;

    //! calculate total energy of a given system
    double cal_energy() override;

    //! calcualte forces for the atoms in the given cell
    void cal_force(BaseCell& basecell, ModuleBase::matrix& force) override;

    //! calcualte stress of given cell
    void cal_stress(BaseCell& basecell, ModuleBase::matrix& stress) override;

  protected:
    // 2d block - cyclic distribution info
    Parallel_Orbitals pv;

    TwoCenterBundle two_center_bundle_;

    // temporary introduced
    LCAO_Orbitals orb_;
};
} // namespace ModuleESolver
#endif
