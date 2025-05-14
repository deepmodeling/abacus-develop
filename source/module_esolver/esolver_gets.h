#ifndef ESOLVER_GETS_H
#define ESOLVER_GETS_H

#include "module_basis/module_nao/two_center_bundle.h"
#include "module_cell/unitcell.h"
#include "module_esolver/esolver_ks.h"

#include <memory>

namespace ModuleESolver
{

class ESolver_GetS : public ESolver_KS<std::complex<double>>
{
  public:
    ESolver_GetS();
    ~ESolver_GetS();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

    void after_all_runners(UnitCell& ucell) override;

    void runner(UnitCell& ucell, const int istep) override;

    //! calculate total energy of a given system
    double cal_energy() override;

    //! calcualte forces for the atoms in the given cell
    void cal_force(UnitCell& ucell, ModuleBase::matrix& force) override;

    //! calcualte stress of given cell
    void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) override;

  protected:
    // 2d block - cyclic distribution info
    Parallel_Orbitals pv;

    TwoCenterBundle two_center_bundle_;

    // temporary introduced during removing GlobalC::ORB
    LCAO_Orbitals orb_;
};
} // namespace ModuleESolver
#endif
