#ifndef ESOLVER_GETS_H
#define ESOLVER_GETS_H

#include "module_basis/module_nao/two_center_bundle.h"
#include "module_cell/unitcell.h"
#include "module_esolver/esolver_ks.h"
#include "module_hamilt_lcao/module_gint/gint_gamma.h"
#include "module_hamilt_lcao/module_gint/gint_k.h"

#include <memory>

namespace ModuleESolver
{

class ESolver_GetS : public ESolver_KS<std::complex<double>>
{
  public:
    ESolver_GetS();
    ~ESolver_GetS();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

    void after_all_runners(UnitCell& ucell){};

    void runner(UnitCell& ucell, const int istep) override;

    //! calculate total energy of a given system
    double cal_energy() {};

    //! calcualte forces for the atoms in the given cell
    void cal_force(UnitCell& ucell, ModuleBase::matrix& force) {};

    //! calcualte stress of given cell
    void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) {};

  protected:
    // we will get rid of this class soon, don't use it, mohan 2024-03-28
    Record_adj RA;

    // 2d block - cyclic distribution info
    Parallel_Orbitals pv;

    // used for k-dependent grid integration.
    Gint_k GK;

    // used for gamma only algorithms.
    Gint_Gamma GG;

    TwoCenterBundle two_center_bundle_;

    // temporary introduced during removing GlobalC::ORB
    LCAO_Orbitals orb_;
};
} // namespace ModuleESolver
#endif
