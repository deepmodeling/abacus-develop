#ifndef ESOLVER_KS_PW_H
#define ESOLVER_KS_PW_H
#include "./esolver_ks.h"
#include "source_psi/setup_psi_pw.h" // mohan add 20251012
#include "source_pw/module_pwdft/vsep_pw.h"
#include "source_pw/module_pwdft/exx_helper_base.h"
#include "source_pw/module_pwdft/op_pw_vel.h"

#include <memory>
#include <source_base/macros.h>

namespace ModuleESolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class ESolver_KS_PW : public ESolver_KS
{
  private:
    using Real = typename GetTypeReal<T>::type;

  public:
    ESolver_KS_PW();

    ~ESolver_KS_PW();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

    double cal_energy() override;

    void cal_force(UnitCell& ucell, ModuleBase::matrix& force) override;

    void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) override;

    void after_all_runners(UnitCell& ucell) override;

    Exx_HelperBase* exx_helper = nullptr;

  protected:
    virtual void before_scf(UnitCell& ucell, const int istep) override;

    virtual void iter_init(UnitCell& ucell, const int istep, const int iter) override;

    virtual void iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver) override;

    virtual void after_scf(UnitCell& ucell, const int istep, const bool conv_esolver) override;

    virtual void others(UnitCell& ucell, const int istep) override;

    virtual void hamilt2rho_single(UnitCell& ucell, const int istep, const int iter, const double ethr) override;

    virtual void allocate_hamilt(const UnitCell& ucell);
    virtual void deallocate_hamilt();

    // Electronic wave function psi
    Setup_Psi_pw stp;

    // DFT-1/2 method
    VSep* vsep_cell = nullptr;

};
} // namespace ModuleESolver
#endif
