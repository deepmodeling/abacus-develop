#ifndef ESOLVER_KS_PW_H
#define ESOLVER_KS_PW_H
#include "./esolver_ks.h"
#include "source_psi/setup_psi.h" // mohan add 20251012
#include "source_pw/module_pwdft/VSep_in_pw.h"
#include "source_pw/module_pwdft/global.h"
#include "source_pw/module_pwdft/module_exx_helper/exx_helper.h"
#include "source_pw/module_pwdft/operator_pw/velocity_pw.h"

#include <memory>
#include <source_base/macros.h>

namespace ModuleESolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class ESolver_KS_PW : public ESolver_KS<T, Device>
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

    Exx_Helper<T, Device> exx_helper;

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
    Setup_Psi<T, Device> stp;

    // DFT-1/2 method
    VSep* vsep_cell = nullptr;

    // for get_pchg and get_wf, use ctx as input of fft
    Device* ctx = {};

    // for device to host data transformation
    base_device::AbacusDevice_t device = {};

};
} // namespace ModuleESolver
#endif
