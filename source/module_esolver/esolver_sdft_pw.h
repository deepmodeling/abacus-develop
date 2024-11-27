#ifndef ESOLVER_SDFT_PW_H
#define ESOLVER_SDFT_PW_H

#include "esolver_ks_pw.h"
#include "module_hamilt_pw/hamilt_stodft/hamilt_sdft_pw.h"
#include "module_hamilt_pw/hamilt_stodft/sto_che.h"
#include "module_hamilt_pw/hamilt_stodft/sto_iter.h"
#include "module_hamilt_pw/hamilt_stodft/sto_wf.h"

namespace ModuleESolver
{

template <typename T, typename Device = base_device::DEVICE_CPU>
class ESolver_SDFT_PW : public ESolver_KS_PW<T, Device>
{
  private:
    using Real = typename GetTypeReal<T>::type;
  public:
    ESolver_SDFT_PW();
    ~ESolver_SDFT_PW();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

    double cal_energy() override;

    void cal_force(UnitCell& ucell, ModuleBase::matrix& force) override;

    void cal_stress(UnitCell& ucell, ModuleBase::matrix& stress) override;

  public:
    Stochastic_WF<T, Device> stowf;
    StoChe<Real, Device> stoche;
    hamilt::HamiltSdftPW<T, Device>* p_hamilt_sto = nullptr;

  protected:
    virtual void before_scf(UnitCell& ucell, const int istep) override;

    virtual void hamilt2density_single(UnitCell& ucell, const int istep, const int iter, const double ethr) override;

    virtual void others(UnitCell& ucell, const int istep) override;

    virtual void iter_finish(UnitCell& ucell, const int istep, int& iter) override;

    virtual void after_scf(UnitCell& ucell, const int istep) override;

    virtual void after_all_runners(UnitCell& ucell) override;

  private:
    int nche_sto;   ///< norder of Chebyshev
    int method_sto; ///< method of SDFT
};

} // namespace ModuleESolver
#endif
