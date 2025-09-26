#ifndef ESOLVER_OF_TDDFT_H
#define ESOLVER_OF_TDDFT_H

#include "esolver_of.h"

namespace ModuleESolver
{
class ESolver_OF_TDDFT : public ESolver_OF
{
  public:
    ESolver_OF_TDDFT();
    ~ESolver_OF_TDDFT();

    virtual void runner(UnitCell& ucell, const int istep) override;

  protected:
    std::complex<double>** pphi_td = nullptr;                     // pphi[i] = ppsi.get_pointer(i), which will be freed in ~Psi().
    psi::Psi<std::complex<double>>* psi_td = nullptr;   
    //std::complex<double>** pdEdphi_phi_ = nullptr;             // (dE/dphi)*phi
    // ==================== main process of TDOFDFT ======================
     void before_opt(const int istep, UnitCell& ucell);
    void get_Hpsi(UnitCell& ucell, const std::complex<double>* const* psi_, ModulePW::PW_Basis* pw_rho, std::complex<double>** Hpsi);
    void get_tf_potential(const double* const* prho, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpot);
    void get_vw_potential_phi(const std::complex<double>* const* pphi, ModulePW::PW_Basis* pw_rho, std::complex<double>** Hpsi); // -1/2 \nabla^2 \phi
    void get_CD_potential(const std::complex<double>* const* psi_, ModulePW::PW_Basis* pw_rho, ModuleBase::matrix& rpot);
    void propagate_psi(UnitCell& ucell, std::complex<double>** pphi_, ModulePW::PW_Basis* pw_rho);
};
} // namespace ModuleESolver

#endif
