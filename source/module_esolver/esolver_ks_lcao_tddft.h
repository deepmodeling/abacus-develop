#ifndef ESOLVER_KS_LCAO_TDDFT_H
#define ESOLVER_KS_LCAO_TDDFT_H
#include "esolver_ks.h"
#include "esolver_ks_lcao.h"
#include "module_basis/module_ao/ORB_control.h"
#include "module_psi/psi.h"
#include "module_hamilt_lcao/hamilt_lcaodft/LCAO_hamilt.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_charge.h"
#include "module_hamilt_lcao/hamilt_lcaodft/local_orbital_wfc.h"
#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
#include "module_elecstate/elecstate_lcao_tddft.h"


namespace ModuleESolver
{

class ESolver_KS_LCAO_TDDFT : public ESolver_KS_LCAO
{
  public:
    ESolver_KS_LCAO_TDDFT();
    ~ESolver_KS_LCAO_TDDFT();
    void Init(Input& inp, UnitCell& cell) override;

    psi::Psi<std::complex<double>>* psi_laststep = nullptr;
    std::complex<double>** Hk_laststep = nullptr;
    std::complex<double>** Sk_laststep = nullptr;
    //same as pelec
    elecstate::ElecStateLCAO_TDDFT* pelec_td = nullptr;
    int td_htype = 1;

  protected:
    virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
    virtual void updatepot(const int istep, const int iter) override;
    virtual void afterscf(const int istep) override;
    void cal_edm_tddft();
};

} // namespace ModuleESolver
#endif
