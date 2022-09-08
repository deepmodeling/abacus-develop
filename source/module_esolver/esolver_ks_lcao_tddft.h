#ifndef ESOLVER_KS_LCAO_TDDFT_H
#define ESOLVER_KS_LCAO_TDDFT_H
#include "./esolver_ks.h"
#include "./esolver_ks_lcao.h"
#include "module_orbital/ORB_control.h"
#include "module_psi/psi.h"
#include "src_lcao/LCAO_hamilt.h"
#include "src_lcao/local_orbital_charge.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_lcao/record_adj.h"
#include "module_elecstate/elecstate_lcao_tddft.h"


namespace ModuleESolver
{

class ESolver_KS_LCAO_TDDFT : public ESolver_KS_LCAO
{
  public:
    ESolver_KS_LCAO_TDDFT();
    ~ESolver_KS_LCAO_TDDFT();
    void Init(Input& inp, UnitCell_pseudo& cell) override;

    psi::Psi<std::complex<double>>* psi_laststep = nullptr;
    elecstate::ElecStateLCAO_TDDFT* pelec_td = nullptr;

  protected:
    virtual void hamilt2density(const int istep, const int iter, const double ethr) override;
    virtual void eachiterinit(const int istep, const int iter) override;
    virtual void updatepot(const int istep, const int iter) override;
    virtual void afterscf() override;
    void cal_edm_tddft();
};

} // namespace ModuleESolver
#endif
