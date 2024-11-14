#ifndef W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_ESOLVER_ESOLVER_KS_LCAO_TDDFT_H
#define W_ABACUS_DEVELOP_ABACUS_DEVELOP_SOURCE_MODULE_ESOLVER_ESOLVER_KS_LCAO_TDDFT_H
#include "esolver_ks.h"
#include "esolver_ks_lcao.h"
#include "module_elecstate/elecstate_lcao_tddft.h"
#include "module_hamilt_lcao/hamilt_lcaodft/record_adj.h"
#include "module_psi/psi.h"

namespace ModuleESolver
{

class ESolver_KS_LCAO_TDDFT : public ESolver_KS_LCAO<std::complex<double>, double>
{
  public:
    ESolver_KS_LCAO_TDDFT();

    ~ESolver_KS_LCAO_TDDFT();

    void before_all_runners(const Input_para& inp, UnitCell& cell) override;

    psi::Psi<std::complex<double>>* psi_laststep = nullptr;

    std::complex<double>** Hk_laststep = nullptr;

    std::complex<double>** Sk_laststep = nullptr;

    // same as pelec
    elecstate::ElecStateLCAO_TDDFT* pelec_td = nullptr;

    int td_htype = 1;

  protected:
    virtual void hamilt2density_single(const int istep, const int iter, const double ethr) override;

    virtual void update_pot(const int istep, const int iter) override;

    virtual void iter_finish(const int istep, int& iter) override;

    virtual void after_scf(const int istep) override;
};

} // namespace ModuleESolver
#endif
