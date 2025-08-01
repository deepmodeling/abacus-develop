#ifndef ESOLVER_DOUBLE_XC_H
#define ESOLVER_DOUBLE_XC_H

#include "source_esolver/esolver_ks_lcao.h"

namespace ModuleESolver
{
// used in deepks, run target and base xc functional simultaneously
template <typename TK, typename TR>
class ESolver_DoubleXC : public ESolver_KS_LCAO<TK, TR>
{
  public:
    ESolver_DoubleXC();
    ~ESolver_DoubleXC();

    void before_all_runners(UnitCell& ucell, const Input_para& inp) override;

    void cal_force(UnitCell& ucell, ModuleBase::matrix& force) override;

  protected:

    void before_scf(UnitCell& ucell, const int istep) override;

    void iter_finish(UnitCell& ucell, const int istep, int& iter, bool& conv_esolver) override;

    //! Hamiltonian
    hamilt::Hamilt<TK>* p_hamilt_base = nullptr;

    //! Electronic wavefunctions
    psi::Psi<TK>* psi_base = nullptr;

    //! Electronic states
    elecstate::ElecState* pelec_base = nullptr;

    //! Electorn charge density
    Charge chr_base;
};
} // namespace ModuleESolver
#endif
