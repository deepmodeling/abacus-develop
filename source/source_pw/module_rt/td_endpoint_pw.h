#ifndef TD_ENDPOINT_PW_H
#define TD_ENDPOINT_PW_H

#include "source_estate/elecstate_pw.h"
#include "source_pw/module_rt/td_orbital_ops.h"
#include "source_pw/module_pwdft/hamilt_pw.h"

namespace pw
{
/** @brief Clean initial orbitals and rebuild occupations and density if orthogonalization changes them. */
template <typename T, typename Device>
void orth_init_pw(hsolver::TDOrbOps<T, Device>& orbitals,
                  hsolver::HSOperator<T, Device>& op,
                  psi::Psi<T, Device>& current,
                  elecstate::ElecStatePW<T, Device>& estate,
                  Charge& charge,
                  UnitCell& cell,
                  const ModulePW::PW_Basis& rhod,
                  int nspin,
                  int nbands,
                  int iter,
                  bool do_orth);

/** @brief Rebuild an accepted endpoint from its orbitals while retaining all occupations. */
template <typename T, typename Device>
void rebuild_td_end(psi::Psi<T, Device>& current,
                    elecstate::ElecStatePW<T, Device>& estate,
                    Charge& charge,
                    UnitCell& cell,
                    hamilt::HamiltPW<T, Device>& hamilt,
                    hsolver::TDOrbOps<T, Device>& orbitals,
                    const ModulePW::PW_Basis_K& basis,
                    const ModulePW::PW_Basis& rhod,
                    const ModuleBase::Vector3<double>& vecpot,
                    int nspin,
                    int imp_sol,
                    bool sc_mag,
                    int dft_u,
                    const std::string& assume_iso);
} // namespace pw
#endif
