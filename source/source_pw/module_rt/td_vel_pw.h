#ifndef TD_VEL_PW_H
#define TD_VEL_PW_H

#include "source_estate/elecstate.h"
#include "source_pw/module_pwdft/op_pw_vel.h"
#include "source_pw/module_rt/td_obs_pw.h"

#include <complex>
#include <memory>
#include <vector>

namespace pw
{
/** @brief Charge and spin currents sharing one velocity action per k point. */
template <typename Real, typename Device>
class TDVelPW
{
  private:
    using Complex = std::complex<Real>;
    const UnitCell* cell_ = nullptr;
    const ModulePW::PW_Basis_K* basis_ = nullptr;
    const K_Vectors* kv_ = nullptr;
    pseudopot_cell_vnl* pp_ = nullptr;
    std::unique_ptr<hamilt::Velocity<Real, Device>> vel_;
    ct::Tensor vpsi_;
    ct::Tensor dots_;
    std::vector<Complex> diag_;
    std::vector<Complex> cross_;

  public:
    /** @brief Bind once; dependencies must outlive this object and retain their layout.
     * Ionic positions may change after the caller updates the structure factors.
     */
    void initialize(const UnitCell& cell, const ModulePW::PW_Basis_K& basis, const K_Vectors& kv, pseudopot_cell_vnl& pp);
    /** @brief Sample an endpoint on all ranks, restoring the wavefunction's k-point cursor.
     * Spin current requires npol=2; requests must agree across ranks.
     */
    TDVelResultPW calculate(const psi::Psi<Complex, Device>& psi,
                            const elecstate::ElecState& state,
                            const ModuleBase::Vector3<double>& a,
                            const CurrentRequest& request);
};
} // namespace pw
#endif
